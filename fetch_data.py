import json
from typing import List
import cdsapi
import numpy as np
import pandas as pd
from pydantic import BaseModel, model_validator
import xarray as xr
from shapely.geometry import shape
from datetime import datetime, timedelta
import logging
import sys
import os
import json
from pathlib import Path

from seasonal import SeasonalForecastHandler, SeasonalForecastHandlerConfig
from utils import generate_hash, increment_months

logging.basicConfig(level=logging.DEBUG)

SCRIPT_DIR = Path(__file__).parent
DEFAULT_OUTPUT_FOLDER = SCRIPT_DIR

class FetchCopernicusDataConfig(BaseModel):

    class Config:
        arbitrary_types_allowed = True

    originating_centre: str = "ukmo" or "ecmwf" or "meteo_france" or "dwd" or "cmcc" or "ncep" or "jma" or "eccc"
    features : List[object]
    indicator : str = "2m_temperature" or "total_precipitation"
    output_folder : str = DEFAULT_OUTPUT_FOLDER
    file_name_postfix : str = ""  # NOTE: applies to result file only
    years : List[int]
    forecast_length : int = 3 * 31 * 24  # default is 3x 31-day months in hours

    @model_validator(mode='before')
    def coerce_input_types(cls, data):
        # default to current year
        if not data.get('years', None):
            data['years'] = [datetime.today().year]

        # ensure list of years
        if not isinstance(data['years'], list):
            data['years'] = [data['years']]

        return data

indicator_dict = {
    "2m_temperature" : "t2m",
    "total_precipitation" : "tp"
}

leadtime_intervals = {
    "2m_temperature" : 6,
    "total_precipitation" : 24,
}

measurement_values =  {
    "2m_temperature" : "K",
    "total_precipitation" : "m"
}

is_total_sum_value =  {
    "2m_temperature" : False,
    "total_precipitation" : True
}

class BoundingBox(BaseModel):
    north: float
    south: float
    east: float
    west: float

class FetchCopernicusData():

    def __init__(self, config : FetchCopernicusDataConfig):
        self.originating_centre = config.originating_centre
        self.features = config.features
        self.output_folder = config.output_folder
        self.file_name_postfix = config.file_name_postfix # not used for now... 
        self.indicator = config.indicator
        self.years = config.years
        self.forecast_length = config.forecast_length

    def check_output_folders(self):
        Path(f'{self.output_folder}/grib').mkdir(parents=True, exist_ok=True)
        Path(f'{self.output_folder}/netcdf').mkdir(parents=True, exist_ok=True)
        Path(f'{self.output_folder}/results').mkdir(parents=True, exist_ok=True)

    def get_data(self):
        # init
        self.check_output_folders()
        request_config = self._get_request_config_for_originating_centre(self.indicator)

        # download data
        self.fetch_data(request_config)

    def _convert_from_grib_to_netcdf(self):
        ds = xr.open_dataset(f'{self.output_folder}/{self.grib_file_name}', engine="cfgrib")
        print("Converting GRIB to netcdf..")
        ds.to_netcdf(f'{self.output_folder}/{self.netcdf_file_name}')

    def _validate_config(self, config):
        if len(config) == 0:
            raise Exception("No configuration found for the given originating centre and variable")
        
        if len(config) > 1:
            raise Exception("Provided combination or originating centre and variable returned more than one result. Check forecast_sources.json file for duplicates.")

        return config[0]
    
    def _get_request_config_for_originating_centre(self, variable: str):
        with open(f"{SCRIPT_DIR}/forecast_sources.json") as file:
            sources = json.load(file)
        
        config = [source for source in sources if source['originating_centre'] == self.originating_centre and source['variable'][0] == variable]
        return self._validate_config(config)

    def get_forecast_handler(self, forecast_date, period_type : str, period_count : int = None):
        if not isinstance(forecast_date, datetime):
            forecast_date = datetime.fromisoformat(forecast_date)

        config = SeasonalForecastHandlerConfig(
            netcdf_file=f'{self.output_folder}/{self.netcdf_file_name}',
            variable=indicator_dict[self.indicator],
            features=self.features,
            forecast_date=forecast_date,
            period_type=period_type,
            period_count=period_count,
            measurement_unit=measurement_values[self.indicator],
            total_sum_value=is_total_sum_value[self.indicator]
        )

        sfh = SeasonalForecastHandler(config=config)
        return sfh

    def _get_all_leadtime_hours(self, leadtime_interval : int, max_leadtime : int):
        '''
        Getting all leadtime hours given a leadtime interval, max forecast length, and period type.
        '''
        lead_time_hours = []
        lead_time_hour = 0

        while lead_time_hour < max_leadtime:
            lead_time_hour += leadtime_interval
            lead_time_hours.append(str(lead_time_hour))

        return lead_time_hours
    
    def create_request_body(self, request_config, bounding_box : BoundingBox, request_years : List[int]):
        today = datetime.today()
        if len(request_years) == 1 and request_years[0] == today.year:
            # only requesting current year, limit the nr of months
            request_months = list(range(1, today.month + 1))
        else:
            # requesting historical years, all months required
            # TODO: this wont work if combining with latest year since it will include nonexistant months
            request_months = list(range(1, 12 + 1))

        return {
            "originating_centre": request_config["originating_centre"],
            "data_format": request_config["data_format"], 
            "variable": request_config["variable"],
            "system": str(request_config["system"]),
            "year": [str(yr) for yr in request_years],
            "month": [str(mn).zfill(2) for mn in request_months],
            "day": ["01"],
            "leadtime_hour": request_config["leadtime_hour"],
            "area":  [bounding_box.north, bounding_box.west, bounding_box.south, bounding_box.east],
        }

    # def _get_dataset_issued_date(self, issued_date : datetime) -> datetime:
    #     today = datetime.today()

    #     if issued_date.year == today.year & issued_date.month == today.month:
    #         # requesting latest available dataset (same as current month)
    #         if today.day > 11:
    #             # we have passed the 11th of the current month, which means the current month dataset should be available 
    #             return datetime(today.year, today.month, 1)
    #         else:
    #             # current month's dataset isn't available until after 11th, revert to the previous month's dataset instead
    #             previous_month = increment_months(today, -1)
    #             return datetime(previous_month.year, previous_month.month, 1)
    #     else:
    #         # requesting dataset from a historical month
    #         # datasets are always issued for the 1st of each month
    #         return datetime(issued_date.year, issued_date.month, 1)

    def fetch_data(self, request_config):
        copernicus_client = cdsapi.Client(timeout=300, quiet=False)

        # add the bounding box to the request
        bounding_box = self._getBoundingBox(self.features)
        print("bounding box: ",  bounding_box.model_dump())

        leadtime_interval = leadtime_intervals[self.indicator]
        request_config['leadtime_hour'] = self._get_all_leadtime_hours(leadtime_interval, self.forecast_length)

        # set api request params
        request_body = self.create_request_body(request_config, bounding_box, self.years)
        print(request_body)

        # determine file names based on request input
        request_hash = generate_hash(request_body)
        self.file_name_base = f'request_hash_{request_hash}'
        self.grib_file_name = f"grib/{self.file_name_base}.grib"
        self.netcdf_file_name = f"netcdf/{self.file_name_base}.nc"

        # only fetch data if not previously downloaded
        if not os.path.exists(f'{self.output_folder}/{self.netcdf_file_name}'):
            print('downloading data')
            copernicus_client.retrieve('seasonal-original-single-levels', request_body, f'{self.output_folder}/{self.grib_file_name}')
            self._convert_from_grib_to_netcdf()
        
        else:
            print('data download already exists, skipping')

    def _getBoundingBox(self, features) -> BoundingBox:
        north = -90
        west = 180
        south = 90
        east = -180

        # Go through each feature in the GeoJSON file
        for feature in features:
            # Convert the feature's geometry to a Shapely geometry
            geom = shape(feature['geometry'])
            
            # This returns (minx, miny, maxx, maxy)
            bounds = geom.bounds 
            
            lat_north = bounds[3] #maxY
            lat_south = bounds[1] #minY
            lon_east = bounds[2] #maxX
            lon_west = bounds[0] #minX

            if lat_north > north:
                north = lat_north

            if lat_south < south:
                south = lat_south

            if lon_east > east:
                east = lon_east

            if lon_west < west:
                west = lon_west
        
        return BoundingBox(north=north, south=south, east=east, west=west)


    
if __name__ == "__main__":

    
    if len(sys.argv) < 3:
        print("""
              Usage: python fetch_data.py [file_path] [indicator] [skipDownload] \n\n 
              • file_path:                  path to geojson-file\n 
              • indicator:                  '2m_temperature' or 'total_precipitation'\n 
              • (optional) years:           year for which forecasts will be fetched, one or more separated by comma, format 'YYYY', default is current year\n
              • (optional) forecastLength:  how many hours into the future to fetch forecast for, default is 3x 31-day months\n

              example: python fetch_data.py data/orgUnitsSingleSierra.geojson total_precipitation
              """)
        sys.exit(1)

    # get cmd args

    file_path = sys.argv[1]
    indicator = sys.argv[2]

    try:
        years = sys.argv[3]
        years = [yr.strip() for yr in years.split(',')]
    except (IndexError):
        years = None

    try:
        forecast_length = sys.argv[4]
    except (IndexError):
        forecast_length = None

    # load geojson features
    with open(file_path) as file:
        geojson_data = json.load(file)
    features = geojson_data['features']

    # find filename
    file_name_geojson = os.path.splitext(os.path.basename(file_path))[0]

    # get output folder as cmd working dir
    output_dir = os.pwd()

    # create config
    config = FetchCopernicusDataConfig(
        originating_centre="ecmwf",
        features=features,
        file_name_postfix="-"+file_name_geojson,
        indicator=indicator,
        years=years,
        forecast_length=forecast_length,
    )

    # get data
    fetch_data = FetchCopernicusData(config)
    fetch_data.get_data()



