import json
from typing import List
import cdsapi
import numpy as np
import pandas as pd
from pydantic import BaseModel, model_validator
import xarray as xr
from shapely.geometry import shape
from forecast_sources import ECMWF, uk_met_office
from seasonal import SeasonalForecastHandler, SeasonalForecastHandlerConfig
from datetime import datetime, date, timedelta
from calendar import monthrange
import logging
import sys
import os
import json
import hashlib
from pathlib import Path

logging.basicConfig(level=logging.DEBUG)

DEFAULT_OUTPUT_FOLDER = Path(__file__).parent

def generate_file_name_base(obj):
    obj_string = json.dumps(obj).encode('utf8')
    obj_hash = hashlib.md5(obj_string).hexdigest()
    print('generating file name base from input args:', obj_string, '->', obj_hash)
    return obj_hash

class FetchCopernicusDataConfig(BaseModel):

    class Config:
        arbitrary_types_allowed = True

    originating_centre: str = "ukmo" or "ecmwf" or "meteo_france" or "dwd" or "cmcc" or "ncep" or "jma" or "eccc"
    features : List[object]
    indicator : str = "2m_temperature" or "total_precipitation"
    output_folder : str = DEFAULT_OUTPUT_FOLDER
    file_name_postfix : str = ""  # NOTE: applies to result file only
    period_type : str = "M" or "W-MON" or "D" or "W-SUN"
    forecast_issued : datetime
    forecast_length : int

    @model_validator(mode='before')
    def coerce_input_types(cls, data):
        # default datetime to today
        if not data.get('forecast_issued', None):
            data['forecast_issued'] = datetime.today()

        # create datetime from str
        elif isinstance(data['forecast_issued'], str):
            data['forecast_issued'] = datetime.fromisoformat(data['forecast_issued'])

        # default forecast length depending on period type
        if not data.get('forecast_length', None):
            data['forecast_length'] = {'M':3, 'W':8, 'D':14}[data['period_type'][0]]

        return data

indicator_dict = {
    "2m_temperature" : "t2m",
    "total_precipitation" : "tp"
}

measurement_values =  {
    "2m_temperature" : "K",
    "total_precipitation" : "m"
}

is_total_sum_value =  {
    "2m_temperature" : False,
    "total_precipitation" : True
}

def increment_months(start_date, months):
    # Compute the total number of months since year 0
    total_months = start_date.year * 12 + (start_date.month - 1) + months
    new_year = total_months // 12
    new_month = total_months % 12 + 1

    # Handle end-of-month cases by adjusting the day if necessary
    last_day_of_new_month = monthrange(new_year, new_month)[1]
    new_day = min(start_date.day, last_day_of_new_month)

    return datetime(new_year, new_month, new_day)

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
        self.file_name_postfix = config.file_name_postfix
        self.period_type = config.period_type
        self.indicator = config.indicator
        self.forecast_issued = config.forecast_issued
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
        
        # calculate results
        df = self._calculate_per_period_and_time(indicator_dict[self.indicator])

        # save
        self._save_calculated_results(df)

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
        with open("forecast_sources.json") as file:
            sources = json.load(file)
        
        config = [source for source in sources if source['originating_centre'] == self.originating_centre and source['variable'][0] == variable]
        return self._validate_config(config)

    def _calculate_per_period_and_time(self, variable):

        config = SeasonalForecastHandlerConfig(
            netcdf_file=f'{self.output_folder}/{self.netcdf_file_name}',
            variable=variable,
            features=self.features,
            period_type=self.period_type,
            measurement_unit=measurement_values[self.indicator],
            total_sum_value=is_total_sum_value[self.indicator]
        )

        sfh = SeasonalForecastHandler(config=config)
        df = sfh.calculate()
        return df
    
    def _save_calculated_results(self, df):
        df.to_csv(
            f"{self.output_folder}/{self.results_file_name}",  
            sep=";",
            index=False
        )

    def _get_snapshot_leadtime_hours(self, dataset_starting_date : datetime, forecast_length : int, period_type : str, leadtime_interval : int):
        '''
        Getting all necessary leadtime hours when the forecast represents snapshot valuse (eg temperature).
        All available leadtime hours are required since we need to calculate some aggregate statistics for a given period. 
        '''
        lead_time_hours = []
        lead_time_hour = 0

        if period_type == 'M':
            dataset_ending_date = increment_months(dataset_starting_date, forecast_length)
        elif period_type[0] == 'W':
            dataset_ending_date = dataset_starting_date + timedelta(weeks=forecast_length)
        elif period_type == 'D':
            dataset_ending_date = dataset_starting_date + timedelta(days=forecast_length)
        
        forecast_length_hours = (dataset_ending_date - dataset_starting_date).days * 24

        while lead_time_hour < forecast_length_hours:
            lead_time_hour += leadtime_interval
            lead_time_hours.append(str(lead_time_hour))

        return lead_time_hours

    def _get_cumulative_leadtime_hours(self, dataset_starting_date : datetime, forecast_length : int, period_type : str):
        '''
        Getting only the necessary leadtime hours when the forecast represents cumulative valuse (eg precipitaiton).
        Returns list of hours since starting date into the future to forecast, at intervals specified by period type.
        '''
        lead_time_hours = []
        current_date = dataset_starting_date

        while len(lead_time_hours) < forecast_length:
            if period_type == 'M':
                next_date = increment_months(current_date, 1)
            elif period_type[0] == 'W':
                next_date = current_date + timedelta(weeks=1)
            elif period_type == 'D':
                next_date = current_date + timedelta(days=1)

            number_of_days_since_starting_date = (next_date - dataset_starting_date).days

            print(number_of_days_since_starting_date)

            lead_time_hour = 24 * int(number_of_days_since_starting_date)

            lead_time_hours.append(str(lead_time_hour))
            current_date = next_date

        return lead_time_hours
    
    def create_request_body(self, request_config, bounding_box : BoundingBox, request_dataset_issued : datetime):
        return {
            "originating_centre": request_config["originating_centre"],
            "data_format": request_config["data_format"], 
            "variable": request_config["variable"],
            "system": str(request_config["system"]),
            "year": [str(request_dataset_issued.year)],
            "month": [str(request_dataset_issued.month).zfill(2)],
            "day": ["01"],
            "leadtime_hour": request_config["leadtime_hour"],
            "area":  [bounding_box.north, bounding_box.west, bounding_box.south, bounding_box.east],
        }

    def _get_dataset_issued_date(self, issued_date : datetime) -> datetime:
        today = datetime.today()

        if issued_date.year == today.year & issued_date.month == today.month:
            # requesting latest available dataset (same as current month)
            if today.day > 11:
                # we have passed the 11th of the current month, which means the current month dataset should be available 
                return datetime(today.year, today.month, 1)
            else:
                # current month's dataset isn't available until after 11th, revert to the previous month's dataset instead
                previous_month = increment_months(today, -1)
                return datetime(previous_month.year, previous_month.month, 1)
        else:
            # requesting dataset from a historical month
            # datasets are always issued for the 1st of each month
            return datetime(issued_date.year, issued_date.month, 1)

    def fetch_data(self, request_config):
        copernicus_client = cdsapi.Client(timeout=300, quiet=False)

        # add the bounding box to the request
        bounding_box = self._getBoundingBox(self.features)
        print("bounding box: ",  bounding_box.model_dump())

        request_dataset_issued : datetime = self._get_dataset_issued_date(self.forecast_issued)

        if is_total_sum_value[self.indicator]:
            # cumulative forecast, only requires leadtime hours between each period type
            request_config['leadtime_hour'] = self._get_cumulative_leadtime_hours(request_dataset_issued, self.forecast_length, self.period_type)
        else:
            # snapshot forecast, requires all available leadtime hours to calculate aggregate stats
            leadtime_interval = 6 # hardcoded to 2m temperature for now
            request_config['leadtime_hour'] = self._get_snapshot_leadtime_hours(request_dataset_issued, self.forecast_length, self.period_type, leadtime_interval)

        # set api request params
        request_body = self.create_request_body(request_config, bounding_box, request_dataset_issued)
        print(request_body)

        # determine file names based on input
        self.file_name_base = generate_file_name_base(request_body)
        self.grib_file_name = f"grib/{self.file_name_base}.grib"
        self.netcdf_file_name = f"netcdf/{self.file_name_base}.nc"
        self.results_file_name = f"results/{self.file_name_base}{self.file_name_postfix}.csv"

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
              • periodType:                 'M' or 'W-MON' or 'D'\n 
              • (optional) date:            date of forecast issued, format 'YYYY-MM-DD' default is today\n
              • (optional) forecastLength:  how far into the future to fetch forecast for, default is 3 months, 8 weeks, or 14 days

              example: python fetch_data.py data/orgUnitsSingleSierra.geojson total_precipitation M
              """)
        sys.exit(1)

    # get cmd args

    file_path = sys.argv[1]
    indicator = sys.argv[2]
    period_type = sys.argv[3]

    try:
        forecast_issued = sys.argv[4]
    except (IndexError):
        forecast_issued = None

    try:
        forecast_length = sys.argv[5]
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
        period_type="M",
        indicator=indicator,
        forecast_issued=forecast_issued,
        forecast_length=forecast_length,
    )

    # get data
    fetch_data = FetchCopernicusData(config)
    fetch_data.get_data()



