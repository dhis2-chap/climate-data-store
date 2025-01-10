import json
from typing import List
import cdsapi
import numpy as np
import pandas as pd
from pydantic import BaseModel
import xarray as xr
from shapely.geometry import shape
from forecast_sources import ECMWF, uk_met_office
from seasonal import SeasonalForecastHandler, SeasonalForecastHandlerConfig
import sys
import os

file_name_base = "seasonal-forecast"

class FetchCopernicusDataConfig(BaseModel):

    class Config:
        arbitrary_types_allowed = True

    originating_centre: str = "ukmo" or "ecmwf" or "meteo_france" or "dwd" or "cmcc" or "ncep" or "jma" or "eccc"
    features : List[object]
    indicator : str = "2m_temperature" or "total_precipitation"
    file_name_postfix : str = ""
    periode_type : str = "M" or "W-MON" or "D" or "W-SUN"
    skip_download : bool = False
    forecast_issued : np.datetime64

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

class BoundingBox(BaseModel):
    north: float
    south: float
    east: float
    west: float

class FetchCopernicusData():

    def __init__(self, config : FetchCopernicusDataConfig):
        self.originating_centre = config.originating_centre
        self.features = config.features
        self.grib_file_name = f"grib/{file_name_base}{config.file_name_postfix}.grib"
        self.netcdf_file_name = f"{file_name_base}{config.file_name_postfix}.nc"
        self.periode_type = config.periode_type
        self.indicator = config.indicator
        self.skip_download = config.skip_download
        self.forecast_issued = config.forecast_issued

    def get_data(self):
        
        request_config = self._get_request_config_for_originating_centre(self.indicator)

        self.fetch_data(request_config, is_value_type_sum=(self.indicator == "total_precipitation"), skip_download=self.skip_download)
        
        self._calculate_per_periode_and_time(indicator_dict[self.indicator])

    def _convert_from_grib_to_netcdf(self):
        ds = xr.open_dataset(self.grib_file_name, engine="cfgrib")
        print("Converting GRIB to netcdf..")
        ds.to_netcdf(self.netcdf_file_name)

    def _validate_config(self, config):
        if len(config) == 0:
            raise Exception("No configuration found for the given originating centre and variable")
        
        if len(config) > 1:
            raise Exception("Provied combination or originating centre and variable returned more than one result. Check forecast_sources.json file for duplicates.")

        return config[0]
    
    def _get_request_config_for_originating_centre(self, variable: str):
        with open("forecast_sources.json") as file:
            sources = json.load(file)
        
        config = [source for source in sources if source['originating_centre'] == self.originating_centre and source['variable'][0] == variable]
        return self._validate_config(config)
       

    def _calculate_per_periode_and_time(self, variable):

        config = SeasonalForecastHandlerConfig(
            netcdf_file_name=self.netcdf_file_name,
            variable=variable,
            features=self.features,
            periode_type=self.periode_type,
            output_file_postfix=self.originating_centre,
            measurement_unit=measurement_values[self.indicator],
            total_sum_value=is_total_sum_value[self.indicator]
        )

        sfh = SeasonalForecastHandler(config=config)
        sfh.calucate()

    #This only support precipitation for monthly data for now
    def _get_leadtime_hours_for_sum_indicator(self, dataset_starting_date : np.datetime64, maximum_lead_time_hours : int):

        lead_time_out_of_range = False
        temp_date = dataset_starting_date
        lead_time_hours = []

        while(not lead_time_out_of_range):
            current_date = np.datetime64(temp_date, 'M')
            start_date_next_periode = current_date + np.timedelta64(1, 'M')

            next_date = np.datetime64(f'{start_date_next_periode.item().year}-{start_date_next_periode.item().month:02d}-01')

            #calculate the number of days to the next period, minus one, since we want last day in previous period
            number_of_days_to_next_periode = next_date.astype('int') - dataset_starting_date.astype('int')# - 1

            print(number_of_days_to_next_periode)

            lead_time_houre = 24 * int(number_of_days_to_next_periode)

            if(lead_time_houre > maximum_lead_time_hours):
                lead_time_out_of_range = True
            else:
                lead_time_hours.append(str(lead_time_houre))
                temp_date = next_date

        return lead_time_hours
    
    def create_request_body(self, request_config, bounding_box : BoundingBox, request_dataset_issued : np.datetime64):
        return {
            "originating_centre": request_config["originating_centre"],
            "data_format": request_config["data_format"], 
            "variable": request_config["variable"],
            "system": str(request_config["system"]),
            "year": [str(request_dataset_issued.item().year)],
            "month": [str(request_dataset_issued.item().month)] if (request_dataset_issued.item().month > 9) else ["0"+str(request_dataset_issued.item().month)],
            "day": ["01"],#[str(request_dataset_issued.item().day)],
            "leadtime_hour": request_config["leadtime_hour"],
            "area":  [bounding_box.north, bounding_box.west, bounding_box.south, bounding_box.east],
        }

    def _get_dataset_issued_date(self, today : np.datetime64) -> np.datetime64:

        if today.item().day > 11:
            return np.datetime64(f'{today.item().year}-{today.item().month:02d}-01', "D")
        else:
            current_month = np.datetime64(today, 'M')
            previous_month = current_month - np.timedelta64(1, 'M')

            return np.datetime64(f'{previous_month.item().year}-{previous_month.item().month:02d}-01', "D")

    def fetch_data(self, request_config, is_value_type_sum=False, skip_download=False):
        copernicus_client = cdsapi.Client()

        #add the bouding box to the request
        bounding_box = self._getBoundingBox(self.features)

        print("bounding box: ",  bounding_box.model_dump())

        
        request_dataset_issued : np.datetime64 = self._get_dataset_issued_date(self.forecast_issued)


        if(is_value_type_sum):
            request_config["leadtime_hour"] = self._get_leadtime_hours_for_sum_indicator(request_dataset_issued, int(request_config["max_leadtime_hour"]))


        request_body = self.create_request_body(request_config, bounding_box, request_dataset_issued)
        print(request_body)
        
        if not skip_download:      
            copernicus_client.retrieve('seasonal-original-single-levels', request_body, f'{self.grib_file_name}')
            self._convert_from_grib_to_netcdf()

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
              • (optional) skipDownload:    skip download of netCDF-file from Copernicus, useful when you already have the file downloaded. (True/False)\n
              • (optional) date:            date of forecast issued, format 'YYYY-MM-DD' default is today\n
              """)
        sys.exit(1)

    file_path = sys.argv[1]
    indicator = sys.argv[2]

    try:
        skip_download = sys.argv[3]
    except (IndexError):
        skip_download = False

    try:
        forecast_issued = np.datetime64(sys.argv[4], 'D')
    except (IndexError):
        forecast_issued = np.datetime64('today', 'D')

    with open(file_path) as file:
        geojson_data = json.load(file)

    features = geojson_data['features']


    #find filename
    file_name_geojson = os.path.splitext(os.path.basename(file_path))[0]

    # Create directories if they do not exist
    if not os.path.exists("grib"):
        os.makedirs("grib")

    if not os.path.exists("results"):
        os.makedirs("results")

    config = FetchCopernicusDataConfig(
        originating_centre="ecmwf",
        features=features,
        file_name_postfix="-"+file_name_geojson,
        periode_type="M",
        indicator=indicator,
        skip_download=skip_download,
        forecast_issued=forecast_issued
    )

    fetch_data = FetchCopernicusData(config)
    fetch_data.get_data()



