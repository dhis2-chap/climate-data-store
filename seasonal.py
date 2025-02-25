#import cdsapi
from datetime import datetime
import json
from typing import Dict, List
from xmlrpc.client import DateTime
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd
import xarray as xr
import cdsapi
from pyproj import Geod
from shapely.geometry import shape
import rioxarray
import geojson
import shapely
from pydantic import BaseModel, model_validator
from pydantic.json_schema import JsonSchemaValue
from shapely import LineString, MultiPoint, Polygon
import matplotlib.pyplot as plt
from dateutil.relativedelta import relativedelta
from datetime import timedelta

from utils import increment_months

class PeriodValue(BaseModel):
    nameOfPeriod : str
    value : float

class PointValue(BaseModel):

    class Config:
        arbitrary_types_allowed = True

    date : np.datetime64
    value : object
    org_unit_id : str
    org_unit_name : str = ""

class SeasonalForecastHandlerConfig(BaseModel):
    variable : str
    netcdf_file : str
    features : List[object]
    forecast_date : datetime
    period_type : str = "M" or "W-MON" or "D" or "W-SUN"
    period_count : int  # number of periods to forecast as defined by period_type
    aggregation_method : str = "mean" or "sum" or "max" or "min"
    measurement_unit : str = "kelvin" or "m"
    total_sum_value : bool = False

    @model_validator(mode='before')
    def coerce_input_types(cls, data):
        # default period count depending on period type
        if not data.get('period_count', None):
            data['period_count'] = {'D':14, 'W':8, 'M':3}[data['period_type']]

        return data

converters = {
    "K" : lambda x: x - 273.15, #converts kelvin to celsius
    "m" : lambda x: x * 1000 #converts meter to millimeter
}

class SeasonalForecastHandler():

    class Config:
        arbitrary_types_allowed=True

    def __init__(self, config : SeasonalForecastHandlerConfig):
        self.variable = config.variable
        self.netcdf_file = config.netcdf_file
        self.features = config.features
        self.forecast_date = config.forecast_date
        self.period_type = config.period_type
        self.period_count = config.period_count
        self.aggregation_method = config.aggregation_method
        self.measurement_unit = config.measurement_unit
        self.total_sum_value = config.total_sum_value

    def kelvin_to_celsius(self, value):
        return value - 273.15
    
    def m_to_mm(self, value):
        return value * 1000

    def _group_by_month(self, df : List[PointValue]):
        
        df['date'] = pd.to_datetime(df['date'])

        # Set 'time' as the index of the DataFrame
        df.set_index('date', inplace=True)

        # Group by month and calculate mean of all values
        monthly_stats = df.resample(self.period_type).agg({'value': ["mean"]})

        return monthly_stats
    
    def _find_center_of_coordinates(self, geometry : List[List[int]]):
        p : Polygon = shape(geometry)
        center = shapely.centroid(p)
        
        return center

    def _get_idx_for_nearest_point(self, ds, center):
        long_values = ds.longitude.values
        lat_values = ds.latitude.values
        long_idx = abs(long_values - center.x).argmin()
        lat_idx = abs(lat_values - center.y).argmin()
        return long_idx, lat_idx

    def _distance_between_two_points_in_km(self, point1, point2):
        line_string = LineString([shapely.geometry.Point(point1), shapely.geometry.Point(point2)])
        geod = Geod(ellps="WGS84")
        return geod.geometry_length(line_string) / 1000

    #based on https://stackoverflow.com/questions/72786576/how-to-scale-polygon-using-shapely
    def _find_nearest_point(self, geometry, ds, variable, feature_name):
        
        center = self._find_center_of_coordinates(geometry)
        long_idx, lat_idx = self._get_idx_for_nearest_point(ds, center)
 
        cropped_ds = ds[variable].isel(latitude=lat_idx, longitude=long_idx)

        point1 = (cropped_ds.longitude.values, cropped_ds.latitude.values)
        point2 = (center.x, center.y)

        distance_km = self._distance_between_two_points_in_km(point1, point2)

        print(f"WARNING: Did not find any points inside '{feature_name}', using nearest point instead, distance to point from feature center: {round(distance_km, 2)} km, using point '{cropped_ds.latitude.values}, {cropped_ds.longitude.values}'")

        return cropped_ds
        
    def _crop_dataset(self, ds, variable : str, feature):
        try:
            geometry = feature["geometry"]

            #cropped_ds would be a three-dimensional array
            #the first dimension is the ensembles, the second dimension represent each time step, the third dimension contain a 1-item array with the value for the given step
            cropped_ds : xr.core.dataarray.DataArray = ds[variable].rio.clip(geometries=[geometry])
       
            print(f"Feature '{feature['properties']['name']}' match on total {len(cropped_ds.latitude.values)*len(cropped_ds.longitude.values)} datapoints ({cropped_ds.latitude.values}, {cropped_ds.longitude.values}).")
            
            #all or a subset of the dimensions of length 1 would be removed
            squeezed_ds = cropped_ds.squeeze()

            #Crop the dataset to the geometry
            return squeezed_ds

        except rioxarray.exceptions.NoDataInBounds as e:
            return self._find_nearest_point(feature["geometry"], ds, variable, feature['properties']['name'])

    def _get_snapshot_leadtime_hours(self, dataset_starting_date : datetime, period_type : str, period_count : int, leadtime_interval : int):
        '''
        Getting all necessary leadtime hours when the forecast represents snapshot valuse (eg temperature).
        All available leadtime hours are required since we need to calculate some aggregate statistics for a given period. 
        '''
        lead_time_hours = []
        lead_time_hour = 0

        if period_type == 'M':
            dataset_ending_date = increment_months(dataset_starting_date, period_count)
        elif period_type[0] == 'W':
            dataset_ending_date = dataset_starting_date + timedelta(weeks=period_count)
        elif period_type == 'D':
            dataset_ending_date = dataset_starting_date + timedelta(days=period_count)
        
        forecast_length_hours = (dataset_ending_date - dataset_starting_date).days * 24

        while lead_time_hour < forecast_length_hours:
            lead_time_hour += leadtime_interval
            lead_time_hours.append(lead_time_hour)

        return lead_time_hours

    def _get_cumulative_leadtime_hours(self, dataset_starting_date : datetime, period_type : str, period_count : int):
        '''
        Getting only the necessary leadtime hours when the forecast represents cumulative values (eg precipitaiton).
        Returns list of hours since starting date into the future to forecast, at intervals specified by period type.
        '''
        lead_time_hours = []
        current_date = dataset_starting_date

        while len(lead_time_hours) < period_count:
            if period_type == 'M':
                next_date = increment_months(current_date, 1)
            elif period_type[0] == 'W':
                next_date = current_date + timedelta(weeks=1)
            elif period_type == 'D':
                next_date = current_date + timedelta(days=1)

            number_of_days_since_starting_date = (next_date - dataset_starting_date).days

            print(number_of_days_since_starting_date)

            lead_time_hour = 24 * int(number_of_days_since_starting_date)

            lead_time_hours.append(lead_time_hour)
            current_date = next_date

        return lead_time_hours

    def _get_mean_value_for_dimension_for_step_for_geometry(self, cropped_ds, feature, forecast_date : np.datetime64, step : int, value_converter):
        #returns all eseambles for the given forecast date and step
        points = cropped_ds.sel(time=forecast_date, step=step)

        #calculate the mean for all ensembles
        ensamble_mean = points.mean(keep_attrs=True)

        return PointValue(
            date = ensamble_mean.coords["valid_time"].values,
            value = value_converter(ensamble_mean),
            org_unit_id = str(feature["id"]),
            org_unit_name = feature["properties"]["name"]
        )

    # def save_calculated_results(self, df):
    #     # not sure if needed...
    #     df.to_csv(
    #         f"{self.output_folder}/{self.results_file_name}",  
    #         sep=";",
    #         index=False
    #     )

    def calculate(self):
        # open the seasonal forecast file downloaded from copernicus
        ds = xr.open_dataset(self.netcdf_file)

        print("\n--- STATS ---")
        print(f"reading ds file {self.netcdf_file}")
        print("number of datapoints: "+str(len(ds.longitude.values)*len(ds.latitude.values)))
        print("number of time-steps: "+str(len(ds.step)))
        print("number of ensembles: "+str(len(ds.number)))
        print("longitude values: "+str(ds.longitude.values))
        print("latitude values: "+str(ds.latitude.values))
        print("\n")

        ds.rio.write_crs("epsg:4326", inplace=True)

        result : List[PointValue] = []

        # get leadtime hours
        if self.total_sum_value:
            # cumulative forecast, only requires leadtime hours between each period type
            leadtime_hours = self._get_cumulative_leadtime_hours(self.forecast_date, self.period_type, self.period_count)
        else:
            # snapshot forecast, requires all available leadtime hours to calculate aggregate stats
            leadtime_interval = 6 # FIXME: hardcoded to 2m temperature for now
            leadtime_hours = self._get_snapshot_leadtime_hours(self.forecast_date, self.period_type, self.period_count, leadtime_interval)

        # convert to step values
        steps = [np.timedelta64(hour, 'h').astype('timedelta64[ns]') for hour in leadtime_hours]

        # loop over every feature
        for f in self.features:

            # crop dataset to this feature
            cropped_ds = self._crop_dataset(ds, self.variable, f,)

            # for every time-step
            for step in steps:
                r = self._get_mean_value_for_dimension_for_step_for_geometry(
                    cropped_ds=cropped_ds,
                    feature=f,
                    forecast_date=self.forecast_date,
                    step=step,
                    value_converter=converters[self.measurement_unit],
                )
                result.append(r)

        df = pd.DataFrame([ob.__dict__ for ob in result])

        if (self.total_sum_value):
            # the values at leadtime hour are cumulative, so have to take the current value minus the previous value
            df['diff'] = df.groupby('org_unit_id')['value'].transform(lambda x: x.diff())
            
            # for the first month entry of each 'org_unit_id', we set the orginal value as the original value
            df['diff'] = df['diff'].fillna(df['value'])
            df = df.drop(columns=['value'])
            df = df.rename(columns={'diff': 'value'})

            # since the leadtime hour corresponds to the amount of climate up to that hour, the previous month/week should be used
            if (self.period_type == "M"):
                df['period'] = (df['date'] - pd.DateOffset(months=1)).dt.to_period(self.period_type)
            elif (self.period_type[0] == "W"):
                df['period'] = (df['date'] - pd.DateOffset(weeks=1)).dt.to_period(self.period_type[0])
            elif (self.period_type == "D"):
                df['period'] = (df['date'] - pd.DateOffset(days=1)).dt.to_period(self.period_type)
                
            df = df.groupby(['org_unit_id', 'period', 'org_unit_name'])['value'].mean().reset_index()

        else:
            df['period'] = df['date'].dt.to_period(self.period_type[0])
            df = df.groupby(['org_unit_id', 'period', 'org_unit_name'])['value'].mean().reset_index()

        df['period'] = df['period'].apply(lambda p: str(p)) # convert to iso string

        df.sort_values(['org_unit_id', 'period'], inplace=True)

        print(df)

        return df

            
if __name__ == "__main__":
    #fetchData()
    seasonal = SeasonalForecastHandler('...')
    df = seasonal.calculate()
    print(df)
