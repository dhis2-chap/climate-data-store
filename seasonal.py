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
from pydantic import BaseModel
from pydantic.json_schema import JsonSchemaValue
from shapely import LineString, MultiPoint, Polygon
import matplotlib.pyplot as plt
from dateutil.relativedelta import relativedelta

class PeriodeValue(BaseModel):
    nameOfPeriode : str
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
    netcdf_file_name : str
    features : List[object]
    periode_type : str = "M" or "W-MON" or "D" or "W-SUN"
    aggreation_method : str = "mean" or "sum" or "max" or "min"
    output_file_postfix : str = ""
    measurement_unit : str = "kelvin" or "m"
    total_sum_value : bool = False


converters = {
    "K" : lambda x: x - 273.15, #converts kelvin to celsius
    "m" : lambda x: x * 1000 #converts meter to millimeter
}

class SeasonalForecastHandler():

    class Config:
        arbitrary_types_allowed=True

    def __init__(self, config : SeasonalForecastHandlerConfig):
        self.variable = config.variable
        self.netcdf_file_name = config.netcdf_file_name
        self.features = config.features
        self.periode_type = config.periode_type
        self.aggreation_method = config.aggreation_method
        self.output_file_postfix = config.output_file_postfix
        self.measurement_unit = config.measurement_unit
        self.total_sum_value = config.total_sum_value


    def kelvin_to_celsius(self, value):
        return value - 273.15
    
    def m_to_mm(self, value):
        return value * 1000

    def _group_by_moth(self, df : List[PointValue]):
        
        df['date'] = pd.to_datetime(df['date'])

        # Set 'time' as the index of the DataFrame
        df.set_index('date', inplace=True)

        # Group by month and calculate mean of all values
        monthly_stats = df.resample(self.periode_type).agg({'value': ["mean"]})

        
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
        

    def crop_dataset(self, ds, variable : str, feature):

        try:      
            geometry = feature["geometry"]

            #croppe_ds would be a three-dimensional array
            #the first dimension is the ensembles, the second dimension represent each time step, the third dimension contain a 1-item array with the value for the given step
            cropped_ds : xr.core.dataarray.DataArray = ds[variable].rio.clip(geometries=[geometry])
       
            print(f"Feature '{feature['properties']['name']}' match on total {len(cropped_ds.latitude.values)*len(cropped_ds.longitude.values)} datapoints ({cropped_ds.latitude.values}, {cropped_ds.longitude.values}).")
            
            #all or a subset of the dimensions of length 1 would be removed
            squeezed_ds = cropped_ds.squeeze()

            #Crop the dataset to the geometry
            return squeezed_ds

        except rioxarray.exceptions.NoDataInBounds as e:
            return self._find_nearest_point(feature["geometry"], ds, variable, feature['properties']['name'])
            
        


    def _get_mean_value_for_dimension_for_step_for_geometry(self, cropped_ds, feature, step : int, value_converter):
    
        
        #returns all eseambles for the given step
        points = cropped_ds.isel(step=step)

        #calculate the mean for all ensembles
        ensamble_mean = points.mean(keep_attrs=True)

        return PointValue(
            date = ensamble_mean.coords["valid_time"].values,
            value = value_converter(ensamble_mean),
            org_unit_id = str(feature["id"]),
            org_unit_name= feature["properties"]["name"]
        )

    def extract_previous_periode():
        pass

    def calculate(self):

        # open the seasonal forecast file downloaded from copernicus
        ds = xr.open_dataset(self.netcdf_file_name)

        print("\n--- STATS ---")
        print(f"reading ds file {self.netcdf_file_name}")
        print("number of datapoints: "+str(len(ds.longitude.values)*len(ds.latitude.values)))
        print("number of time-steps: "+str(len(ds.step)))
        print("number of ensembles: "+str(len(ds.number)))
        print("longitude values: "+str(ds.longitude.values))
        print("latitude values: "+str(ds.latitude.values))
        print("\n")

        ds.rio.write_crs("epsg:4326", inplace=True)


        result : List[PointValue] = []
        
        #loop over every featrue
        for f in self.features:

            #crop dataset to this feature
            cropped_ds = self.crop_dataset(ds, self.variable, f,)

            #for every time-step
            for i in range(len(ds.step)):
                r = self._get_mean_value_for_dimension_for_step_for_geometry(
                    cropped_ds=cropped_ds,
                    value_converter=converters[self.measurement_unit],
                    step=i,
                    feature=f
                )
                result.append(r)
        

        df = pd.DataFrame([ob.__dict__ for ob in result])


        if(self.total_sum_value):
            df['diff'] = df.groupby('org_unit_id')['value'].transform(lambda x: x.diff())
            # For the first month entry of each 'org_unit_id', we set the orginal value as the original value
            df['diff'] = df['diff'].fillna(df['value'])
            df = df.drop(columns=['value'])
            df = df.rename(columns={'diff': 'value'})

            #since the leadtime hour correspond amount of climate up to that houre, the previous month/week should be used
            if(self.periode_type == "M"):
                df['year_month'] = (df['date'] - pd.DateOffset(months=1)).dt.to_period(self.periode_type)
            else:
                df['year_month'] = (df['date'] - pd.DateOffset(weeks=1)).dt.to_period(self.periode_type)
                
            df = df.groupby(['org_unit_id', 'year_month', 'org_unit_name'])['value'].mean().reset_index()

        else:
            df['year_month'] = df['date'].dt.to_period(self.periode_type)
            df = df.groupby(['org_unit_id', 'year_month', 'org_unit_name'])['value'].mean().reset_index()

        df.sort_values(['org_unit_id', 'year_month'], inplace=True)

        print(df)

        df.to_csv(f"results/result_{datetime.today().strftime('%Y%m%d-%H-%M-%S')}_{self.variable}_{self.output_file_postfix}.csv",  sep=";")

            
if __name__ == "__main__":
    #fetchData()
    seasonal = SeasonalForecastHandler()
    seasonal.calculate()



#import matplotlib.pyplot as plt
