from fetch_data import FetchCopernicusData, FetchCopernicusDataConfig
import json

if __name__ == '__main__':

    file_path = 'data/sierraLeone.geojson'
    with open(file_path) as file:
        geojson_data = json.load(file)
    features = geojson_data['features']

    config = FetchCopernicusDataConfig(
        originating_centre="ecmwf",
        features=features,
        file_name_postfix="_2024",
        indicator='total_precipitation', #'2m_temperature',
        years=[2024],
        #forecast_length=...,
    )

    fetch_data = FetchCopernicusData(config)
    
    #fetch_data.get_data()
    fetch_data.netcdf_file_name = 'netcdf/request_hash_cb8b2f753d86_2025.nc'
    #fetch_data.netcdf_file_name = 'netcdf/request_hash_8ce608b0b18a.nc'

    forecast_handler = fetch_data.get_forecast_handler('2025-02-01', 'D') # period_count=3
    df = forecast_handler.calculate()
