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
        file_name_postfix="_temp_2025",
        indicator='2m_temperature', # 'total_precipitation'
        year=2025,
    )

    fetch_data = FetchCopernicusData(config)
    
    #fetch_data.get_data()
    #fetch_data.netcdf_file_name = 'netcdf/request_hash_cb8b2f753d86_precip_2025.nc'
    fetch_data.netcdf_file_name = 'netcdf/request_hash_9c0f81957200_temp_2025.nc'

    forecast_handler = fetch_data.get_forecast_handler('2025-02-01', 'W', 8) # period_count=3
    df = forecast_handler.calculate()
