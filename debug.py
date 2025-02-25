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
        #file_name_postfix="...",
        indicator='total_precipitation', # 2m_temperature
        year=2025,
    )

    fetch_data = FetchCopernicusData(config)
    
    #fetch_data.get_data()
    #fetch_data.netcdf_file_name = 'netcdf/ecmwf_total_precipitation_2025.nc'

    forecast_handler = fetch_data.get_forecast_handler('2025-02-01', 'M', 6) # period_count=3
    df = forecast_handler.calculate()
