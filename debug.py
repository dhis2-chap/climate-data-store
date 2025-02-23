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
        file_name_postfix="-debug",
        period_type="D",
        indicator='total_precipitation', #'2m_temperature',
        #skip_download=False,
        #forecast_issued=None,
        #forecast_length=None,
    )

    fetch_data = FetchCopernicusData(config)
    fetch_data.get_data()