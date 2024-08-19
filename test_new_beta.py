import cdsapi

import cdsapi

dataset = 'reanalysis-era5-pressure-levels'
request = {
    'product_type': ['reanalysis'],
    'variable': ['geopotential'],
    'year': [2024],
    'month': [3],
    'day': [1],
    'time': ['13:00'],
    'pressure_level': ['1000'],
    'data_format': 'grib',
}
target = 'download.grib'

client = cdsapi.Client()

print(request, client.key, target)



client.retrieve(dataset, request)#.download()



