## Get started

## Create a ECMWF/Copernicus-account
- Create ECMWF account: https://www.ecmwf.int/user/login
- Follow CDSAPI setup: https://cds.climate.copernicus.eu/how-to-api

## Conda (recommended)

### Install Miniconda
https://docs.anaconda.com/miniconda/install/#quick-command-line-install

### Create Conda environment
Create conda environment named climate-data-store and required dependencies, by running:"
```bash
conda env create -f environment.yml
```

Active environment by running:
```bash
conda activate climate-data-store
```

## Supported seasonal forecast sources

Supported seasonal forecast sources is located in file: [./forecast_sources.json](./forecast_sources.json)

You could add more from: https://cds.climate.copernicus.eu/datasets/seasonal-original-single-levels?tab=download

## Run script
- [fetch_data.py](fetch_data.py), takes care of fetching data from ECMWF
- [seasonal.py](seasonal.py) is called from fetch_data.py and takes care of calcuating the response from ECMWF.


## Example
```bash
python fetch_data.py orgUnitsSingleSierra.geojson total_precipitation
```

Results will be located in folder "results"