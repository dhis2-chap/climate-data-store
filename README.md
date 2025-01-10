## Get started

## Create a copernicus-account
Go to...

## Conda

### Install Miniconda

### Create Conda environment
conda env create -f environment.yml
conda activate climate-data-store

### Install climate data store

To update your environment to include updated packages

conda env update -f environment.yml

## Supported seasonal forecast sources

Supported seasonal forecast sources is located in file "./forecast_sources.json"

## Run script


fetch_data.py, takes care of fetching data from ECMWF
seasonal.py, is called from fetch_data.py and takes care of calcuating the response from ECMWF.
