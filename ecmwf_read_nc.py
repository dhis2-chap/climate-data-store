import xarray as xr

#ds = xr.open_dataset("20240731000000-210h-oper-fc.grib2", engine="cfgrib")
print("Converting GRIB to netcdf..")
#ds.to_netcdf("20240731000000-210h-oper-fc.nc")

ds = xr.open_dataset("20240731000000-210h-oper-fc.nc")

print(len(ds))