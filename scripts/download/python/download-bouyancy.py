#!/usr/bin/env python
import cdsapi
import os

datadir = '/n/holyscratch01/kuang_lab/nwong/TroPrecLS/data/reanalysis/bouyancy/'

c = cdsapi.Client()

if not os.path.exists(datadir):
    os.makedirs(datadir)

for yr in range(1979,2020):
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': [
                'geopotential', 'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content',
                'specific_humidity', 'specific_rain_water_content', 'specific_snow_water_content',
                'temperature',
            ],
            'year': yr,
            'month': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            'time': '00:00',
            'pressure_level': '500',
        },
        datadir + 'era5-GLBx0.25-b_air-500hPa-' + str(yr) + '.nc')

for yr in range(1979,2020):
    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': [ '2m_temperature', 'surface_pressure', ],
            'year': yr,
            'month': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            'time': '00:00',
        },
        datadir + 'era5-GLBx0.25-b_sfc-' + str(yr) + '.nc')
