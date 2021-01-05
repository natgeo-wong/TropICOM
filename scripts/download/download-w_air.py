#!/usr/bin/env python
import cdsapi
import os

datadir = '/n/holyscratch01/kuang_lab/nwong/TroPrecLS/data/reanalysis/w_air/'

c = cdsapi.Client()

if not os.path.exists(datadir):
    os.makedirs(datadir)

for yr in range(1979,2020):
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': 'vertical_velocity',
            'pressure_level': [
                1,2,3,5,7,10,20,30,50,70,
                100,125,150,175,200,225,250,300,350,
                400,450,500,550,600,650,700,750,775,
                800,825,850,875,900,925,950,975,1000
            ],
            'year': yr,
            'month': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            'area': [30, 0, -30, 360],
            'time': '00:00',
        },
        datadir + 'era5-TRPx0.25-w_air-' + str(yr) + '.nc')
