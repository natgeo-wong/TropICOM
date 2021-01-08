#!/usr/bin/env python
import cdsapi
import os

datadir = '/n/holyscratch01/kuang_lab/nwong/TroPrecLS/data/reanalysis/u_air/'

c = cdsapi.Client()

if not os.path.exists(datadir):
    os.makedirs(datadir)

for yr in range(1979,2020):
    for mo in [1,2,3,4,5,6,7,8,9,10,11,12]:
        c.retrieve(
            'reanalysis-era5-pressure-levels-monthly-means',
            {
                'format': 'netcdf',
                'product_type': 'monthly_averaged_reanalysis_by_hour_of_day',
                'variable': 'u_component_of_wind',
                'year': yr,
                'pressure_level': [
                    50,70,100,125,150,175,200,225,250,300,
                    350,400,450,500,550,600,650,700,750,775,
                    800,825,850,875,900,925,950,975,1000
                ],
                'month': mo,
                'time': [
                    '00:00', '01:00', '02:00',
                    '03:00', '04:00', '05:00',
                    '06:00', '07:00', '08:00',
                    '09:00', '10:00', '11:00',
                    '12:00', '13:00', '14:00',
                    '15:00', '16:00', '17:00',
                    '18:00', '19:00', '20:00',
                    '21:00', '22:00', '23:00',
                ],
                'area': [30, 0, -30, 360],
            },
            datadir + 'era5-TRPx0.25-u_air-' + str(yr) + str(mo).zfill(2) + '.nc')

for yr in range(1979,2020):
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': 'u_component_of_wind',
            'pressure_level': [
                50,70,100,125,150,175,200,225,250,300,
                350,400,450,500,550,600,650,700,750,775,
                800,825,850,875,900,925,950,975,1000
            ],
            'year': yr,
            'month': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            'area': [30, 0, -30, 360],
            'time': '00:00',
        },
        datadir + 'era5-TRPx0.25-u_air-' + str(yr) + '.nc')
