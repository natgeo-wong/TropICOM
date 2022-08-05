using DrWatson
@quickactivate "TroPrecLS"

using ERA5Reanalysis

e5ds = ERA5Monthly(start=Date(1979),stop=Date(2021),path=datadir(),hours=true)
egeo = ERA5Region(GeoRegion("TRP"))

epre = [
    PressureVariable("cc",hPa=0,throw=false),
    PressureVariable("ciwc",hPa=0,throw=false),
    PressureVariable("clwc",hPa=0,throw=false),
    PressureVariable("q",hPa=0,throw=false),
    PressureVariable("r",hPa=0,throw=false),
    PressureVariable("t",hPa=0,throw=false),
    PressureVariable("u",hPa=0,throw=false),
    PressureVariable("v",hPa=0,throw=false),
    PressureVariable("w",hPa=0,throw=false),
    PressureVariable("z",hPa=0,throw=false)
]

for evar in epre
    download(e5ds,evar,egeo,pall=true,ptop=50)
end