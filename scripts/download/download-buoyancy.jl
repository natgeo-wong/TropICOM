using DrWatson
@quickactivate "TroPrecLS"

using ERA5Reanalysis

e5ds = ERA5Monthly(dtbeg=Date(1979),dtend=Date(2021),eroot=datadir())
egeo = ERA5Region(GeoRegion("GLB"),gres=0.25)

esgl = [SingleVariable("t2m"),SingleVariable("skt"),SingleVariable("d2m")]
epre = [
    PressureVariable("t",hPa=0,throw=false),PressureVariable("q",hPa=0,throw=false),
    PressureVariable("z",hPa=0,throw=false)
]

for evar in esgl
    download(e5ds,evar,egeo)
end

for evar in epre
    download(e5ds,evar,egeo,pall=true)
end