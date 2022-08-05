using DrWatson
@quickactivate "TroPrecLS"

using ERA5Reanalysis

e5ds = ERA5Monthly(start=Date(1979),stop=Date(2021),path=datadir(),hours=true)
egeo = ERA5Region(GeoRegion("TRP"))

esgl = [
    SingleVariable("tsr"),SingleVariable("ttr"),SingleVariable("ssr"),
    SingleVariable("slhf"),SingleVariable("sshf"),SingleVariable("str"),
    SingleVariable("sst"),SingleVariable("t2m"),SingleVariable("skt"),
    SingleVariable("tcw"),SingleVariable("tcwv"),SingleVariable("sp"),
    SingleVariable("hcc"),SingleVariable("mcc"),SingleVariable("lcc"),SingleVariable("cc"),
]

for evar in esgl
    download(e5ds,evar,egeo)
end