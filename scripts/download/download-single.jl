using DrWatson
@quickactivate "TroPrecLS"

using ERA5Reanalysis

e5ds = ERA5Monthly(dtbeg=Date(1979),dtend=Date(2021),eroot=datadir(),hours=true)
egeo = ERA5Region(GeoRegion("TRP"))

esgl = [
    SingleVariable("tsr"),SingleVariable("ttr"),SingleVariable("ssr"),
    SingleVariable("slhf"),SingleVariable("sshf"),SingleVariable("str"),
    SingleVariable("sst"),SingleVariable("t2m"),SingleVariable("skt"),
    SingleVariable("tcw"),SingleVariable("tcwv"),SingleVariable("sp"),
    SingleVariable("hcc"),SingleVariable("mcc"),SingleVariable("lcc"),SingleVariable("tcc"),
]

for evar in esgl
    download(e5ds,evar,egeo)
end