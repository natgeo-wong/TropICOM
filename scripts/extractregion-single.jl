using DrWatson
@quickactivate "TroPrecLS"

using ERA5Reanalysis

e5ds = ERA5Monthly(dtbeg=Date(1979),dtend=Date(2021),eroot=datadir(),hours=true)
gvec = [
    ERA5Region(GeoRegion("DTP")),ERA5Region(GeoRegion("SEA")),
    ERA5Region(GeoRegion("AMZ")),ERA5Region(GeoRegion("TRA")),
    ERA5Region(GeoRegion("CRB"))
]

esgl = [
    SingleVariable("tsr"),SingleVariable("ttr"),SingleVariable("ssr"),
    SingleVariable("slhf"),SingleVariable("sshf"),SingleVariable("str"),
    SingleVariable("sst"),SingleVariable("t2m"),SingleVariable("skt"),
    SingleVariable("tcw"),SingleVariable("tcwv"),SingleVariable("sp"),
    SingleVariable("hcc"),SingleVariable("mcc"),SingleVariable("lcc"),SingleVariable("cc"),
]

for egeo in gvec, evar in esgl
    extract(e5ds,evar,egeo)
end