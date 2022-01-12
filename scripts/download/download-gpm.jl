using DrWatson
@quickactivate "TroPrecLS"

addGeoRegions(srcdir("addgeorect.txt"))

npd = IMERGFinalHH(dtbeg=Date(2001,1,1),dtend=Date(2020,12,31),sroot=datadir())
geo = GeoRegion("TRP")
download(npd,geo)