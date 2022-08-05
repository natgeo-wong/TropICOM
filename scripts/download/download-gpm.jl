using DrWatson
@quickactivate "TroPrecLS"

addGeoRegions(srcdir("addgeorect.txt"))

npd = IMERGFinalHH(start=Date(2001,1,1),stop=Date(2020,12,31),path=datadir())
geo = GeoRegion("TRP")
download(npd,geo)