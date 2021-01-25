using DrWatson
using JLD2
@quickactivate "TroPrecLS"
using DelimitedFiles
using NCDatasets

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

function kc2int(kc::AbstractString)

    kclist = [
        "Af","Am","As","Aw","BWk","BWh","BSk","BSh","Cfa","Cfb","Cfc",
        "Csa","Csb","Csc","Cwa","Cwb","Cwc","Dfa","Dfb","Dfc","Dfd",
        "Dsa","Dsb","Dsc","Dsd","Dwa","Dwb","Dwc","Dwd","EF","ET"
    ]
    ikclist = kclist .== kc
    return findfirst(ikclist)

end

koppen = readdlm(datadir("Koeppen-Geiger-ASCII.txt"),skipstart=1)
klon = koppen[:,2]; klat = koppen[:,1]; kcls = koppen[:,3]; nk = length(kcls)
kcls2 = kc2int.(kcls)

lon = collect(-179.75:0.5:179.75); nlon = length(lon)
lat = collect(-89.75:0.5:89.75);   nlat = length(lat)
KCC = zeros(nlon,nlat)

for ik = 1 : nk

    ilat = klat[ik]; ilon = klon[ik]
    indlon = argmin(abs.(lon.-ilon))
    indlat = argmin(abs.(lat.-ilat))
    KCC[indlon,indlat] = kcls2[ik]

end
KCC[iszero.(KCC)] .= NaN
KCC_trop = deepcopy(KCC); KCC_trop[KCC_trop.>4] .= NaN;
KCC_arid = deepcopy(KCC); KCC_arid[KCC_arid.>8] .= NaN; KCC_arid[KCC_arid.<5] .= NaN
KCC_temph = deepcopy(KCC); KCC_temph[KCC_temph.>11] .= NaN; KCC_temph[KCC_temph.<9] .= NaN
KCC_temps = deepcopy(KCC); KCC_temps[KCC_temps.>14] .= NaN; KCC_temps[KCC_temps.<12] .= NaN
KCC_tempw = deepcopy(KCC); KCC_tempw[KCC_tempw.>17] .= NaN; KCC_tempw[KCC_tempw.<15] .= NaN

pplt.close(); f,axs = pplt.subplots(axwidth=9,aspect=3)
axs[1].pcolormesh(lon,lat,KCC_trop',levels=-1:5,cmap="Reds_r")
axs[1].pcolormesh(lon,lat,KCC_arid',levels=5:9,cmap="Yellow3_r")
axs[1].pcolormesh(lon,lat,KCC_temph',levels=6:12,cmap="Green2_r")
axs[1].pcolormesh(lon,lat,KCC_temps',levels=11:15,cmap="Brown1_r")
axs[1].pcolormesh(lon,lat,KCC_tempw',levels=12:18,cmap="Blue3_r")
axs[1].format(ylim=(-60,60),xlim=(-150,210))
f.savefig("test.png",transparent=false,dpi=200)

ds = NCDataset("koppen.nc","c")
ds.dim["longitude"] = nlon
ds.dim["latitude"]  = nlat
nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
    "units"                     => "degrees_east",
    "long_name"                 => "longitude",
))

nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
    "units"                     => "degrees_north",
    "long_name"                 => "latitude",
))
nctrop = defVar(ds,"koppenclass_tropics",Float32,("longitude","latitude"))
ncarid = defVar(ds,"koppenclass_arid",Float32,("longitude","latitude"))
nctemph = defVar(ds,"koppenclass_temperatehumid",Float32,("longitude","latitude"))
nctemps = defVar(ds,"koppenclass_temperatedrysummer",Float32,("longitude","latitude"))
nctempw = defVar(ds,"koppenclass_temperatedrywinter",Float32,("longitude","latitude"))
nclon[:] = lon
nclat[:] = lat
nctrop[:] = KCC_trop
ncarid[:] = KCC_arid
nctemph[:] = KCC_temph
nctemps[:] = KCC_temps
nctempw[:] = KCC_tempw

close(ds)
