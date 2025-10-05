using DrWatson
@quickactivate "TropICOM"
using DelimitedFiles, Logging, Statistics
using ERA5Reanalysis, RegionGrids

include(srcdir("common.jl"))

e5ds = ERA5Monthly(start=Date(1980),stop=Date(2024,12,31),path=datadir())
evar = PressureVariable("t")

dtvec = e5ds.start : Year(1) : e5ds.stop; ndt = length(dtvec)
egeo = ERA5Region("TRP",path=srcdir(),native=true)
elon,elat = ERA5Reanalysis.nativelonlat()
ggrd = RegionGrid(egeo,Point2.(elon,elat)); npnts = length(ggrd.ipoint)
vmat = zeros(npnts,37,12)
ipnt = minimum(ggrd.ipoint) : maximum(ggrd.ipoint)

for idt in 1 : ndt
    @info "$(now()) - S2DExploration - Extracting $(evar.name) data for the \"$(egeo.ID)\" GeoRegion from the DKRZ servers for $(dtvec[idt])"; flush(stderr)
    gds = accessdkrz(e5ds,evar,dtvec[idt])
    vmat[:,:,:] .+= gds[evar.ID][ipnt,:,:]
    close(gds)
end

vmat ./= ndt
vmat = dropdims(mean(vmat,dims=3),dims=3)

save_climatology(e5ds,evar,egeo,vmat,ERA5Reanalysis.era5Pressures(),ggrd)