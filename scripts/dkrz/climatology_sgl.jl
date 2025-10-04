using DrWatson
@quickactivate "TropICOM"
using DelimitedFiles, Logging
using ERA5Reanalysis, RegionGrids

include(srcdir("common.jl"))

e5ds = ERA5Hourly(start=Date(1980),stop=Date(2024,12,31),path=datadir())
evar = SingleVariable("skt")

dtvec = e5ds.start : Day(1) : e5ds.stop; ndt = length(dtvec)
geo = GeoRegion("TRP",path=srcdir())
elon,elat = ERA5Reanalysis.nativelonlat()
ggrd = RegionGrid(geo,Point2.(elon,elat)); npnts = length(ggrd.ipoint)
vmat = zeros(npnts,24)
ipnt = minimum(ggrd.ipoint) : maximum(ggrd.ipoint)

for idt in 1 : ndt
    @info "$(now()) - S2DExploration - Extracting $(evar.name) data for the \"$(geo.ID)\" GeoRegion from the DKRZ servers for $(dtvec[idt])"; flush(stderr)
    ibeg = (idt-1) * 24 + 1
    iend = idt * 24
    gds = dkrz(e5ds,evar,dtvec[idt])
    vmat[:,:] .+= gds[evar.ID][ipnt,:]
    close(gds)
end

vmat ./= ndt
vmat = dropdims(mean(vmat,dims=2),dims=2)

save_climatology(geo,e5ds,evar,vmat,ERA5Reanalysis.era5Pressures(),ggrd)