using DrWatson
@quickactivate "TroPrecLS"

using Logging
using ERA5Reanalysis
using NCDatasets

function compile(
    e5ds :: ERA5Monthly,
    evar :: ERA5Variable,
    egeo :: ERA5Region,
)

    dtbeg = e5ds.start
    dtend = e5ds.stop
    dtvec = dtbeg : Day(1) : dtend
    lsd   = getLandSea(egeo,path=datadir("emask"))
    ndt   = length(dtvec)
    nlon  = length(lsd.lon)
    nlat  = length(lsd.lat)

    tint16 = zeros(Int16,nlon,nlat,288)
    tflt32 = zeros(Int16,nlon,nlat,288)
    var = zeros(Float32,nlon,nlat,288)

    for dt in dtvec

        @info "$(now()) - TroPrecLS - Loading $(uppercase(e5ds.lname)) $(evar.vname) in $(ereg.geo.name) (Horizontal Resolution: $(ereg.gres)) for $(year(dt)) ..."
        ids = read(e5ds,evar,egeo,dt)
        isc = ids[evar.varID].attrib["scale_factor"]
        iof = ids[evar.varID].attrib["add_offset"]
        imv = ids[evar.varID].attrib["missing_value"]
        ifv = ids[evar.varID].attrib["_FillValue"]
        NCDatasets.load!(ids[evar.varID].var,tint16,:,:,:)
        close(ids)

        ERA5Reanalysis.int2real!(tflt32,tint16,scale=isc,offset=iof,mvalue=imv,fvalue=ifv)

        for it = 1 : 288, ilat = 1 : nlat, ilon = 1 : nlon
            var[ilon,ilat,it] += tflt32[ilon,ilat,it]
        end

    end

    for it = 1 : 288, ilat = 1 : nlat, ilon = 1 : nlon
        var[ilon,ilat,it] = var[ilon,ilat,it] / ndt
    end

    var = dropdims(mean(reshape(var,nlon,nlat,24,12),dims=4),dims=4)

    @info "$(now()) - TroPrecLS - Saving compiled $(uppercase(e5ds.lname)) $(evar.vname) in $(ereg.geo.name) (Horizontal Resolution: $(ereg.gres)) ..."

    if evar <: SingleLevel
        fnc = datadir("era5mh-$(egeo.gstr)-$(evar.varID)-compiled.nc")
    else
        fnc = datadir("era5mh-$(egeo.gstr)-$(evar.varID)-$(evar.hPa)hPa-compiled.nc")
    end

    if isfile(fnc)
        rm(fnc,force=true)
    end
    ds  = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

	ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["time"]      = 24

    scale,offset = ncoffsetscale(var)

	nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

	ncvar = defVar(ds,evar.varID,Int16,("longitude","latitude","time"),attrib = Dict(
        "long_name"     => evar.lname,
        "full_name"     => evar.vname,
        "units"         => evar.units,
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
    ))

	nclongitude[:] = lsd.lon
	nclatitude[:]  = lsd.lat
	ncvar[:] = var

	close(ds)

end

e5ds_mh  = ERA5Monthly(start=Date(2001,1,1),stop=Date(2020,12,31),path=datadir(),hours=true)
egeo_TRP = ERA5Region(GeoRegion("TRP"))

evar = SingleVariable("skt")

compile(e5ds_mh,evar,egeo_TRP)