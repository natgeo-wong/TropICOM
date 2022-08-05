using DrWatson
@quickactivate "TroPrecLS"

using Logging
using ERA5Reanalysis
using NCDatasets
using Statistics

function compile(
    e5ds :: ERA5Monthly,
    evar :: ERA5Variable,
    egeo :: ERA5Region,
)

    dtbeg = e5ds.start
    dtend = e5ds.stop
    dtvec = dtbeg : Year(1) : dtend
    lsd   = getLandSea(egeo,path=datadir("emask"))
    ndt   = length(dtvec)
    nlon  = length(lsd.lon)
    nlat  = length(lsd.lat)

    tint16 = zeros(Int16,nlon,nlat,288)
    tflt32 = zeros(Float32,nlon,nlat,288)
    var = zeros(Float32,nlon,nlat,288)

    for dt in dtvec

        @info "$(now()) - TroPrecLS - Loading $(uppercase(e5ds.lname)) $(evar.vname) in $(egeo.geo.name) (Horizontal Resolution: $(egeo.gres)) for $(year(dt)) ..."
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

    @info "$(now()) - TroPrecLS - Saving compiled $(uppercase(e5ds.lname)) $(evar.vname) in $(egeo.geo.name) (Horizontal Resolution: $(egeo.gres)) ..."

    fol = "compiled"
    if !isdir(datadir(fol)); mkpath(datadir(fol)) end
    if typeof(evar) <: SingleLevel
        fnc = datadir(fol,"era5mh-$(egeo.gstr)-$(evar.varID)-compiled.nc")
    else
        fnc = datadir(fol,"era5mh-$(egeo.gstr)-$(evar.varID)-$(evar.hPa)hPa-compiled.nc")
    end

    if isfile(fnc)
        rm(fnc,force=true)
    end
    ds  = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

	ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["time"]      = 24

    scale,offset = ERA5Reanalysis.ncoffsetscale(var)

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
	ncvar.var[:] = ERA5Reanalysis.real2int16(var,scale,offset)

	close(ds)

end

e5ds_mh  = ERA5Monthly(start=Date(1979,1,1),stop=Date(2020,12,31),path=datadir(),hours=true)
egeo_TRP = ERA5Region(GeoRegion("TRP"))

esgl = [
    SingleVariable("tsr"),SingleVariable("ttr"),SingleVariable("ssr"),
    SingleVariable("slhf"),SingleVariable("sshf"),SingleVariable("str"),
    SingleVariable("sst"),SingleVariable("t2m"),SingleVariable("skt"),
    SingleVariable("tcw"),SingleVariable("tcwv"),SingleVariable("sp"),
    SingleVariable("hcc"),SingleVariable("mcc"),SingleVariable("lcc"),SingleVariable("cc"),
]

for evar in esgl
    compile(e5ds_mh,evar,egeo_TRP)
end

for ip in era5Pressures()
    epre = [
        PressureVariable("cc",  hPa=ip,throw=false),
        PressureVariable("ciwc",hPa=ip,throw=false),
        PressureVariable("clwc",hPa=ip,throw=false),
        PressureVariable("q",   hPa=ip,throw=false),
        PressureVariable("r",   hPa=ip,throw=false),
        PressureVariable("t",   hPa=ip,throw=false),
        PressureVariable("u",   hPa=ip,throw=false),
        PressureVariable("v",   hPa=ip,throw=false),
        PressureVariable("w",   hPa=ip,throw=false),
        PressureVariable("z",   hPa=ip,throw=false)
    ]
    for evar in epre
        compile(e5ds_mh,evar,egeo_TRP)
    end
end