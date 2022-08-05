using DrWatson
@quickactivate "TroPrecLS"

using Logging
using NASAPrecipitation
using NCDatasets

function compile(
    npd :: NASAPrecipitation.NASAPrecipitationDataset,
    geo :: GeoRegion,
)

    dtbeg = npd.start
    dtend = npd.stop
    dtvec = dtbeg : Day(1) : dtend
    lsd   = getIMERGlsd(geo,path=datadir("imergmask"))
    ndt   = length(dtvec)
    nlon  = length(lsd.lon)
    nlat  = length(lsd.lat)

    tmp  = zeros(Float32,nlon,nlat,48)
    prcp = zeros(Float32,nlon,nlat,48)

    for dt in dtvec

        @info "$(now()) - TroPrecLS - Loading GPM IMERG Half-Hourly data in Tropics GeoRegion for $dt ..."
        ids = read(npd,geo,dt)
        NCDatasets.load!(ids["precipitation"].var,tmp,:,:,:)
        close(ids)

        for it = 1 : 48, ilat = 1 : nlat, ilon = 1 : nlon
            prcp[ilon,ilat,it] += tmp[ilon,ilat,it]
        end

    end

    for it = 1 : 48, ilat = 1 : nlat, ilon = 1 : nlon
        prcp[ilon,ilat,it] = prcp[ilon,ilat,it] / ndt
    end

    @info "$(now()) - TroPrecLS - Saving compiled GPM IMERG Half-Hourly data in Tropics GeoRegion ..."

    fnc = datadir("imergfinalhh-$(geo.regID)-compiled.nc")
    if isfile(fnc)
        rm(fnc,force=true)
    end
    ds  = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));

	ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["time"]      = 48

	nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

	ncprcp = defVar(ds,"precipitation",Float32,("longitude","latitude","time"),attrib = Dict(
	    "units"     => "kg m**-2 s**-1",
	    "long_name" => "precipitation_rate",
		"full_name" => "Precipitation Rate",
    ))

	nclongitude[:] = lsd.lon
	nclatitude[:]  = lsd.lat
	ncprcp[:] = prcp

	close(ds)

end

npdhh = IMERGFinalHH(start=Date(2001,1,1),stop=Date(2020,12,31),path=datadir())
geo_TRP = GeoRegion("TRP")

compile(npdhh,geo_TRP)