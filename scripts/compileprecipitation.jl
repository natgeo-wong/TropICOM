using DrWatson
@quickactivate "TroPrecLS"

using Logging
using NASAPrecipitation

function compile(
    npd :: NASAPrecipitation.NASAPrecipitationDataset,
    geo :: GeoRegion,
)

    dtbeg = npd.start
    dtend = npd.stop
    dtvec = dtbeg : Day(1) : dtend
    lsd   = getLandSea(npd,geo)
    ndt   = length(dtvec)

    tmp  = zeros(Float32,length(lsd.lon),length(lsd.lat),48)
    prcp = zeros(Float32,length(lsd.lon),length(lsd.lat),48)

    for dt in dtvec


        @info "$(now()) - TroPrecLS - Loading GPM IMERG Half-Hourly data in Tropics GeoRegion for $dt"
        ds = read(npd,geo,dt)
        NCDatasets.load!(ds["precipitation"].var,tmp,:,:,:)
        
        for ilat = 1 : nlat, ilon = 1 : nlon
            prcp[ilon,ilat] += tmp[ilon,ilat]
        end

    end

    for ilat = 1 : nlat, ilon = 1 : nlon
        prcp[ilon,ilat] = prcp[ilon,ilat] / ndt
    end

    return prcp

end

npdhh = IMERGFinalHH(start=Date(2001,1,1),stop=Date(2020,12,31),path=datadir())
geo_TRP = GeoRegion("TRP")

prcp = compile(npdhh,geo_TRP)