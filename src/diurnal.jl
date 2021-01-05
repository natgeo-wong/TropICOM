using Interpolations
using NCDatasets

function compileresavediurnalsfc(varname::AbstractString)

    tds = NCDataset(datadir("reanalysis/era5-TRPx0.25-$varname-sfc-hour.nc"))
    lon = tds["longitude"][:]*1; nlon = length(lon)
    lat = tds["latitude"][:]*1;  nlat = length(lat)
    var = tds[varname][:]*1;

    vart = zeros(25)

    for ilat = 1 : nlat, ilon = 1 : nlon

        t = (0:24) .+ lon[ilon]/15
        vart[1:24] = var[ilon,ilat,:]
        vart[end]  = var[ilon,ilat,1]
        itp = interpolate(vart,BSpline(Cubic(Periodic(OnGrid()))))
        stp = scale(itp,t)
        etp = extrapolate(stp,Periodic())

        var[ilon,ilat,:] .= etp[0:23]

    end

    if !isdir(datadir("reanalysis/diurnal"))
        mkpath(datadir("reanalysis/diurnal"))
    end

    fnc = datadir("reanalysis/diurnal/era5-TRPx0.25-$varname-sfc-diurnal.nc")
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"               => "CF-1.6",
        "history"                   => "2020-10-21 23:54:22 GMT by grib_to_netcdf-2.16.0: /opt/ecmwf/eccodes/bin/grib_to_netcdf -S param -o /cache/data4/adaptor.mars.internal-1603323636.0468726-6113-3-fb54412a-f0a5-4783-956d-46233705e403.nc /cache/tmp/fb54412a-f0a5-4783-956d-46233705e403-adaptor.mars.internal-1603323636.0477095-6113-1-tmp.grib",
    ))

    # Dimensions

    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["hour"]      = 24

    # Declare variables

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"                     => "degrees_east",
        "long_name"                 => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"                     => "degrees_north",
        "long_name"                 => "latitude",
    ))

    nchr  = defVar(ds,"level",Int8,("hour",),attrib = Dict(
        "units"                     => "hours",
        "long_name"                 => "Hours past 00:00 (GMT)",
    ))

    ncscale,offset = ncoffsetscale(var)

    varattribs = Dict(
        "scale_factor"  => ncscale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units"         => tds[varname].attrib["units"],
        "long_name"     => tds[varname].attrib["long_name"],
    )

    if haskey(tds[varname].attrib,"standard_name")
        varattribs["standard_name"] = tds[varname].attrib["standard_name"]
    end

    ncvar = defVar(ds,varname,Int16,("longitude","latitude","hour"),attrib = varattribs)

    nclon[:] = lon
    nclat[:] = lat
    nchr[:]  = collect(0:23)
    ncvar[:] = var

    close(ds)

    close(tds)

end

function compileresavediurnalpre(varname::AbstractString)

    tds = NCDataset(datadir("reanalysis/era5-TRPx0.25-$varname-hour.nc"))
    lon = tds["longitude"][:]*1; nlon = length(lon)
    lat = tds["latitude"][:]*1;  nlat = length(lat)
    lvl = tds["level"][:]*1;     nlvl = length(lvl)
    var = tds[varname][:]*1;

    vart = zeros(25)

    for ilvl = 1 : nlvl, ilat = 1 : nlat, ilon = 1 : nlon

        t = (0:24) .+ lon[ilon]/15
        vart[1:24] = var[ilon,ilat,ilvl,:]
        vart[end]  = var[ilon,ilat,ilvl,1]
        itp = interpolate(vart,BSpline(Cubic(Periodic(OnGrid()))))
        stp = scale(itp,t)
        etp = extrapolate(stp,Periodic())

        var[ilon,ilat,ilvl,:] .= etp[0:23]

    end

    if !isdir(datadir("reanalysis/diurnal"))
        mkpath(datadir("reanalysis/diurnal"))
    end

    fnc = datadir("reanalysis/diurnal/era5-TRPx0.25-$varname-diurnal.nc")
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions"               => "CF-1.6",
        "history"                   => "2020-10-21 23:54:22 GMT by grib_to_netcdf-2.16.0: /opt/ecmwf/eccodes/bin/grib_to_netcdf -S param -o /cache/data4/adaptor.mars.internal-1603323636.0468726-6113-3-fb54412a-f0a5-4783-956d-46233705e403.nc /cache/tmp/fb54412a-f0a5-4783-956d-46233705e403-adaptor.mars.internal-1603323636.0477095-6113-1-tmp.grib",
    ))

    # Dimensions

    ds.dim["longitude"] = nlon
    ds.dim["latitude"]  = nlat
    ds.dim["level"]     = nlvl
    ds.dim["hour"]      = 24

    # Declare variables

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"                     => "degrees_east",
        "long_name"                 => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"                     => "degrees_north",
        "long_name"                 => "latitude",
    ))

    ncpre = defVar(ds,"level",Int32,("level",),attrib = Dict(
        "units"                     => "millibars",
        "long_name"                 => "pressure_level",
    ))

    nchr  = defVar(ds,"hour",Int8,("hour",),attrib = Dict(
        "units"                     => "hours",
        "long_name"                 => "Hours past 00:00 (GMT)",
    ))

    ncscale,offset = ncoffsetscale(var)

    varattribs = Dict(
        "scale_factor"  => ncscale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units"         => tds[varname].attrib["units"],
        "long_name"     => tds[varname].attrib["long_name"],
    )

    if haskey(tds[varname].attrib,"standard_name")
        varattribs["standard_name"] = tds[varname].attrib["standard_name"]
    end

    ncvar = defVar(
        ds,varname,Int16,
        ("longitude","latitude","level","hour"),
        attrib = varattribs
    )

    nclon[:] = lon
    nclat[:] = lat
    ncpre[:] = lvl
    nchr[:]  = collect(0:23)
    ncvar[:] = var

    close(ds)

    close(tds)

end
