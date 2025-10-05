using Dates
using ERA5Reanalysis
using NCDatasets
using Printf

expdir(args...) = joinpath(projectdir("exp"), args...)
prmdir(args...) = joinpath(projectdir("exp","prm"), args...)
lsfdir(args...) = joinpath(projectdir("exp","lsf"), args...)
snddir(args...) = joinpath(projectdir("exp","snd"), args...)

rundir(args...) = joinpath(projectdir("run"), args...)

## DateString Aliasing
yr2str(date::TimeType)   = Dates.format(date,dateformat"yyyy")
yrmo2dir(date::TimeType) = Dates.format(date,dateformat"yyyy/mm")
yrmo2str(date::TimeType) = Dates.format(date,dateformat"yyyymm")
ymd2str(date::TimeType)  = Dates.format(date,dateformat"yyyymmdd")

function read_climatology(
    e5ds :: ERA5Monthly,
    evar :: ERA5Variable,
    egeo :: ERA5Region;
    days :: Int = 0
)

    if iszero(days)
        return NCDataset(joinpath(e5ds.path,"climatology",
            e5ds.ID * "-" * egeo.string * "-" * evar.ID * "-" *
            yrmo2str(e5ds.start) * "-" * yrmo2str(e5ds.stop) * ".nc"
        ))
    else
        return NCDataset(joinpath(e5ds.path,"climatology",
            e5ds.ID * "-" * egeo.string * "-" * evar.ID * "-" *
            yrmo2str(e5ds.start) * "-" * yrmo2str(e5ds.stop) * "-" *
            "smooth$(@sprintf("%02d",days))days.nc"
        ))
    end

end

function save_climatology(
    e5ds :: ERA5Monthly,
    evar :: SingleLevel,
    egeo :: ERA5Region,
    data :: Vector{<:Real},
    ggrd :: RegionGrid;
    days :: Int = 0
)

    npts = length(data)
    if iszero(days)
        fnc = joinpath(e5ds.path,"climatology",
            e5ds.ID * "-" * egeo.string * "-" * evar.ID * "-" *
            yrmo2str(e5ds.start) * "-" * yrmo2str(e5ds.stop) * ".nc"
        )
    else
        fnc = joinpath(e5ds.path,"climatology",
            e5ds.ID * "-" * egeo.string * "-" * evar.ID * "-" *
            yrmo2str(e5ds.start) * "-" * yrmo2str(e5ds.stop) * "-" *
            "smooth$(@sprintf("%02d",days))days.nc"
        )
    end

    if isfile(fnc); rm(fnc) end
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now())",
        "comments"    => "These NetCDF files were created in the same format that data is saved on the Climate Data Store",
    ))

    ds.dim["values"] = npts

    nclon = defVar(ds,"longitude",Float64,("values",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float64,("values",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    ncvar = defVar(ds,evar.ID,Float64,("values",),attrib = Dict(
        "long_name" => evar.long,
        "full_name" => evar.name,
        "units"     => evar.units,
    ))

    nclon[:] = ggrd.lon
    nclat[:] = ggrd.lat
    ncvar[:] = data

    close(ds)
    
end

function save_climatology(
    e5ds :: ERA5Monthly,
    evar :: PressureLevel,
    egeo :: ERA5Region,
    data :: Matrix{<:Real},
    lvls :: Vector{Int},
    ggrd :: RegionGrid;
    days :: Int = 0
)

    npts,nlvls = size(data)

    if iszero(days)
        fnc = joinpath(e5ds.path,"climatology",
            e5ds.ID * "-" * egeo.string * "-" * evar.ID * "-" *
            yrmo2str(e5ds.start) * "-" * yrmo2str(e5ds.stop) * ".nc"
        )
    else
        fnc = joinpath(e5ds.path,"climatology",
            e5ds.ID * "-" * egeo.string * "-" * evar.ID * "-" *
            yrmo2str(e5ds.start) * "-" * yrmo2str(e5ds.stop) * "-" *
            "smooth$(@sprintf("%02d",days))days.nc"
        )
    end

    if isfile(fnc); rm(fnc) end
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now())",
        "comments"    => "These NetCDF files were created in the same format that data is saved on the Climate Data Store",
    ))

    ds.dim["values"] = npts
    ds.dim["levels"] = nlvls

    nclon = defVar(ds,"longitude",Float64,("values",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float64,("values",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nclvl = defVar(ds,"pressures",Int64,("levels",),attrib = Dict(
        "units"     => "hPa",
        "long_name" => "Pressure Levels",
    ))

    ncvar = defVar(ds,evar.ID,Float64,("values","levels",),attrib = Dict(
        "long_name" => evar.long,
        "full_name" => evar.name,
        "units"     => evar.units,
    ))

    nclon[:]   = ggrd.lon
    nclat[:]   = ggrd.lat
    nclvl[:]   = lvls
    ncvar[:,:] = data

    close(ds)
    
end