using Glob
using NCDatasets
using NumericalIntegration
using Printf
using Statistics

function outstatname(
    expname :: AbstractString,
    config  :: AbstractString,
    isensemble :: Bool = false,
    member     :: Integer = 0
)

    if isensemble
    	  expname = "$(expname)-member$(@sprintf("%02d",member))"
    end

    return datadir(joinpath(expname,config,"OUT_STAT","DGW_TroPrecLS-$(expname).nc"))

end

function out2Dname(
    expname :: AbstractString,
    config  :: AbstractString,
    isensemble :: Bool = false,
    member     :: Integer = 0
)

    if isensemble
    	  expname = "$(expname)-member$(@sprintf("%02d",member))"
    end

	fnclist = glob(
		"DGW_TroPrecLS-$(expname)*.nc",
		joinpath(datadir(expname,config,"OUT_2D"))
	)

    return datadir(fnclist[1])

end

function retrievedims(
    expname :: AbstractString,
    config  :: AbstractString = "";
    isensemble :: Bool = false,
    member     :: Integer=0
)

    ds  = NCDataset(outstatname(expname,config,isensemble,member))
    z   = ds["z"][:]
    p   = ds["p"][:]
    t   = ds["time"][:]
    close(ds)

    return z,p,t

end

function retrievedims_fnc(fnc::AbstractString)

    ds  = NCDataset(fnc)
    z   = ds["z"][:]
    p   = ds["p"][:]
    t   = ds["time"][:]
    close(ds)

    return z,p,t

end

function retrievedims2D(
    expname :: AbstractString,
    config  :: AbstractString = "",
    isensemble :: Bool = false,
    member     :: Integer = 0
)

    ds = NCDataset(out2Dname(expname,config,isensemble,member))
	x = ds["x"][:]
    t = ds["time"][:]

	if haskey(ds,"y")
		  isy = true; y = ds["y"][:]
	else; isy = false
	end

    close(ds)

	if isy
    	  return x,y,t
	else; return x,t
	end

end

function retrievedims2D_fnc(fnc::AbstractString)

    ds = NCDataset(fnc)
    x = ds["x"][:]
    t = ds["time"][:]

	if haskey(ds,"y")
		  isy = true; y = ds["y"][:]
	else; isy = false
	end

    close(ds)

	if isy
    	  return x,y,t
	else; return x,t
	end

end

function retrievevar(
    varname :: AbstractString,
    expname :: AbstractString,
    config  :: AbstractString = "";
    isensemble :: Bool = false,
    member     :: Integer=0
)

    ds  = NCDataset(outstatname(expname,config,isensemble,member))
    var = ds[varname][:]
    close(ds)

    return var

end

function retrievevar_fnc(variable::AbstractString, fnc::AbstractString)

    ds  = NCDataset(fnc)
    var = ds[variable][:]
    close(ds)

    return var

end

function retrievevar2D(
    varname :: AbstractString,
    expname :: AbstractString,
    config  :: AbstractString = "";
    isensemble :: Bool = false,
    member     :: Integer=0
)

    ds  = NCDataset(out2Dname(expname,config,isensemble,member))
    var = ds[varname][:]
    close(ds)

    return var

end

function retrievevar2D_fnc(variable::AbstractString, fnc::AbstractString)

    ds  = NCDataset(fnc)
    var = ds[variable][:]
    close(ds)

    return var

end

function calcrh(QV,TAIR,P)

    RH = zeros(size(QV)); np = size(RH,1); nt = size(RH,2)

    for it = 1 : nt, ip = 1 : np
    	RH[ip,it] = QV[ip,it] / tair2qsat(TAIR[ip,it],P[ip]*100)
    end

    return RH

end

function tair2qsat(T,P)

    tb = T - 273.15
    if tb <= 0
    	esat = exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
    else
    	esat = exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
    end


    r = 0.622 * esat / max(esat,P-esat)
    return r / (1+r)

end

function calcswp(RH,QV,P)

	pvec = vcat(0,reverse(P)) * 100; nt = size(RH,2)
	QVsat = zeros(length(pvec))
	swp = zeros(nt)

    for it = 1 : nt
		QVsat[2:end] .= reverse(QV[:,it]) ./ reverse(RH[:,it])
		swp[it] = integrate(pvec,QVsat) / 9.81 / 1000
    end

	return swp

end

function t2d(t::Vector{<:Real}, days::Integer)

    tstep = round(Integer,(length(t)-1)/(t[end]-t[1]))
    t = mod.(t[(end-tstep+1):end],1); tmin = argmin(t)
    tshift = tstep-tmin+1; t = circshift(t,tshift)
    t = vcat(t[end]-1,t,t[1]+1)
    beg = days*tstep - 1

    return t*tstep,tstep,tshift,beg

end

function diurnal2D(data::AbstractVector,tstep::Integer,tshift::Integer)

    data = dropdims(mean(reshape(data,tstep,:),dims=2),dims=2);
    data = circshift(data,tshift)

    return vcat(data[end],data,data[1])

end

function diurnal3D(data::AbstractArray{<:Real,2},tstep::Integer,tshift::Integer)

    nz = size(data,1)
    data = dropdims(mean(reshape(data,nz,tstep,:),dims=3),dims=3);
    data = circshift(data,(0,tshift))

    return cat(dims=2,data[:,end],data,data[:,1])

end
