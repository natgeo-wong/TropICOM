### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 8102297a-50b7-11eb-0430-f79371a66174
begin
	using DrWatson
	using Pkg
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 9340fa4e-50b4-11eb-253e-ab01deb80456
begin
	@quickactivate "TroPrecLS"
	Pkg.instantiate()
	using NCDatasets
	using NumericalIntegration
	using Statistics
	
	using ImageShow, FileIO, ImageMagick
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ addc35d6-50b3-11eb-02dc-452ced2a45ef
md"
# X. Plotting Functions for SAM Output

Here I am creating and testing functions that will be used for the extraction and plotting of SAM statistical output data.
"

# ╔═╡ 19174960-50cf-11eb-12a3-cf977e483262
md"
### 1. Retrieving and Manipulating Statistical Output

We first create the function `retrievedims` to extract the dimensional data: (1) height (in both vertical `z` and pressure `p` coordinates) and (2) time `t`.  Next, we define the function `retrievevar` to extract variables of interest (one variable at a time):
"

# ╔═╡ c1489ae0-5114-11eb-3a56-5b75d263ae63
function retrievedims(experiment::AbstractString, config::AbstractString)
	
	rce = NCDataset(datadir(joinpath(
        experiment,config,"OUT_STAT",
        "RCE_TroPrecLS-$(experiment)-test.nc"
    )))
	
    z = rce["z"][:]
    p = rce["p"][:]
	t = rce["time"][:]
	
    close(rce)

    return z,p,t
	
end

# ╔═╡ 891e5992-50cf-11eb-208e-71e78380caf7
function retrievevar(
    variable::AbstractString,
	experiment::AbstractString,
	config::AbstractString
)
	
	rce = NCDataset(datadir(joinpath(
        experiment,config,"OUT_STAT",
        "RCE_TroPrecLS-$(experiment)-test.nc"
    )))
	
    var = rce[variable][:]
	
    close(rce)

    return var
	
end

# ╔═╡ edc78916-5115-11eb-00a0-498face6f531
md"
Relative humidity in SAM is plotted against the liquid saturation point.  However, we are able to calculate relative humidity using the ice-saturation point for air temperatures below 0ºC, using specific humidity and atmospheric temperature data.
"

# ╔═╡ a8331f82-5117-11eb-158b-b3073e67affd
md"
From the relative humidity that we calculate, we then can calculate column saturation fraction `CSF`, or the overall relative humidity of the column
"

# ╔═╡ dacb2ffe-5117-11eb-33bf-ff68ff8c6504
function calccsf(RH,P)
	
	pvec = vcat(0,reverse(P)); nt = size(RH,2)
	pint = integrate(pvec,ones(length(pvec)))
	RHtmp = zeros(length(pvec))
    csf = zeros(nt)
	
    for it = 1 : nt
		RHtmp[2:end] .= RH[:,it]
        csf[it] = integrate(pvec,RHtmp) / pint
    end
	
	return csf
	
end

# ╔═╡ 2fb20134-5119-11eb-1d07-db7f8c326dda
md"
Let us now define functions that will plot the time series 2D and 3D variables on a given axes object `axs`, and axes number `ii`, between days `dbeg` and `dend`.  The function for 3D plot will return information that can be used to plot a colorbar.  The 3D timeseries can also accept optional colormap `cmapname` and levels `lvl` input.
"

# ╔═╡ 4324677a-5119-11eb-24f1-8d80c5a15ceb
function plot2Dtimeseries(axs,ii,t,var;dbeg,dend)
	axs[ii].plot(t,var,lw=1,c="k")
	axs[ii].format(xlim=(dbeg,dend))
end

# ╔═╡ e006175a-5119-11eb-02ba-b50a41b540cb
function plot3Dtimeseries(axs,ii,t,p,var;lvl=[],cmapname="Fire",dbeg,dend)
	
	if isempty(lvl)
		  c = axs[ii].contourf(t,p,var,cmap=cmapname)
	else; c = axs[ii].contourf(t,p,var,cmap=cmapname,levels=lvl)
	end
	
	axs[ii].format(xlim=(dbeg,dend),ylim=(maximum(p),minimum(p)))
	axs[ii].colorbar(c,loc="r")
	
end

# ╔═╡ f7cf698a-511f-11eb-0276-7d4d7ae912e2
md"
Following this, we define functions that allow us to calculate and plot the diurnal cycle of the variables considered based on the last `ndays` of the simulation.
"

# ╔═╡ 23322f5a-5126-11eb-3d55-f3135e8a8bdc
function t2d(t::Vector{<:Real}, days::Integer)

    tstep = round(Integer,(length(t)-1)/(t[end]-t[1]))
    t = mod.(t[(end-tstep+1):end],1); tmin = argmin(t)
    tshift = tstep-tmin+1; t = circshift(t,tshift)
    t = vcat(t[end]-1,t,t[1]+1)
	beg = days*tstep - 1

    return t*tstep,tstep,tshift,beg

end

# ╔═╡ bea41e8c-5127-11eb-0715-e55c3f2a05ab
function diurnal2D(data::AbstractVector,tstep::Integer,tshift::Integer)

    data = dropdims(mean(reshape(data,tstep,:),dims=2),dims=2);
    data = circshift(data,tshift)

    return vcat(data[end],data,data[1])

end

# ╔═╡ 55946b50-5126-11eb-14e4-570a92dee964
function diurnal3D(data::AbstractArray{<:Real,2},tstep::Integer,tshift::Integer)

    nz = size(data,1)
    data = dropdims(mean(reshape(data,nz,tstep,:),dims=3),dims=3);
    data = circshift(data,(0,tshift))

    return cat(dims=2,data[:,end],data,data[:,1])

end

# ╔═╡ e0d83a56-5127-11eb-26a2-9dc0aa1c9822
function plot2Ddiurnal(axs,ii,t,var;subtractm=true)
	
	mvar = mean(var)
	if subtractm
		  axs[ii].plot(t,var .- mvar,lw=1)
	else; axs[ii].plot(t,var,lw=1)
	end
	
	axs[ii].format(xlim=(0,24))
	if !subtractm; axs[ii].format(ylim=(0,2.5)) end
	
end

# ╔═╡ a8d78482-5126-11eb-2dd1-8751cbd884a7
function plot3Ddiurnal(axs,ii,t,p,var;lvl=[],cmapname="Fire")
	
	mvar = dropdims(mean(var,dims=2),dims=2)
	axs[ii].plot(mvar,p,lw=1)
	
	if isempty(lvl)
		  c = axs[ii+1].contourf(t,p,var,cmap=cmapname)
	else; c = axs[ii+1].contourf(t,p,var,cmap=cmapname,levels=lvl)
	end
	
	axs[ii+1].format(xlim=(0,24),ylim=(maximum(p),minimum(p)))
	axs[ii+1].colorbar(c,loc="r")
	
end

# ╔═╡ 6bc27074-5115-11eb-29ed-bf79782322d0
md"
### 2. 2D Statistical Output

The following 2D statistical output are of interest:
* Sea-surface temperature | `SST`
* Surface pressure | `Ps`
* Precipitation | `Prec`
* Precipitable Water | `PW`
* Sensible/Latent Heat Flux | `SHF`/`LHF`
* Surface Net Shortwave/Longwave Flux | `SWNS`/`LWNS`
"

# ╔═╡ 504cbace-5128-11eb-1d27-879bffb48098
md"
### 3. Testing functionality
"

# ╔═╡ 5ffd515e-5128-11eb-3b11-815af069d22f
begin
	exp = "3SWTGamExp0"; config = "damping02d00"
	z,p,t = retrievedims(exp,config)
	v2D = retrievevar("PREC",exp,config)
	v3D = retrievevar("WWTG",exp,config) * 3600
	size(v2D), size(v3D)
end

# ╔═╡ 6b0b697c-5117-11eb-1c0f-13c2840a65d2
function tair2qsat(T,P)

	tb = T - 273.15
	if tb <= 0
		esat = exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
	else
		esat = exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
	end


	r = 0.622 * esat / max(esat,p-esat)
	return r / (1+r)

end

# ╔═╡ 48f8e67e-5116-11eb-0037-3b6630655917
function calcrh(QV,TAIR,P)

	RH = zeros(size(RH)); np = size(RH,1); nt = size(RH,2)

	for it = 1 : nt, ip = 1 : np
		RH[ip,it] = QV[ip,it] / tair2qsat(TAIR[ip,it],p[ip]*100)
	end

	return RH

end

# ╔═╡ f440d5ca-5128-11eb-2574-d58f2f7d8fc3
begin
	
	pplt.close(); fts,axsts = pplt.subplots(nrows=2,aspect=2)
	
	lvls=vcat(-5:-1,1:5)*20
	plot2Dtimeseries(axsts,1,t.-80,v2D,dbeg=30,dend=40)
	plot3Dtimeseries(axsts,2,t.-80,p,v3D,dbeg=30,dend=60,lvl=lvls,cmapname="RdBu_r")
	axsts[1].format(ylim=(0,2))
	
	fts.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
	
end

# ╔═╡ 9eab2f9c-5129-11eb-07fc-2d9b382f5d49
begin
	
# 	td,tstep,tshift,beg = t2d(t,100);
# 	v2Dd = diurnal2D(v2D[(end-beg):end],tstep,tshift);
# 	v3Dd = diurnal3D(v3D[:,(end-beg):end],tstep,tshift);
	
# 	arr = [[0,1,1,1],[2,3,3,3]]
# 	pplt.close(); fdh,axsdh = pplt.subplots(arr,nrows=2,aspect=2)
	
# 	plot2Ddiurnal(axsdh,1,td,v2Dd,subtractm=false)
# 	plot3Ddiurnal(axsdh,2,td,p,v3Dd)
	
# 	fdh.savefig("test2.png",transparent=false,dpi=200)
# 	load("test2.png")
	
end

# ╔═╡ Cell order:
# ╟─addc35d6-50b3-11eb-02dc-452ced2a45ef
# ╟─8102297a-50b7-11eb-0430-f79371a66174
# ╠═9340fa4e-50b4-11eb-253e-ab01deb80456
# ╟─19174960-50cf-11eb-12a3-cf977e483262
# ╠═c1489ae0-5114-11eb-3a56-5b75d263ae63
# ╠═891e5992-50cf-11eb-208e-71e78380caf7
# ╟─edc78916-5115-11eb-00a0-498face6f531
# ╠═48f8e67e-5116-11eb-0037-3b6630655917
# ╠═6b0b697c-5117-11eb-1c0f-13c2840a65d2
# ╟─a8331f82-5117-11eb-158b-b3073e67affd
# ╠═dacb2ffe-5117-11eb-33bf-ff68ff8c6504
# ╟─2fb20134-5119-11eb-1d07-db7f8c326dda
# ╠═4324677a-5119-11eb-24f1-8d80c5a15ceb
# ╠═e006175a-5119-11eb-02ba-b50a41b540cb
# ╟─f7cf698a-511f-11eb-0276-7d4d7ae912e2
# ╠═23322f5a-5126-11eb-3d55-f3135e8a8bdc
# ╠═bea41e8c-5127-11eb-0715-e55c3f2a05ab
# ╠═55946b50-5126-11eb-14e4-570a92dee964
# ╠═e0d83a56-5127-11eb-26a2-9dc0aa1c9822
# ╠═a8d78482-5126-11eb-2dd1-8751cbd884a7
# ╟─6bc27074-5115-11eb-29ed-bf79782322d0
# ╟─504cbace-5128-11eb-1d27-879bffb48098
# ╠═5ffd515e-5128-11eb-3b11-815af069d22f
# ╠═f440d5ca-5128-11eb-2574-d58f2f7d8fc3
# ╠═9eab2f9c-5129-11eb-07fc-2d9b382f5d49
