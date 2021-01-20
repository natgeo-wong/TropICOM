### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 36a4e004-5a26-11eb-1ba4-db83d18c7845
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 3a8bbab2-5a26-11eb-20e0-d551c2a68858
begin
	@quickactivate "TroPrecLS"
	using GeoRegions
	using Statistics
	using NCDatasets
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("common.jl"))
	include(srcdir("samlsf.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 582bcb26-5a22-11eb-3b52-09e3ec08b3f9
md"
# 5. Reanalysis Moisture-Fluxes

In this notebook, we calculate and find the vertical profile of mean moisture-fluxes into our predefined domains.  These will then be used to force the mean-state of our SAM models into wetter/drier states, such that we can reconstruct the P-r curve in Bretherton et al. [2006].

We use mean values of ERA5 specific humidity, zonal and meridional wind to calculate the moisture flux vertical profile at each level, since this is not given as raw data in the ERA5 CDS data store.  By using only the mean-values calculated, and not the raw hourly output, we assume the following

$$\overline{qv}
= \overline{q}\overline{v} + \overline{q'}\overline{v'}
\approx \overline{q}\overline{v}$$

i.e., that $\overline{q'}\overline{v'} \ll \overline{q}\overline{v}$.
"

# ╔═╡ 2b736522-5a33-11eb-0eea-473a5b8566f3
md"
### A. Calculating Moisture Flux Convergence

We calculate the gradients and divergence using finite difference.  Then, we calculate the moisture flux convergence via the following:

$$\text{MFC} = - \vec{u}\cdot\nabla q - q\nabla\cdot\vec{u}$$
"

# ╔═╡ 28ffa404-5a33-11eb-3dc7-7d716d6cf06f
function ddx(
	data::AbstractArray{<:Real,3},
	lon::AbstractVector{<:Real},
	lat::AbstractVector{<:Real}
)
	
	nx = length(lon); ny = length(lat); nlvl = size(data,3)
	rx = 6378e3; dx = rx * (lon[2]-lon[1]) * pi / 180 * cosd.(lat)
	ddatadx = zeros(nx,ny,nlvl)
	
	for ilvl = 1 : nlvl, iy = 1 : ny, ix = 2 : (nx-1)

        ddatadx[ix,iy,ilvl] = (data[ix+1,iy,ilvl] - data[ix-1,iy,ilvl]) / (2 * dx[iy])

    end

    for iy = 1 : ny, ilvl = 1 : nlvl

        ddatadx[1,iy,ilvl] = (data[2,iy,ilvl] - data[1,iy,ilvl]) / dx[iy]
        ddatadx[nx,iy,ilvl] = (data[nx,iy,ilvl] - data[nx-1,iy,ilvl]) / dx[iy]

    end
	
	return ddatadx
	
end

# ╔═╡ 9e01c7ea-5a34-11eb-1c41-8b991a399584
function ddy(
	data::AbstractArray{<:Real,3},
	lat::AbstractVector{<:Real}
)
	
	nx = size(data,1); ny = length(lat); nlvl = size(data,3)
	ry = 6357e3; dy = ry * (lat[2]-lat[1]) * pi / 180
	ddatady = zeros(nx,ny,nlvl)
	
	for ilvl = 1 : nlvl, iy = 2 : (ny-1), ix = 1 : nx

        ddatady[ix,iy,ilvl] = (data[ix,iy+1,ilvl] - data[ix,iy-1,ilvl]) / (2 * dy)

    end

    for ilvl = 1 : nlvl, ix = 1 : nx

        ddatady[ix,1,ilvl] = (data[ix,1,ilvl] - data[ix,2,ilvl]) / dy
        ddatady[ix,ny,ilvl] = (data[ix,ny-1,ilvl] - data[ix,ny,ilvl]) / dy

    end
	
	return ddatady
	
end

# ╔═╡ abde634c-5a83-11eb-00d6-69b233e0c00f
md"
### B. Extracting Moisture Fluxes for Tropical REgion
"

# ╔═╡ 4927b4c2-5a8a-11eb-37e4-53efb9445b23
function saveqfc(qfc,rlon,rlat,lvl,regID)
	
	fnc = datadir("reanalysis/era5-$(regID)x0.25-qfc_air.nc")
	if isfile(fnc); rm(fnc) end
	ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));
	
	ds.dim["longitude"] = length(rlon)
    ds.dim["latitude"]  = length(rlat)
    ds.dim["level"]  = length(lvl)
	
	scale,offset = ncoffsetscale(qfc)
	
	nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))
	
	nclevel = defVar(ds,"level",Int32,("level",),attrib = Dict(
        "units"     => "millibars",
        "long_name" => "pressure_level",
    ))
	
	ncqfc = defVar(ds,"qfc_air",Int16,("longitude","latitude","level"),attrib = Dict(
        "scale_factor"  => scale,
        "add_offset"    => offset,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
        "units" 	    => "kg kg**-1 s**-1",
		"long_name"     => "moisture_flux_convergence",
    ))
	
	nclongitude[:] = rlon
	nclatitude[:]  = rlat
	nclevel[:]     = lvl
	ncqfc[:]       = qfc
	
	close(ds)
	
end

# ╔═╡ 87a146d4-5a29-11eb-1571-b9d9266ee54b
function calcqfc(regID::AbstractString)
	
	ds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-q_air.nc"))
    lon = ds["longitude"][:]; nlon = length(lon)
    lat = ds["latitude"][:];  nlat = length(lat)
	lvl = ds["level"][:];	  nlvl = length(lvl)
    q_air = ds["q_air"][:]*1
    close(ds)
	
	ds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-u_air.nc"))
    u_air = ds["u_air"][:]*1
    close(ds)
	
	ds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-v_air.nc"))
    v_air = ds["v_air"][:]*1
    close(ds)
	
	tmp = zeros(nlon,nlat,nlvl)
	rlon,rlat,rinfo = gregiongridvec(regID,lon,lat)
	
	ru = regionextractgrid!(u_air,rinfo,tmp)
	rv = regionextractgrid!(v_air,rinfo,tmp)
	rq = regionextractgrid!(q_air,rinfo,tmp)
	
	ru_x = ddx(ru,rlon,rlat); rv_y = ddy(rv,rlat)
	rq_x = ddx(rq,rlon,rlat); rq_y = ddy(rq,rlat)
	
	qfc = -((ru_x + rv_y) .* rq .+ ru .* rq_x .+ rv .* rq_y)
	
	saveqfc(qfc,rlon,rlat,lvl,regID)
	
	return rlon,rlat,lvl,qfc
	
end

# ╔═╡ 795a822a-5a29-11eb-37bb-3b6d62cec6be
regID = "TRP"

# ╔═╡ e3c519fc-5a35-11eb-0d34-cd4ca8b50bf9
if isfile(datadir("reanalysis/era5-$(regID)x0.25-qfc_air.nc"))
	qfcds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-qfc_air.nc"))
	lon = qfcds["longitude"][:]
	lat = qfcds["latitude"][:]
	lvl = qfcds["level"][:]
	qfc = qfcds["qfc_air"][:] * 1
	close(qfcds)
else lon,lat,lvl,qfc = calcqfc(regID)
end

# ╔═╡ e0251f16-5a91-11eb-1bc8-a70691df76e2
begin
	lds = NCDataset(datadir("reanalysis/era5-TRPx0.25-lsm-sfc.nc"))
	lsm = lds["lsm"][:]*1
	close(lds)
	
md"Loading Land-Sea Mask for ERA5 data ..."
end

# ╔═╡ 05f8a0be-5a92-11eb-0187-e19a897d9bf9
begin
	lqfc_SEA,sqfc_SEA = getmean([10,-10,135,90],qfc,lon,lat,length(lvl),lsm)
	lqfc_TRA,sqfc_TRA = getmean([10,-10,40,-10],qfc,lon,lat,length(lvl),lsm)
	lqfc_CRB,sqfc_CRB = getmean([25,15,-60,-90],qfc,lon,lat,length(lvl),lsm)
	lqfc_AMZ,sqfc_AMZ = getmean([10,-10,-45,-75],qfc,lon,lat,length(lvl),lsm)
	lqfc_DTP,sqfc_DTP = getmean([10,-10,360,0],qfc,lon,lat,length(lvl),lsm)
end

# ╔═╡ 46671d58-5a92-11eb-2a04-975238148a1b
begin
	pplt.close(); fqfc,aqfc = pplt.subplots(ncols=2,aspect=0.5,axwidth=1.5);
	
	lgd = Dict("frame"=>false,"ncols"=>1)
	
	aqfc[1].plot(lqfc_SEA*1e8,lvl,c="b")
	aqfc[1].plot(lqfc_TRA*1e8,lvl,c="r")
	aqfc[1].plot(lqfc_CRB*1e8,lvl,c="blue3")
	aqfc[1].plot(lqfc_AMZ*1e8,lvl,c="g")
	
	aqfc[2].plot(sqfc_DTP*1e8,lvl,c="k",label="DTP",legend="r")
	aqfc[2].plot(sqfc_SEA*1e8,lvl,c="b",label="SEA",legend="r",legend_kw=lgd)
	aqfc[2].plot(sqfc_TRA*1e8,lvl,c="r",label="TRA",legend="r")
	aqfc[2].plot(sqfc_CRB*1e8,lvl,c="blue3",label="CRB",legend="r")
	aqfc[2].plot(lqfc_AMZ*1e8*NaN,lvl,c="g",label="AMZ",legend="r")
	
	aqfc[1].format(
		xlim=(-12.5,12.5),ylim=(1000,50),urtitle="Land",
		ylabel="Pressure / hPa",
		xlabel=L"$-\nabla\cdot(q\vec{u})$ / 10$^{-8}$ kg kg$^{-1}$"
	)
	aqfc[2].format(xlim=(-4,4),ylim=(1000,50),urtitle="Ocean")
	
	fqfc.savefig(plotsdir("qfc.png"),transparent=false,dpi=200)
	load(plotsdir("qfc.png"))
end

# ╔═╡ 8774c3bc-5aa0-11eb-2dd9-d1f3daf9d815
md"We then proceed to replot this by height"

# ╔═╡ 06eebd8c-5aa1-11eb-12d9-7314d15a98c9
begin
	zds = NCDataset(datadir("reanalysis/era5-$(regID)x0.25-z_air.nc"))
	zlon = zds["longitude"][:]
	zlat = zds["latitude"][:]
	zlvl = zds["level"][8:end]
	zair = zds["z_air"][:,:,8:end] /1e3/9.81
	close(zds)
end

# ╔═╡ c4795296-5aa0-11eb-2038-1be17d8ce5cd
begin
	lzair_SEA,szair_SEA = getmean([10,-10,135,90],zair,zlon,zlat,length(zlvl),lsm)
	lzair_TRA,szair_TRA = getmean([10,-10,40,-10],zair,zlon,zlat,length(zlvl),lsm)
	lzair_CRB,szair_CRB = getmean([25,15,-60,-90],zair,zlon,zlat,length(zlvl),lsm)
	lzair_AMZ,szair_AMZ = getmean([10,-10,-45,-75],zair,zlon,zlat,length(zlvl),lsm)
	lzair_DTP,szair_DTP = getmean([10,-10,360,0],zair,zlon,zlat,length(zlvl),lsm)
end

# ╔═╡ f3bd1c4a-5aa0-11eb-049c-8d6783190671
begin
	pplt.close(); fzair,azair = pplt.subplots(ncols=2,aspect=0.5,axwidth=1.5);
	
	azair[1].plot(lqfc_SEA*1e8,lzair_SEA,c="b")
	azair[1].plot(lqfc_TRA*1e8,lzair_TRA,c="r")
	azair[1].plot(lqfc_CRB*1e8,lzair_CRB,c="blue3")
	azair[1].plot(lqfc_AMZ*1e8,lzair_AMZ,c="g")
	
	azair[2].plot(sqfc_DTP*1e8,szair_DTP,c="k",label="DTP",legend="r")
	azair[2].plot(sqfc_SEA*1e8,szair_SEA,c="b",label="SEA",legend="r",legend_kw=lgd)
	azair[2].plot(sqfc_TRA*1e8,szair_TRA,c="r",label="TRA",legend="r")
	azair[2].plot(sqfc_CRB*1e8,szair_CRB,c="blue3",label="CRB",legend="r")
	azair[2].plot(sqfc_AMZ*1e8*NaN,szair_AMZ/1e3,c="g",label="AMZ",legend="r")
	
	azair[1].format(
		xlim=(-12.5,12.5),ylim=(0.1,20),urtitle="Land",
		ylabel="Height / km",
		xlabel=L"$-\nabla\cdot(q\vec{u})$ / 10$^{-8}$ kg kg$^{-1}$"
	)
	azair[2].format(xlim=(-4,4),urtitle="Ocean",)#yscale="log")
	
	fzair.savefig(plotsdir("zair.png"),transparent=false,dpi=200)
	load(plotsdir("zair.png"))
end

# ╔═╡ 0f6c5ccc-5aa5-11eb-001f-75a4b37304d6
begin
	lsf = lsfinit(length(zlvl)); lsf[:,2] .= -999.0
	
	lsf[:,1] .= reverse(lzair_SEA)
	lsf[:,4] .= reverse(lqfc_SEA)
	lsfprint(projectdir("exp/lsf/SEA-land"),lsf,1009.32)
	
	lsf[:,1] .= reverse(szair_SEA)
	lsf[:,4] .= reverse(sqfc_SEA)
	lsfprint(projectdir("exp/lsf/SEA-ocean"),lsf,1009.32)
end

# ╔═╡ Cell order:
# ╟─582bcb26-5a22-11eb-3b52-09e3ec08b3f9
# ╟─36a4e004-5a26-11eb-1ba4-db83d18c7845
# ╟─3a8bbab2-5a26-11eb-20e0-d551c2a68858
# ╟─2b736522-5a33-11eb-0eea-473a5b8566f3
# ╠═28ffa404-5a33-11eb-3dc7-7d716d6cf06f
# ╠═9e01c7ea-5a34-11eb-1c41-8b991a399584
# ╠═87a146d4-5a29-11eb-1571-b9d9266ee54b
# ╟─abde634c-5a83-11eb-00d6-69b233e0c00f
# ╟─4927b4c2-5a8a-11eb-37e4-53efb9445b23
# ╠═795a822a-5a29-11eb-37bb-3b6d62cec6be
# ╠═e3c519fc-5a35-11eb-0d34-cd4ca8b50bf9
# ╟─e0251f16-5a91-11eb-1bc8-a70691df76e2
# ╠═05f8a0be-5a92-11eb-0187-e19a897d9bf9
# ╟─46671d58-5a92-11eb-2a04-975238148a1b
# ╟─8774c3bc-5aa0-11eb-2dd9-d1f3daf9d815
# ╠═06eebd8c-5aa1-11eb-12d9-7314d15a98c9
# ╠═c4795296-5aa0-11eb-2038-1be17d8ce5cd
# ╠═f3bd1c4a-5aa0-11eb-049c-8d6783190671
# ╟─0f6c5ccc-5aa5-11eb-001f-75a4b37304d6
