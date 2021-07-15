### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 270f908f-59ed-4a2c-b362-dc644f1220dd
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 2b8d1f00-c20b-4a46-a400-b8718755a719
begin
	@quickactivate "TroPrecLS"
	using ClimateERA
	using CLIMAParameters
	using CLIMAParameters.Planet
	using Dates
	using GeoRegions
	using Insolation
	using NCDatasets
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 82c17f18-d863-11eb-15c5-79aa6ea7b5b5
md"
# 10a. Insolation and Regional Averages

In this notebook, we just take a quick look at the seasonal cycle of skin temperature and Insolation, using `Insolation.jl`.
"

# ╔═╡ 6303f62c-6e13-40b7-bbec-2ad87975c4cf
md"
### A. Let's do some `Insolation.jl` stuff ...
"

# ╔═╡ 63bb4e71-36b6-4f4a-8e14-66c748a94779
begin
	struct EarthParameterSet <: AbstractEarthParameterSet end
	const param_set = EarthParameterSet()
end

# ╔═╡ 01415268-650e-43a0-9604-f6d10b165c23
begin
	ndy = 14610
	dys = collect(range(0,stop=14609,length=ndy))
	tlt = collect(-90:0.25:90); nlat = length(tlt)
	ins = zeros(ndy,nlat)
	for (j, latii) in enumerate(tlt), (i, dayii) in enumerate(dys)
		date = DateTime(1980,1,1) + Day(dayii)
		θ,dist = daily_zenith_angle(date,latii,EarthParameterSet())
		ins[i,j] = insolation(θ,dist,EarthParameterSet())
	end
end

# ╔═╡ b68aa1c3-acea-4575-864d-897e85a99653
md"
### B. What about regional insolation stuff?
"

# ╔═╡ fc592dc3-1fdd-4579-823d-c05a95171a49
begin
	dsl = NCDataset(datadir("reanalysis","era5-TRPx0.25-lsm-sfc.nc"))
	lon = dsl["longitude"][:]
	lat = dsl["latitude"][:]
	lsm = dsl["lsm"][:] * 1
	close(dsl)
	md"Loading land-sea mask data ..."
end

# ╔═╡ a37b6d52-6255-4c05-b9a3-317fbe247a93
function regioninsol(geo,lon,lat,lsm)
	
	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
    rins = zeros(nlon,nlat)
    rlsm = zeros(nlon,nlat)
    rwgt = zeros(nlon,nlat)
	oins = zeros(14610)
	lins = zeros(14610)

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for it = 0 : 14609
		date = DateTime(1980,1,1) + Day(it)
		for glat in 1 : nlat
			θ,dist = daily_zenith_angle(date,ggrd.glat[glat],EarthParameterSet())
			ins = insolation(θ,dist,EarthParameterSet())
			for glon in 1 : nlon
				rlsm[glon,glat] = lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
				rins[glon,glat] = ins * cosd.(ggrd.glat[glat])
				rwgt[glon,glat] = cosd.(ggrd.glat[glat])
			end
		end
		oins[it+1] = sum(rins[rlsm.<0.5]) ./ sum(rwgt[rlsm.<0.5])
		lins[it+1] = sum(rins[rlsm.>0.5]) ./ sum(rwgt[rlsm.>0.5])
	end
	
	return oins,lins
	
end

# ╔═╡ 9b2157e0-8b0f-4f8b-9e37-9a8f435527e9
oinstst,linstst = regioninsol(GeoRegion("SEA"),lon,lat,lsm)

# ╔═╡ eeb1b9b7-4fea-46f0-a483-c01315c0137f
begin
	pplt.close(); f,a = pplt.subplots(aspect=4,axwidth=6)
	
	a[1].plot((0:14609)/365.25,oinstst)
	a[1].plot((0:14609)/365.25,linstst)
	a[1].format(xlim=(0,4))
	
	f.savefig("testinsol.png",transparent=false,dpi=200)
	PNGFiles.load("testinsol.png")
end

# ╔═╡ Cell order:
# ╟─82c17f18-d863-11eb-15c5-79aa6ea7b5b5
# ╟─270f908f-59ed-4a2c-b362-dc644f1220dd
# ╟─2b8d1f00-c20b-4a46-a400-b8718755a719
# ╟─6303f62c-6e13-40b7-bbec-2ad87975c4cf
# ╠═63bb4e71-36b6-4f4a-8e14-66c748a94779
# ╠═01415268-650e-43a0-9604-f6d10b165c23
# ╟─b68aa1c3-acea-4575-864d-897e85a99653
# ╠═fc592dc3-1fdd-4579-823d-c05a95171a49
# ╠═a37b6d52-6255-4c05-b9a3-317fbe247a93
# ╠═9b2157e0-8b0f-4f8b-9e37-9a8f435527e9
# ╠═eeb1b9b7-4fea-46f0-a483-c01315c0137f
