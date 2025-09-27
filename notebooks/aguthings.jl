### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 6f00b8fc-530c-11eb-2242-99d8544f6e14
begin
	using Pkg; Pkg.activate()
	using DrWatson

md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 8f30c56c-530c-11eb-2782-33f3c4ed9e89
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using Interpolations
	using NASAPrecipitation
	using NCDatasets
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 90fffbc8-524d-11eb-232a-1bada28d5505
md"
# 2. Investigation of Precipitation Data

We first investigate GPM IMERGv6 precipitation data from 2001-2019 over the domains listed in the notebook `01-domains`.  The data was retrieved from the website of NASA's Precipitation Measurement Mission (PMM), and is found in 0.1º x 0.1º horizontal resolution every half-hour.

The data was downloaded using the `ClimateSatellite.jl` package, which downloads the data in monthly batches in the date range specified.  More information can be found in the package documentation, and instructions on how to use the package to download the data can be found there as well.
"

# ╔═╡ d82366b0-53b1-11eb-26c1-ff1bb6ccb027
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ c1fafcba-530f-11eb-1cc2-67d10a3fa606
md"
### A. Characteristics of the Diurnal Cycle of Precipitation

We wish to find the following characteristics of the diurnal cycle of precipitation rate:
* The mean $\mu$ rate of precipitation
* The amplitude $A$ of the diurnal cycle (max-min)/2 of precipitation
* The hour $\theta$ at which the rate of precipitation is a maximum
"

# ╔═╡ 3565af3c-5311-11eb-34c4-2d228b05b17c
md"
GPS IMERGv6 output time is defined according to UTC.  Therefore, for every grid point, we must also correct to the local time before fitting the data to the model.
"

# ╔═╡ 577a59d8-5311-11eb-2e2e-53bbeff12648
longitude2timeshift(longitude::Real) = longitude / 180 * 12

# ╔═╡ a6a688ca-53ab-11eb-2776-b5380ffb26c1
function gpm2model(data,longitude)

	nlon,nlat,nt = size(data)
	θ = zeros(nlon,nlat)
	A = zeros(nlon,nlat)
	μ = zeros(nlon,nlat)

	p0 = [0,0.5,0]; t = ((1:(nt+1)).-0.5)/2
	ts = longitude2timeshift.(longitude)
	it = 0:0.1:24; it = it[1:(end-1)]; nit = length(it)

	idat = zeros(nt+1); var  = zeros(nit)

	for ilat = 1 : nlat, ilon = 1 : nlon

        tl = t .+ ts[ilon]

        idat[1:nt] = data[ilon,ilat,:]
		idat[end]  = data[ilon,ilat,1]

		itp = interpolate(idat,BSpline(Cubic(Periodic(OnGrid()))))
        stp = scale(itp,tl)
        etp = extrapolate(stp,Periodic())
		var[:] = etp[it]
		var[var.<0] .= 0

		μ[ilon,ilat] = mean(@view data[ilon,ilat,:])
		A[ilon,ilat] = (maximum(idat) .- minimum(idat))/2
		θ[ilon,ilat] = argmax(var) / nit * 24

    end

	return μ,A,θ

end

# ╔═╡ b4eec397-798d-4332-9365-c5ad374d9b4f
md"
### B. Defining GeoRegions
"

# ╔═╡ aa05317e-530b-11eb-2ec1-93aff65659dd
md"
### B. Retrieving GPM Precipitation Data
"

# ╔═╡ bb90be66-554c-11eb-05de-a5db553ad4b1
md"
We fit the model to GPM IMERGv6 data from 2001 to 2018 for each individual spatial point (i.e. at 0.1º resolution).  Our results are plotted below:
"

# ╔═╡ 103f85e8-530c-11eb-047d-a537aa60075d
function retrievegpm()

    ds  = NCDataset(datadir("compiled","imergfinalhh-TRP-compiled.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	var = ds["precipitation"][:,:,:] * 3600
	close(ds)

	μ,A,θ = gpm2model(var,lon)

    return lon,lat,μ,A,θ

end

# ╔═╡ e8141e20-53af-11eb-1a23-81d34293c5eb
begin
	lon,lat,μ,A,θ = retrievegpm()
	md"Extracting diurnal information for GPM IMERGv6 precipitation ..."
end

# ╔═╡ 271c15e6-9426-11eb-3c64-912620b1a8ce
begin
	geo1 = RectRegion("","GLB","",[7,-4,119.5,108.5],save=false)
	geo2 = RectRegion("","GLB","",[6,-6,107,95],save=false)
	md"Defining region coordinates ..."
end

# ╔═╡ ea7f0956-575b-11eb-3e3f-a1ba3e08b771
begin
	ggrd1 = RegionGrid(geo1,lon,lat)
	θ1 = extractGrid(θ,ggrd1)
	md"Extracting information for region ..."
end

# ╔═╡ c3897f33-c7db-4b72-aa2a-a659428ef307
begin
	ggrd2 = RegionGrid(geo2,lon,lat)
	θ2 = extractGrid(θ,ggrd2)
	md"Extracting information for region ..."
end

# ╔═╡ 5714c13c-575c-11eb-06d4-838b4e8dbcd7
begin
	pplt.close(); fig,axs = pplt.subplots(nrows=2,sharex=0,sharey=0)
	
	c = axs[1].pcolormesh(ggrd1.lon,ggrd1.lat,θ1',cmap="romaO",levels=0:24)
	axs[2].pcolormesh(ggrd2.lon,ggrd2.lat,θ2',cmap="romaO",levels=0:24)

	axs[1].format(xlim=(geo1.W,geo1.E),ylim=(geo1.S,geo1.N))
	axs[2].format(xlim=(geo2.W,geo2.E),ylim=(geo2.S,geo2.N))

	for ax in axs
		ax.plot(x,y,lw=0.5,c="k")
	end

	fig.colorbar(c,loc="b",locator=0:3:24,label="Hour of Day")
	fig.savefig("test.png",transparent=true,dpi=800)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─90fffbc8-524d-11eb-232a-1bada28d5505
# ╟─6f00b8fc-530c-11eb-2242-99d8544f6e14
# ╟─8f30c56c-530c-11eb-2782-33f3c4ed9e89
# ╟─d82366b0-53b1-11eb-26c1-ff1bb6ccb027
# ╟─c1fafcba-530f-11eb-1cc2-67d10a3fa606
# ╟─3565af3c-5311-11eb-34c4-2d228b05b17c
# ╠═577a59d8-5311-11eb-2e2e-53bbeff12648
# ╠═a6a688ca-53ab-11eb-2776-b5380ffb26c1
# ╟─b4eec397-798d-4332-9365-c5ad374d9b4f
# ╟─aa05317e-530b-11eb-2ec1-93aff65659dd
# ╟─bb90be66-554c-11eb-05de-a5db553ad4b1
# ╠═103f85e8-530c-11eb-047d-a537aa60075d
# ╟─e8141e20-53af-11eb-1a23-81d34293c5eb
# ╠═271c15e6-9426-11eb-3c64-912620b1a8ce
# ╟─ea7f0956-575b-11eb-3e3f-a1ba3e08b771
# ╟─c3897f33-c7db-4b72-aa2a-a659428ef307
# ╠═5714c13c-575c-11eb-06d4-838b4e8dbcd7
