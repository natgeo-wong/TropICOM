### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ abdaeba6-56a5-11eb-222b-a12819fce07c
begin
	using Pkg; Pkg.activate()
	using DrWatson

md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ ae7e6c26-56a5-11eb-2444-694b368de12c
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using ImageFiltering
	using NASAPrecipitation
	using NCDatasets
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ f838a206-56a1-11eb-35dc-f1489317cd8a
md"
# 2a. GPM land-sea mask data

For GPM IMERGv6 precipitation, NASA has provided a land-sea mask denoting the percentage of the grid-box that is considered to be water.  In this notebook, we do a quick exploration of this dataset and extract the data for certain regions.  We also convert the land-sea mask from percentage of water, into the 0-1 range similar to that found in ECMWF reanalysis land-sea mask datasets.
"

# ╔═╡ a907ca40-56dd-11eb-07e6-dd780174001e
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ b28b2cbe-56a7-11eb-3a3c-09b80ad4e0f2
md"
### A. Loading and presenting GPM IMERG land-sea mask

We start by loading and plotting the GPM IMERG land-sea mask ...
"

# ╔═╡ c43c087a-56a7-11eb-2c12-7db44d0a7eae
lsd = getIMERGlsd(GeoRegion("GLB"),path=datadir("imergmask"))

# ╔═╡ ed63a2e2-56a7-11eb-10f2-2b291133bb4e
begin
	asp1 = (maximum(lsd.lon)-minimum(lsd.lon))/(maximum(lsd.lat)-minimum(lsd.lat))
	pplt.close(); f,axs = pplt.subplots(aspect=6,axwidth=6)

	c = axs[1].pcolormesh(
		lsd.lon,lsd.lat,lsd.lsm',levels=0:10:100,
		cmap="Delta_r",extend="both"
	)
	axs[1].colorbar(c,loc="r")

	for ax in axs
		ax.format(xlim=(0,360),ylim=(-30,30),xlocator=-0:60:360)
	end

	f.savefig(plotsdir("02a-GPMlsm-TRP.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("02a-GPMlsm-TRP.png"))
end

# ╔═╡ b5841fc6-56ab-11eb-0c0d-cdf45571f6d7
md"
### B. Extracting and Regridding the ETOPO1 data

Here, we define the functions that are used to extract and regrid the land-sea mask, using the `GeoRegions.jl` package.  The land-sea mask, is recalculated from percentage water, into fraction which is land, similar to ECMWF land-sea mask data.
"

# ╔═╡ e82fb766-56f1-11eb-34aa-f35e027df299
begin
	lsm01 = lsd.lsm / 100
	lsm01 = 1 .- lsm01
	lsm01[lsm01.>0.1] .= 1
	lsm01[lsm01.<0.1] .= 0

md"Binning land-sea mask into 0s and 1s"
end

# ╔═╡ 3443392c-56ad-11eb-3a44-137601990567
md"Using the functions `regiongridvec` and `regionextractgrid`, we are able to extract the land-sea mask for a given region."

# ╔═╡ 0b25e87e-56f0-11eb-1cb9-4574948fd72a
geo = GeoRegion("TRP")

# ╔═╡ 0b9d45f8-56dd-11eb-3928-2b11e9d7bdf7
begin
	ggrd = RegionGrid(geo,lsd.lon,lsd.lat)
	rlsm = extractGrid(lsm01,ggrd)
	md"Extracting Land-Sea Mask for the region ..."
end

# ╔═╡ 4536cd18-56dd-11eb-0a9b-a1530a862491
begin
	asp = (maximum(ggrd.lon)-minimum(ggrd.lon))/(maximum(ggrd.lat)-minimum(ggrd.lat))
	pplt.close(); freg,areg = pplt.subplots(aspect=asp,axwidth=asp)

	creg = areg[1].pcolormesh(
		ggrd.lon,ggrd.lat,rlsm',levels=(0:10)/10,
		cmap="Delta",extend="both"
	)
	areg[1].colorbar(creg,loc="r")

	for ax in areg
		ax.format(
			xlim=(minimum(ggrd.lon),maximum(ggrd.lon)),
			ylim=(minimum(ggrd.lat),maximum(ggrd.lat)),
			xlocator=0:60:360,
		)
	end

	freg.savefig(plotsdir("02a-GPMlsmreg.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("02a-GPMlsmreg.png"))
end

# ╔═╡ 13a9227b-bbec-4987-8436-fca7a9f047ea
md"
### C. Smoothing the Land-Sea Mask?
"

# ╔═╡ 4d386a7d-503c-4fa1-ae59-e79e9690580b
function filterlsm(olsm;iterations=1,smooth=1)
	
	it = 0
	nlsm = deepcopy(olsm)
	while it < iterations
		nlsm = log10.(imfilter(10. .^nlsm, Kernel.gaussian(smooth),"circular"));
		nlsm = (nlsm.+olsm)/2
		it += 1
	end
	
	nlsm[nlsm.<0] .= minimum(nlsm[nlsm.>0])
	
	return nlsm
	
end

# ╔═╡ 9f4165e6-1203-4626-8f4b-3d1d278cfb4a
begin
	nlsm1d0 = filterlsm(lsm01,iterations=10,smooth=2.5)
	nlsm1d5 = filterlsm(lsm01,iterations=10,smooth=3.75)
	nlsm2d0 = filterlsm(lsm01,iterations=10,smooth=5)
	md"Performing gaussian filtering/smoothing on land-sea mask ..."
end

# ╔═╡ 7f72bf8a-998f-4060-bbfa-afc830afc423
begin
	N,S,E,W = geo.N,geo.S,geo.E,geo.W
	flsm1d0 = extractGrid(nlsm1d0,ggrd)
	flsm1d5 = extractGrid(nlsm1d5,ggrd)
	flsm2d0 = extractGrid(nlsm2d0,ggrd)
	md"Extracting information for region ..."
end

# ╔═╡ 9870db8d-ea2b-419e-ad2d-0c5f7e268ce0
begin
	pplt.close(); asp2 = (E-W+2)/(N-S+2)
	if asp2 >= 2.9
		fflt,aflt = pplt.subplots(nrows=4,axwidth=asp2*1.2,aspect=asp2)
	else
		fflt,aflt = pplt.subplots(nrows=2,ncols=2,axwidth=asp2*1.5,aspect=asp2)
	end
	
	lvls = vcat(10. .^(-4:-1),0.2,0.5,0.8,0.9,0.99,0.999,0.9999)

	cflt = aflt[1].pcolormesh(
		ggrd.lon,ggrd.lat,rlsm',
		levels=lvls,cmap="delta",extend="both"
	)
	aflt[1].format(urtitle="Raw")
	
	aflt[2].pcolormesh(
		ggrd.lon,ggrd.lat,flsm1d0',
		levels=lvls,cmap="delta",extend="both"
	)
	aflt[2].format(urtitle="Smooth = 2.5")
	
	aflt[3].pcolormesh(
		ggrd.lon,ggrd.lat,flsm1d5',
		levels=lvls,cmap="delta",extend="both"
	)
	aflt[3].format(urtitle="Smooth = 3.75")
	
	aflt[4].pcolormesh(
		ggrd.lon,ggrd.lat,flsm2d0',
		levels=lvls,cmap="delta",extend="both"
	)
	aflt[4].format(urtitle="Smooth = 5")

	for ax in aflt
		ax.plot(x,y,c="k",lw=0.5)
		ax.format(
			xlim=(ggrd.lon[1].-1,ggrd.lon[end].+1),ylim=(S-1,N+1),
			ylabel=L"Latitude / $\degree$",xlabel=L"Longitude / $\degree$",
			suptitle="Land-Sea Mask",grid=true
		)
	end

	if asp2 >=2.9
		fflt.colorbar(cflt,loc="r",length=0.4)
	else
		fflt.colorbar(cflt,loc="r",length=0.8)
	end
	fflt.savefig(plotsdir("02a-fgpmlsm_$(geo.regID).png"),transparent=false,dpi=200)
	load(plotsdir("02a-fgpmlsm_$(geo.regID).png"))
end

# ╔═╡ Cell order:
# ╟─f838a206-56a1-11eb-35dc-f1489317cd8a
# ╟─abdaeba6-56a5-11eb-222b-a12819fce07c
# ╟─ae7e6c26-56a5-11eb-2444-694b368de12c
# ╟─a907ca40-56dd-11eb-07e6-dd780174001e
# ╟─b28b2cbe-56a7-11eb-3a3c-09b80ad4e0f2
# ╟─c43c087a-56a7-11eb-2c12-7db44d0a7eae
# ╟─ed63a2e2-56a7-11eb-10f2-2b291133bb4e
# ╟─b5841fc6-56ab-11eb-0c0d-cdf45571f6d7
# ╟─e82fb766-56f1-11eb-34aa-f35e027df299
# ╟─3443392c-56ad-11eb-3a44-137601990567
# ╠═0b25e87e-56f0-11eb-1cb9-4574948fd72a
# ╟─0b9d45f8-56dd-11eb-3928-2b11e9d7bdf7
# ╟─4536cd18-56dd-11eb-0a9b-a1530a862491
# ╟─13a9227b-bbec-4987-8436-fca7a9f047ea
# ╠═4d386a7d-503c-4fa1-ae59-e79e9690580b
# ╟─9f4165e6-1203-4626-8f4b-3d1d278cfb4a
# ╟─7f72bf8a-998f-4060-bbfa-afc830afc423
# ╟─9870db8d-ea2b-419e-ad2d-0c5f7e268ce0
