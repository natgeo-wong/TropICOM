### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 8fea5fc0-4f97-11eb-295b-f98a3c89c4cf
begin
	using DrWatson
	using Pkg
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 9b9c0b36-4f97-11eb-32a1-0d1f15b34aa5
begin
	@quickactivate "TroPrecLS"
	Pkg.instantiate()
	using DelimitedFiles
	using GeoRegions
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase
	
	using ImageShow, FileIO
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot");
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ b1798d06-4f94-11eb-0132-4370fe8d3a0f
md"
# X. Finding representative Skin Temperatures over the Maritime Continent

In this notebook, I plotted the spatial distribution of average skin temperature distribution within the region of the western Maritime Continent in order to find a representative skin-temperature that we can use for our baseline RCE simulations.  I downloaded data from monthly-averaged ERA5 reanalysis dataset for the [TropicalRCE](https://github.com/natgeo-wong/TropicalRCE) experiment, which I duplicated for this project ([TroPrecLS]()).

Instructions on how to download the data for yourself can be found in the notebook *'XXX. Notebook'*
"

# ╔═╡ 032b77fc-4f99-11eb-2eb9-af412c4f1d25
md"
### 1. Function Setup

Drawing upon the region-extraction capabilities of `GeoRegions.jl`, we define the functions `bindatalnd()` and `bindatasea()` that bins data within a user-defined rectilinear region, and according to the land-sea mask provided, bins data over land and ocean respectively."

# ╔═╡ aaf0ea76-4f99-11eb-3e92-63bd3414e7f7
begin
	function bindatalnd(coords,bins,var,lon,lat,lsm)
	
	    tlon,tlat,rinfo = regiongridvec(coords,lon,lat);
	    rvar = regionextractgrid(var,rinfo)
	    rlsm = regionextractgrid(lsm,rinfo)
	    rwgt = ones(size(rlsm)) .* cosd.(reshape(tlat,1,:))
	
	    lvar = rvar[rlsm.>0.5];
		lvar = lvar[.!ismissing.(lvar)]; lvar = lvar[.!isnan.(lvar)]
	    lbin = fit(Histogram,lvar,bins).weights;
		lbin = lbin ./ sum(lbin) * (length(bins) - 1)
	
	    rvar = rvar .* cosd.(reshape(tlat,1,:))
	    lvar = rvar[rlsm.>0.5]; 
		lvar = lvar[.!ismissing.(lvar)]; lvar = lvar[.!isnan.(lvar)]
	    lvar = lvar / mean(rwgt[rlsm.>0.5])
	
	    return lbin,mean(lvar)
	
	end
	
	function bindatasea(coords,bins,var,lon,lat,lsm)
	
	    tlon,tlat,rinfo = regiongridvec(coords,lon,lat);
	    rvar = regionextractgrid(var,rinfo)
	    rlsm = regionextractgrid(lsm,rinfo)
	    rwgt = ones(size(rlsm)) .* cosd.(reshape(tlat,1,:))
	
	    svar = rvar[rlsm.<0.5];
		svar = svar[.!ismissing.(svar)]; svar = svar[.!isnan.(svar)]
	    sbin = fit(Histogram,svar,bins).weights;
		sbin = sbin ./ sum(sbin) * (length(bins) - 1)
	
	    rvar = rvar .* cosd.(reshape(tlat,1,:))
	    svar = rvar[rlsm.<0.5];
		svar = svar[.!ismissing.(svar)]; svar = svar[.!isnan.(svar)]
	    svar = svar / mean(rwgt[rlsm.<0.5])
	
	    return sbin,mean(svar)
	
	end
end

# ╔═╡ 9531925e-4f9d-11eb-0e09-853a81c5f83e
md"N : $(@bind N PlutoUI.Slider(10:5:20))"

# ╔═╡ b65293a2-4f9d-11eb-2cce-05f276357434
md"S : $(@bind S PlutoUI.Slider(-20:5:-5))"

# ╔═╡ b8c1f15a-4f9d-11eb-1089-95cdd781aa0f
md"W : $(@bind W PlutoUI.Slider(80:5:100))"

# ╔═╡ bc01b6de-4f9d-11eb-3094-1b24d8af464b
md"E : $(@bind E PlutoUI.Slider(135:5:180))"

# ╔═╡ 485ba5b2-4f9d-11eb-3860-c18a35e3f1c2
md"### 2. Region of Interest

*Note: use the sliders to adjust the bounding box for the region of interest*

We define our region of interest to be [$N, $S, $W, $E], and show the geographical domain in the figure below:"

# ╔═╡ ec2265ca-4f9d-11eb-0a9a-7b4b9aaa2763
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
	aratio = (E-W)/(N-S)
	
	pplt.close(); f1,axs1 = pplt.subplots(aspect=aratio,axwidth=4,sharey=0);
	axs1[1].plot(x,y,c="k",lw=0.2)
	axs1[1].format(
		xlim=(W,E),xlocator=80:10:180,
    	ylim=(S,N),ylocator=-30:10:30
	)
	
	f1.savefig("domain.png",transparent=false,dpi=200)
	load("domain.png")
end

# ╔═╡ 0b9f29f2-4f9d-11eb-0680-8b2dbabf879c
md"### 3. Plotting and Analysis

First, we load the following datsets: (1) Skin surface temperature data, and (2) Land-sea mask"

# ╔═╡ ca513f5c-4f98-11eb-2f12-87c79bfb38d1
begin
	ds = NCDataset(datadir("reanalysis/era5-TRPx0.25-skt-sfc.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	skt = ds["skt"][:]*1
	close(ds)
	
	ds = NCDataset(datadir("reanalysis/era5-TRPx0.25-lsm-sfc.nc"))
	lsm = ds["lsm"][:]*1
	close(ds)
end

# ╔═╡ 48403866-4f9c-11eb-1de4-bb0223eae1bd
md"lower bound (sea) : $(@bind smin PlutoUI.Slider(295:300))"

# ╔═╡ 534e72f4-4f9c-11eb-0cba-2501d13478ce
md"upper bound (sea) : $(@bind smax PlutoUI.Slider(303:0.5:305))"

# ╔═╡ 8efd83be-4f9c-11eb-1211-0fcd19ef5976
md"lower bound (land) : $(@bind lmin PlutoUI.Slider(290:0.5:300))"

# ╔═╡ 94d00a08-4f9c-11eb-372c-214e62a0687e
md"upper bound (land) : $(@bind lmax PlutoUI.Slider(300:305))"

# ╔═╡ 85bd3c62-4f9b-11eb-3523-6f6a181001c8
md"
Next, we define our bin ranges over Land ( $lmin, $lmax ) and Sea: ( $smin, $smax )

*Note: use the sliders to adjust the range of the bins*"

# ╔═╡ ca3e0c4c-4f9e-11eb-0c7b-3b6875674f64
md"We then proceed to bin the data and plot it"

# ╔═╡ de0a8f18-4f9e-11eb-2cfb-35a503c455e9
begin
	sbinvec = smin:0.05:smax; psbin = sbinvec[2:end] .- 0.025
	lbinvec = lmin:0.05:lmax; plbin = lbinvec[2:end] .- 0.025
	
	sbin,savg = bindatasea([N,S,E,W],sbinvec,skt,lon,lat,lsm)
	lbin,lavg = bindatalnd([N,S,E,W],lbinvec,skt,lon,lat,lsm)
	
	arr = [[1,1],[1,1],[2,3]]
	pplt.close(); f2,axs2 = pplt.subplots(arr,aspect=aratio,axwidth=4,sharey=3);
	axs2[1].plot(x,y,c="k",lw=0.2)
	c = axs2[1].contourf(lon,lat,skt',levels=298:0.5:303,extend="both")
	axs2[1].contour(lon,lat,skt',levels=301:303,colors="k",lw=0.5,linestyle="--")
	axs2[1].format(
		xlim=(W,E),xlocator=80:10:180,
    	ylim=(S,N),ylocator=-30:10:30
	)
	f2.colorbar(c,loc="r")
	
	axs2[2].plot([minimum(lbinvec),maximum(lbinvec)],[2,2],c="gray6",lw=0.5)
	axs2[3].plot([minimum(sbinvec),maximum(sbinvec)],[2,2],c="gray6",lw=0.5)
	
	axs2[2].plot(plbin,lbin,c="k",lw=0.5); axs2[2].plot([1,1]*lavg,[0.1,50],c="k")
	axs2[3].plot(psbin,sbin,c="k",lw=0.5); axs2[3].plot([1,1]*savg,[0.1,50],c="k")
	
	lavgs = @sprintf("%0.2f",lavg)
	savgs = @sprintf("%0.2f",savg)
	
	axs2[2].format(
        xlim=(minimum(lbinvec),maximum(lbinvec)),
		ylim=(0.1,10),yscale="log",ylabel="Density",
		ultitle="Land: $lavgs K"
    )
	
	axs2[3].format(
        xlim=(minimum(sbinvec),maximum(sbinvec)),
		ylim=(0.1,10),yscale="log",
		ultitle="Sea: $savgs K"
    )
	
	f2.savefig("era5skt.png",transparent=false,dpi=200)
	load("era5skt.png")
end

# ╔═╡ 94884a3c-4fa2-11eb-340a-15d16998969d
md"
### 4. Saving the figures
"

# ╔═╡ cb733eb4-4fa2-11eb-011b-83632ab01544
md"Save Plot? $(@bind dosave PlutoUI.Slider(0:1))"

# ╔═╡ f6eddeb2-4fa2-11eb-0389-bb88c0547f2c
if isone(dosave)
	if !isdir(plotsdir()); mkpath(plotsdir()) end
	f1.savefig(plotsdir("domain.png"),transparent=false,dpi=200)
	f2.savefig(plotsdir("era5skt.png"),transparent=false,dpi=200)

md"We have decided to save the plots"
else
md"We are not saving the plots yet"
end

# ╔═╡ Cell order:
# ╟─b1798d06-4f94-11eb-0132-4370fe8d3a0f
# ╟─8fea5fc0-4f97-11eb-295b-f98a3c89c4cf
# ╟─9b9c0b36-4f97-11eb-32a1-0d1f15b34aa5
# ╟─032b77fc-4f99-11eb-2eb9-af412c4f1d25
# ╠═aaf0ea76-4f99-11eb-3e92-63bd3414e7f7
# ╟─485ba5b2-4f9d-11eb-3860-c18a35e3f1c2
# ╟─9531925e-4f9d-11eb-0e09-853a81c5f83e
# ╟─b65293a2-4f9d-11eb-2cce-05f276357434
# ╟─b8c1f15a-4f9d-11eb-1089-95cdd781aa0f
# ╟─bc01b6de-4f9d-11eb-3094-1b24d8af464b
# ╟─ec2265ca-4f9d-11eb-0a9a-7b4b9aaa2763
# ╟─0b9f29f2-4f9d-11eb-0680-8b2dbabf879c
# ╠═ca513f5c-4f98-11eb-2f12-87c79bfb38d1
# ╟─85bd3c62-4f9b-11eb-3523-6f6a181001c8
# ╟─48403866-4f9c-11eb-1de4-bb0223eae1bd
# ╟─534e72f4-4f9c-11eb-0cba-2501d13478ce
# ╟─8efd83be-4f9c-11eb-1211-0fcd19ef5976
# ╟─94d00a08-4f9c-11eb-372c-214e62a0687e
# ╟─ca3e0c4c-4f9e-11eb-0c7b-3b6875674f64
# ╟─de0a8f18-4f9e-11eb-2cfb-35a503c455e9
# ╟─94884a3c-4fa2-11eb-340a-15d16998969d
# ╟─cb733eb4-4fa2-11eb-011b-83632ab01544
# ╟─f6eddeb2-4fa2-11eb-0389-bb88c0547f2c
