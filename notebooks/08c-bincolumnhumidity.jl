### A Pluto.jl notebook ###
# v0.19.23

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8f2172e0-5b81-420a-8301-dbb8ed0c290a
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ bb51d00a-ee75-4bd8-9744-381b87c6441b
begin
	@quickactivate "TroPrecLS"
	using Dierckx
	using ERA5Reanalysis
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
	include(srcdir("samlsf.jl"))
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ e0fee7ac-fc9b-11ec-0741-9f810d7bcd4c
md"
# 08b. Preliminary Data Exploration

A preliminary exploration of the 2D output data for column relative humidity and outgoing longwave radiation
"

# ╔═╡ 2f1ee85d-cef1-4fd2-ba06-7c5ec1a3d4a5
md"
### A. Checking the self-aggregated nature of our Spinup
"

# ╔═╡ cf2e61e6-b46f-4b77-960a-e1d8a69aeefe
md"Day $(@bind day PlutoUI.Slider(40:(1/24):50,default=50,show_value=true))"

# ╔═╡ 6c4ec7a4-bac6-4033-81c0-26af845cc11f
begin
	ds  = NCDataset(datadir("Aggregation","OUT_2D","RCE_TroPrecLS-Aggregation-control_64_0000$(@sprintf("%06d",day*2880)).2Dbin_1.nc"))
	x   = ds["x"][:] / 1000
	y   = ds["y"][:] / 1000
	pwv = ds["PW"][:,:,1]
	swv = ds["SWVP"][:,:,1]
	olr = ds["LWNT"][:,:,1]
	csf = pwv ./ swv * 100
	# olr = csf
	close(ds)
end

# ╔═╡ 2244daac-8c12-4111-a883-94637d601c3a
begin
	pplt.close(); fig,axs = pplt.subplots(nrows=4,aspect=8,axwidth=6,hspace=1.25)

	lvls = 90 : 10 : 270
	cmp = "Blues"
	
	c = axs[1].pcolormesh(cmap=cmp,levels=lvls,x[1:512],y,olr[1:512,:]',extend="both")
	axs[2].pcolormesh(cmap=cmp,levels=lvls,x[1:512],y,olr[513:1024,:]' ,extend="both")
	axs[3].pcolormesh(cmap=cmp,levels=lvls,x[1:512],y,olr[1025:1536,:]',extend="both")
	axs[4].pcolormesh(cmap=cmp,levels=lvls,x[1:512],y,olr[1537:2048,:]',extend="both")

	for ax in axs
		ax.format(xlocator=0:128:1024,ylocator=0:64:128,ylim=(0,128),xlim=(0,1024))
	end
	fig.colorbar(c,length=0.75)
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ c163e44b-fd33-4389-9f44-77cef9c75b55
md"Create Animation? $(@bind isanim PlutoUI.Slider(0:1))"

# ╔═╡ 0560a557-70a8-4a69-95f1-39f73921b146
if isone(isanim)

	for idy = 40 : (1/24) : 50

		ids = NCDataset(datadir("Aggregation","OUT_2D","RCE_TroPrecLS-Aggregation-control_64_0000$(@sprintf("%06d",idy*2880)).2Dbin_1.nc"))
		ilw = ids["LWNT"][:,:,1]
		close(ids)

		pplt.close(); f2,a2 = pplt.subplots(nrows=4,aspect=8,axwidth=6)
	
		c = a2[1].pcolormesh(cmap="Blues",levels=120:15:270,x[1:512],y,ilw[1:512,:]',extend="both")
		a2[2].pcolormesh(cmap="Blues",levels=120:15:270,x[1:512],y,ilw[513:1024,:]' ,extend="both")
		a2[3].pcolormesh(cmap="Blues",levels=120:15:270,x[1:512],y,ilw[1025:1536,:]',extend="both")
		a2[4].pcolormesh(cmap="Blues",levels=120:15:270,x[1:512],y,ilw[1537:2048,:]',extend="both")
	
		for ax in a2
			ax.format(xlocator=0:128:1024,ylocator=0:32:128,ylim=(0,128),xlim=(0,1024))
		end
		f2.colorbar(c,length=0.75)
	
		f2.savefig(
			plotsdir("noshear-anim","$(@sprintf("%06d",idy*2880)).png"),
			transparent=false,dpi=150
		)
		
	end
	
end

# ╔═╡ 612f386a-6c7a-429d-9622-d346bc3cc798
md"
### B. Doing some spatial smoothing (0.25$\degree$, or ~28 km)
"

# ╔═╡ 5954b93e-3d1a-4299-99c4-00e5c7bc6286
begin
	pplt.close(); f3,a3 = pplt.subplots(nrows=4,aspect=8,axwidth=6)

	lvl2 = 10 : 5 : 90

	n = 14
	csftem = zeros(length(x),length(y))
	csfnew = zeros(length(x),length(y))
	for xshift = 1 : n, yshift = 1 : n
		circshift!(csftem,csf,(xshift-1,yshift-1))
		for iy = 1 : length(y), ix = 1 : length(x)
			csfnew[ix,iy] += csftem[ix,iy]
		end
	end
	csfnew /= n^2
	
	c3 = a3[1].pcolormesh(cmap=cmp,levels=lvl2,x[1:512],y,csfnew[1:512,:]',extend="both")
	a3[2].pcolormesh(cmap=cmp,levels=lvl2,x[1:512],y,csfnew[513:1024,:]' ,extend="both")
	a3[3].pcolormesh(cmap=cmp,levels=lvl2,x[1:512],y,csfnew[1025:1536,:]',extend="both")
	a3[4].pcolormesh(cmap=cmp,levels=lvl2,x[1:512],y,csfnew[1537:2048,:]',extend="both")

	for ax in a3
		ax.format(xlocator=0:128:1024,ylocator=0:32:128,ylim=(0,128),xlim=(0,1024))
	end
	f3.colorbar(c3,length=0.75)
	
	f3.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╠═e0fee7ac-fc9b-11ec-0741-9f810d7bcd4c
# ╟─8f2172e0-5b81-420a-8301-dbb8ed0c290a
# ╟─bb51d00a-ee75-4bd8-9744-381b87c6441b
# ╟─2f1ee85d-cef1-4fd2-ba06-7c5ec1a3d4a5
# ╟─cf2e61e6-b46f-4b77-960a-e1d8a69aeefe
# ╟─6c4ec7a4-bac6-4033-81c0-26af845cc11f
# ╟─2244daac-8c12-4111-a883-94637d601c3a
# ╟─c163e44b-fd33-4389-9f44-77cef9c75b55
# ╟─0560a557-70a8-4a69-95f1-39f73921b146
# ╟─612f386a-6c7a-429d-9622-d346bc3cc798
# ╟─5954b93e-3d1a-4299-99c4-00e5c7bc6286
