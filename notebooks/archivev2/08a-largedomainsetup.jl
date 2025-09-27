### A Pluto.jl notebook ###
# v0.19.22

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
# 08a. Setting up Large Domains

We have setup a large-domain simulation such that convective self-aggregation occurs. And then we set up zonal wind-shear so that we can see how those spots of convective self-aggregation move and evolve over time. This notebook here is simply to:

1. Double-check that self-aggregation occurs after 150 days of model run
2. Find the vertical temperature and pressure profiles
3. Superimpose the tropical wind shear from ERA5 upon it
"

# ╔═╡ 2f1ee85d-cef1-4fd2-ba06-7c5ec1a3d4a5
md"
### A. Checking the self-aggregated nature of our Spinup
"

# ╔═╡ cf2e61e6-b46f-4b77-960a-e1d8a69aeefe
md"Day $(@bind day PlutoUI.Slider(100:150,default=150,show_value=true))"

# ╔═╡ 6c4ec7a4-bac6-4033-81c0-26af845cc11f
begin
	ds  = NCDataset(datadir("Aggregation","OUT_2D","RCE_TroPrecLS-Aggregation-spinup_64_0000$(@sprintf("%06d",day*2880)).2Dbin_1.nc"))
	x   = ds["x"][:] / 1000
	y   = ds["y"][:] / 1000
	olr = ds["LWNT"][:,:,1]
	close(ds)
end

# ╔═╡ 2244daac-8c12-4111-a883-94637d601c3a
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=8,axwidth=8)
	
	axs[1].pcolormesh(x,y,olr')
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ 0f4a925e-f772-43d7-937d-2f868cac2e19
begin
	sds = NCDataset(datadir("Aggregation","OUT_STAT","RCE_TroPrecLS-Aggregation-spinup.nc"))
	t   = sds["time"][:]
	z   = sds["z"][:] ./ 1000
	p   = sds["p"][:]
	pp  = sds["PRES"][:]
	pwv = sds["PW"][:]
	close(sds)
end

# ╔═╡ f11b57d1-ca84-4287-9352-1a965b406fc4
begin
	pplt.close(); f2,a2 = pplt.subplots(nrows=2,aspect=4,axwidth=6)
	
	c = a2[1].pcolormesh(t,z,pp.-p,levels=vcat(-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5),extend="both",cmap="RdBu_r")
	a2[1].colorbar(c,label="Pressure Perturbation / hPa")
	a2[1].format(xlim=(0,150),xlabel="Time / Days",ylabel="Height / km")

	a2[2].plot(t,pwv)
	a2[2].format(ylim=(20,40))
	
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ e6f6f18b-0816-4c72-8319-e75651b0030f
md"
### B. Loading the ERA5 Data again ...
"

# ╔═╡ 3a6d9c73-5336-4319-a33e-486befa25bd7
e5ds = ERA5Monthly(start=Date(1979),stop=Date(2021),path=datadir())

# ╔═╡ 10e704a9-d93f-4ba8-b484-00b42981c13e
evar = PressureVariable("u",hPa=0,throw=false)

# ╔═╡ 58fc29ed-b84e-4f15-80e3-8897befa554d
egeo = ERA5Region(GeoRegion("TRP"))

# ╔═╡ f86de49e-6a25-4b36-91f7-08a2a1ee56a8
lsd = getLandSea(egeo,path=datadir("emask"))

# ╔═╡ 50f205ec-ebb0-42e8-93cc-d56bb268439d
begin
	umat = zeros(length(lsd.lon),length(lsd.lat),37)
	plist = era5Pressures()[6:end]
	ip = 0
	for p in plist
		global ip += 1
		evar_ii = PressureVariable(evar.varID,hPa=p)
		for yr = 1979 : 2021
			ds = read(e5ds,evar_ii,egeo,Date(yr))
			umat[:,:,ip] += dropdims(mean(ds["u"][:],dims=3),dims=3)
			close(ds)
		end
	end
	umat = umat[:,:,1:ip]
	umat = umat ./ 42
	md"Find the climatological mean from 1979 to 2021 ..."
end

# ╔═╡ 6c3181d2-f0e1-4547-bdb0-67d0cf93b507
begin
	uzon = dropdims(mean(umat[:,121 .+ (-40:40),:],dims=2),dims=2)
	md"Let us do the meridional mean ..."
end

# ╔═╡ b4fca708-a325-489d-9cc3-73261c9d2b9b
begin
	uprofile = dropdims(mean(uzon[241:361,:],dims=1),dims=1)
	uprofile[1:4] .= 0
	md"And let us take a zonal-mean profile from 60-90º"
end

# ╔═╡ f655fff3-e85f-4de1-afab-e058aa4139aa
md"
### C. Finding the vertical profile of zonal-wind
"

# ╔═╡ a571be3e-0abe-4d68-a46c-17d5db82a03f
uu = vcat(uprofile[vcat(1:3,5:9,17:23)],0)

# ╔═╡ 36fafe0e-4455-44c7-be68-5bb9745aefe0
ep = log10.(vcat(plist[vcat(1:3,5:9,17:23)],1009.32))

# ╔═╡ 5372afc4-83eb-4674-a3c6-53b4244f3067
spl = Spline1D(ep,uu,k=2)

# ╔═╡ c617470c-236e-43c2-9930-e4ddcebed69f
begin
	u = evaluate(spl,log10.(dropdims(mean(pp[:,561:600],dims=2),dims=2)))
	u[p.<30] .= 0
	lsfdata = lsfinit(length(z))
	lsfdata[:,1] .= z * 1000
	lsfdata[:,2] .= p
	lsfdata[:,5] .= u

	lsfprint(projectdir("exp","lsf","uforcing"),lsfdata,1009.32)
end

# ╔═╡ 87494580-ea7d-4f08-b35f-587d092aaca7
begin
	pplt.close(); f3,a3 = pplt.subplots(aspect=2/5,axwidth=1)
	
	a3[1].plot(u,z)
	a3[1].format(
		ylim=(0,25),ylabel="Height / km",
		xlim=(-15,15),xlabel=L"u / m s$^{-1}$"
	)
	
	f3.savefig(plotsdir("08a-zonalwindSAM.png"),transparent=false,dpi=200)
	load(plotsdir("08a-zonalwindSAM.png"))
end

# ╔═╡ Cell order:
# ╟─e0fee7ac-fc9b-11ec-0741-9f810d7bcd4c
# ╟─8f2172e0-5b81-420a-8301-dbb8ed0c290a
# ╟─bb51d00a-ee75-4bd8-9744-381b87c6441b
# ╟─2f1ee85d-cef1-4fd2-ba06-7c5ec1a3d4a5
# ╟─cf2e61e6-b46f-4b77-960a-e1d8a69aeefe
# ╟─6c4ec7a4-bac6-4033-81c0-26af845cc11f
# ╟─2244daac-8c12-4111-a883-94637d601c3a
# ╟─0f4a925e-f772-43d7-937d-2f868cac2e19
# ╟─f11b57d1-ca84-4287-9352-1a965b406fc4
# ╟─e6f6f18b-0816-4c72-8319-e75651b0030f
# ╟─3a6d9c73-5336-4319-a33e-486befa25bd7
# ╟─10e704a9-d93f-4ba8-b484-00b42981c13e
# ╟─58fc29ed-b84e-4f15-80e3-8897befa554d
# ╟─f86de49e-6a25-4b36-91f7-08a2a1ee56a8
# ╟─50f205ec-ebb0-42e8-93cc-d56bb268439d
# ╟─6c3181d2-f0e1-4547-bdb0-67d0cf93b507
# ╟─b4fca708-a325-489d-9cc3-73261c9d2b9b
# ╟─f655fff3-e85f-4de1-afab-e058aa4139aa
# ╟─a571be3e-0abe-4d68-a46c-17d5db82a03f
# ╟─36fafe0e-4455-44c7-be68-5bb9745aefe0
# ╟─5372afc4-83eb-4674-a3c6-53b4244f3067
# ╟─c617470c-236e-43c2-9930-e4ddcebed69f
# ╟─87494580-ea7d-4f08-b35f-587d092aaca7
