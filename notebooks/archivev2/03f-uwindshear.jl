### A Pluto.jl notebook ###
# v0.19.22

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
	using Dates
	using Dierckx
	using DelimitedFiles
	using ERA5Reanalysis
	using NCDatasets
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	include(srcdir("sam.jl"))
	include(srcdir("samsnd.jl"))
	include(srcdir("samlsf.jl"))

md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 90fffbc8-524d-11eb-232a-1bada28d5505
md"
# 3f. Zonal Wind Shear

The first variable that we want to explore among the ERA5 reanalysis variables is Skin Temperature.  More details about the retrieval of reanalysis data can be found in notebook `03-reanalysis.jl`.
"

# ╔═╡ d82366b0-53b1-11eb-26c1-ff1bb6ccb027
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ 3565af3c-5311-11eb-34c4-2d228b05b17c
md"
ERA5 reanalysis is defined according to UTC.  Therefore, for every grid point, we must also correct to the local time before fitting the data to the model.
"

# ╔═╡ 409521c5-a78e-424d-81ae-21e8df7a581c
md"
### B. Defining Datasets and GeoRegions
"

# ╔═╡ aaf4de88-66dd-41c2-8b89-5ed9af777f05
e5ds = ERA5Monthly(start=Date(1979),stop=Date(2021),path=datadir())

# ╔═╡ 12c9c7c3-eeec-40f2-beeb-cd43b3dafeef
pressure = 250

# ╔═╡ 3a1cd5ea-1828-48f8-b50d-2c3b933926ec
evar = PressureVariable("u",hPa=pressure,throw=false)

# ╔═╡ 37d2e058-22dc-4bd3-b992-b1037d2b8ada
egeo = ERA5Region(GeoRegion("TRP"))

# ╔═╡ 1fadf4ca-5755-11eb-1ece-a99313019785
lsd = getLandSea(egeo,path=datadir("emask"))

# ╔═╡ aa05317e-530b-11eb-2ec1-93aff65659dd
md"
### C. Retrieving ERA5 Skin Temperature Data
"

# ╔═╡ 5c0e5bae-554e-11eb-3f83-a364ae0a2485
md"
Test
"

# ╔═╡ 1eff88e0-0d85-4d4e-9dfa-b140aacaf2fd
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

# ╔═╡ 4d73a001-5df5-4231-8f8b-f7614b92522e
begin
	asp1 = (maximum(lsd.lon)-minimum(lsd.lon))/(maximum(lsd.lat)-minimum(lsd.lat))
	pplt.close(); f1,a1 = pplt.subplots(aspect=18,axwidth=asp1*1.2,nrows=4)

	i100 = argmin(abs.(plist.-100))
	i250 = argmin(abs.(plist.-250))
	i500 = argmin(abs.(plist.-500))
	i850 = argmin(abs.(plist.-850))
	
	c = a1[1].pcolormesh(
		lsd.lon,lsd.lat,umat[:,:,i100]',levels=vcat(-12:2:-2,2:2:12),
		cmap="RdBu_r",extend="both"
	)

	a1[2].pcolormesh(
		lsd.lon,lsd.lat,umat[:,:,i250]',levels=vcat(-12:2:-2,2:2:12),
		cmap="RdBu_r",extend="both"
	)

	a1[3].pcolormesh(
		lsd.lon,lsd.lat,umat[:,:,i500]',levels=vcat(-12:2:-2,2:2:12),
		cmap="RdBu_r",extend="both"
	)

	a1[4].pcolormesh(
		lsd.lon,lsd.lat,umat[:,:,i850]',levels=vcat(-12:2:-2,2:2:12),
		cmap="RdBu_r",extend="both"
	)
	
	a1[1].format(urtitle="(a) 100 hPa")
	a1[2].format(urtitle="(b) 250 hPa")
	a1[3].format(urtitle="(c) 500 hPa")
	a1[4].format(urtitle="(d) 850 hPa")
	
	for ax in a1
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(
			xlim=(0,360),xlabel=L"Longitude / $\degree$",xlocator=-0:60:360,
			ylim=(-10,10),ylabel=L"Latitude / $\degree$"
		)
	end
	
	f1.colorbar(c,loc="r",length=0.8,label=L"u / m s$^{-1}$")
	f1.savefig(plotsdir("03f-uwind_$(pressure)hPa-$(egeo.geoID).png"),transparent=false,dpi=400)
	PNGFiles.load(plotsdir("03f-uwind_$(pressure)hPa-$(egeo.geoID).png"))
end

# ╔═╡ e2295eb9-3038-4aba-a5bb-e2cfcc8bcd04
begin
	uzon = dropdims(mean(umat[:,121 .+ (-40:40),:],dims=2),dims=2)
	md"Let us do the meridional mean ..."
end

# ╔═╡ ccf0f9b5-f9e2-4dd1-8f9e-ede9a67246af
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=3,axwidth=6)

	c2 = a2[1].contourf(
		lsd.lon,era5Pressures()[6:end],uzon',
		extend="both",levels=vcat(-12:2:-2,2:2:12)
	)
	a2[1].format(
		xlim=(0,360),xlabel=L"Longitude / $\degree$",xlocator=0:60:360,
		ylim=(1000,10),ylabel="Pressure / hPa",yscale="log"
	)
	
	f2.colorbar(c2,label=L"u / m s$^{-1}$")
	f2.savefig(plotsdir("03f-uwindshear_zonalplot.png"),transparent=false,dpi=400)
	load(plotsdir("03f-uwindshear_zonalplot.png"))
end

# ╔═╡ 65447745-2535-41d6-bfbc-90c0a1a1c73c
begin
	uprofile = dropdims(mean(uzon[241:361,:],dims=1),dims=1)
	uprofile[1:4] .= 0
	md"And let us take a zonal-mean profile from 60-90º"
end

# ╔═╡ 4b2535c8-ac84-423c-87c3-25b91a4ba790
begin
	pplt.close(); f3,a3 = pplt.subplots(aspect=1/3,axwidth=1)

	a3[1].plot(uprofile,plist)

	uu = vcat(uprofile[vcat(1:3,5:9,17:23)],0)
	pp = log10.(vcat(plist[vcat(1:3,5:9,17:23)],1009.32))
	a3[1].scatter(uu,10. .^pp)
	spl = Spline1D(pp,uu,k=2)
	uu = evaluate(spl,vcat(1:0.01:3,log10(1009.32)))
	a3[1].plot(uu,10. .^vcat(1:0.01:3,log10(1009.32)))
	a3[1].format(ylim=(1000,10),yscale="log",xlim=(-15,15))

	f3.savefig(plotsdir("03f-uwindshear_profile.png"),transparent=false,dpi=200)
	load(plotsdir("03f-uwindshear_profile.png"))
end

# ╔═╡ 68cfc46c-5755-11eb-1702-373942539652
md"
### D. Add to SAM Large-Scale Forcing
"

# ╔═╡ b73feeaf-a98f-4815-b6d0-8604172f4353
begin
	z,p,_,_,_,_ = readsnd("islandsize3D")
	u = evaluate(spl,log10.(p)); u[p.<30] .= 0
	lsfdata = lsfinit(length(z))
	lsfdata[:,1] .= z
	lsfdata[:,2] .= p
	lsfdata[:,5] .= u

	lsfprint(projectdir("exp","lsf","uforcing"),lsfdata,1009.32)
end

# ╔═╡ Cell order:
# ╟─90fffbc8-524d-11eb-232a-1bada28d5505
# ╟─6f00b8fc-530c-11eb-2242-99d8544f6e14
# ╟─8f30c56c-530c-11eb-2782-33f3c4ed9e89
# ╟─d82366b0-53b1-11eb-26c1-ff1bb6ccb027
# ╟─3565af3c-5311-11eb-34c4-2d228b05b17c
# ╟─409521c5-a78e-424d-81ae-21e8df7a581c
# ╟─aaf4de88-66dd-41c2-8b89-5ed9af777f05
# ╠═12c9c7c3-eeec-40f2-beeb-cd43b3dafeef
# ╟─3a1cd5ea-1828-48f8-b50d-2c3b933926ec
# ╟─37d2e058-22dc-4bd3-b992-b1037d2b8ada
# ╟─1fadf4ca-5755-11eb-1ece-a99313019785
# ╟─aa05317e-530b-11eb-2ec1-93aff65659dd
# ╟─5c0e5bae-554e-11eb-3f83-a364ae0a2485
# ╟─1eff88e0-0d85-4d4e-9dfa-b140aacaf2fd
# ╟─4d73a001-5df5-4231-8f8b-f7614b92522e
# ╟─e2295eb9-3038-4aba-a5bb-e2cfcc8bcd04
# ╟─ccf0f9b5-f9e2-4dd1-8f9e-ede9a67246af
# ╟─65447745-2535-41d6-bfbc-90c0a1a1c73c
# ╟─4b2535c8-ac84-423c-87c3-25b91a4ba790
# ╟─68cfc46c-5755-11eb-1702-373942539652
# ╠═b73feeaf-a98f-4815-b6d0-8604172f4353
