### A Pluto.jl notebook ###
# v0.19.14

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

# ╔═╡ b82d16fd-c9a9-4e3f-ac00-3f189a5eb249
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 8ca6782b-b5f8-4908-8c40-419f976b6ff2
begin
	@quickactivate "TroPrecLS"
	using Dates
	using DelimitedFiles
	using DSP
	using ERA5Reanalysis
	using NCDatasets
	using PlutoUI
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 1be6a70c-4f5e-11ed-15c0-2db8088f3b9d
md"
# 04d. Exploring the Time Series for Skin Temperature Data

In some work we have done, we find that using WTG Approximation schemes in cloud-resolving models result in the two end-states of self-aggregated RCE: (1) A very wet, precipitating state and (2) a very dry, non-precipitating state. However, analysis of even smaller tropical domains (meaning in regions small enough that the time-average does not vary much) we see that the majority of humidity profiles actually would fall somewhere between the two.

We postulate that this is because of the fact that the domain-mean sea-surface temperature and air-temperature profiles are constantly varying. From previous WTG experiments, we know that even very small perturbations can be used to deliberately force out wet and dry states. Perhaps if the reference sounding profiles (or the sea surface temperature) were varying (intraseasonally of course), then we could force out these states.

But first, we must decompose out the signal from both sea-surface temperature and skin-temperature data, and separate out the intraseasonal/monthly signals (if any) from the diurnal and seasonal variability. We do this via fourier transform using the `DSP.jl` package over the next few notebooks.
"

# ╔═╡ 50e38500-5574-4aa8-9a7f-8aca1bcc7a29
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ 749d7a0e-ba5c-4fde-9323-0714783d3252
md"
### A. The Domain of Interest

Our domain of interest is the deep-tropical region of Southeast Asia and the Indo-Pacifc warmpool, as this merges both our domain of interest (with lots of islands) along with the latitudes ($<10\degree$) where the weak-temperature gradient approximation holds most strongly.
"

# ╔═╡ 1d7a6a61-4c3d-40aa-bdbb-8f84ed44ac1d
geo = GeoRegion("DTP_IPW")

# ╔═╡ 4a161cef-b374-41f1-b44a-10c77abbcd75
lsd = getLandSea(ERA5Region(geo),path=datadir("emask"))

# ╔═╡ e6c078e1-8c52-4398-98f1-5230d4b37c55
begin
	pplt.close()
	asp = (geo.E - geo.W + 2) / (geo.N - geo.S + 2)
	lvl = vcat(10. .^(-2:-1),0.2,0.5,0.9,0.95,0.99)
	fig,axs = pplt.subplots(aspect=asp,axwidth=asp)
	
	c = axs[1].contourf(
		lsd.lon,lsd.lat,lsd.lsm',levels=lvl,
		cmap="delta",extend="both"
	)
	axs[1].plot(x,y,lw=1,c="k")
	axs[1].format(
		xlim=(geo.W-1,geo.E+1),xlabel=L"Longitude / $\degree$",
		ylim=(geo.S-1,geo.N+1),ylabel=L"Latitude / $\degree$"
	)
	axs[1].colorbar(c)
	
	fig.savefig(plotsdir("04d-DTP_IPW.png"),transparent=false,dpi=400)
	load(plotsdir("04d-DTP_IPW.png"))
end

# ╔═╡ 920e88d9-5f07-405a-86f7-121db3c430cb
md"
### B. Climatological Averages for the Domain ...
"

# ╔═╡ 91145e0f-a44e-4fa3-b442-bb27cf982cff
lsd_TRP = getLandSea(ERA5Region(GeoRegion("TRP")),path=datadir("emask"))

# ╔═╡ 1e665edb-a4b8-4d0b-a348-aa0fef8a354e
begin
	skd = NCDataset(datadir("compiled","era5mh-TRPx0.25-skt-compiled.nc"))
	skt = dropdims(mean(nomissing(skd["skt"][:],NaN),dims=3),dims=3)
	close(skd)
	md"Loading Skin Temperature data ..."
end

# ╔═╡ 351c941e-7f1a-4e22-ab41-2d56258de20c
begin
	ssd = NCDataset(datadir("compiled","era5mh-TRPx0.25-sst-compiled.nc"))
	lon = ssd["longitude"][:]
	lat = ssd["latitude"][:]
	sst = dropdims(mean(nomissing(ssd["sst"][:],NaN),dims=3),dims=3)
	close(ssd)
	md"Loading Sea Surface Temperature data ..."
end

# ╔═╡ 68728fde-74da-484c-afd6-fba573021b05
begin
	pplt.close()
	f2,a2 = pplt.subplots(nrows=3,aspect=asp,axwidth=asp)
	
	c2 = a2[2].pcolormesh(
		lsd_TRP.lon,lsd_TRP.lat,dropdims(mean(skt,dims=3),dims=3)',
		levels=300:0.5:305,extend="both"
	)
	a2[1].pcolormesh(
		lsd_TRP.lon,lsd_TRP.lat,dropdims(mean(sst,dims=3),dims=3)',
		levels=300:0.5:305,extend="both"
	)
	c2_2 = a2[3].pcolormesh(
		lsd_TRP.lon,lsd_TRP.lat,dropdims(mean(sst,dims=3),dims=3)'.-dropdims(mean(skt,dims=3),dims=3)',
		levels=0.02:0.02:0.2,extend="both"
	)

	a2[1].format(urtitle="Sea Surface Temperature (SST) / K")
	a2[2].format(urtitle="Skin Temperature (SKT) / K")
	a2[3].format(urtitle="SST - SKT / K")

	for ax in a2
		ax.plot(x,y,lw=1,c="k")
		ax.format(
			xlim=(geo.W,geo.E),xlabel=L"Longitude / $\degree$",
			ylim=(geo.S,geo.N),ylabel=L"Latitude / $\degree$"
		)
	end
	
	f2.colorbar(c2,locator=300:305,length=0.6,row=(1,2))
	f2.colorbar(c2_2,row=3)
	f2.savefig(plotsdir("04d-DTP_IPW-sfcclimateavg.png"),transparent=false,dpi=400)
	load(plotsdir("04d-DTP_IPW-sfcclimateavg.png"))
end

# ╔═╡ 534bf1be-b842-42cc-a872-335fb4bf0e15
md"So we estimate that the sea surface temperature is roughly 0.2 K warmer than the calculated ocean skin temperature over the Maritime Continent / Indo-Pacific Warmpool. Otherwise howevver, we see that there isn't much of a difference between the two variables. Of course, this is still enough to cause significant issue when the WTG approximation is used. But then again, we're not looking for differences between means, but rather the impact of the seasonal variability. So let's do a month-by-month comparison."

# ╔═╡ f3df206a-49d1-4125-b291-8c99a5621f72
md"
### C. But what about the Seasonal?
"

# ╔═╡ a6bfffd0-8b2d-4c76-af25-986e53b20dbe
begin
	for imo = 1 : 12
		pplt.close()
		f3,a3 = pplt.subplots(nrows=3,aspect=asp,axwidth=asp)
		
		c3 = a3[2].contourf(
			lsd_TRP.lon,lsd_TRP.lat,skt[:,:,imo]',
			levels=300:0.5:305,extend="both"
		)
		a3[1].contourf(
			lsd_TRP.lon,lsd_TRP.lat,sst[:,:,imo]',
			levels=300:0.5:305,extend="both"
		)
		c3_2 = a3[3].contourf(
			lsd_TRP.lon,lsd_TRP.lat,sst[:,:,imo]'.-skt[:,:,imo]',
			levels=0.02:0.02:0.2,extend="both"
		)
	
		a3[1].format(urtitle="Sea Surface Temperature (SST) / K")
		a3[2].format(urtitle="Skin Temperature (SKT) / K")
		a3[3].format(urtitle="SST - SKT / K")
	
		for ax in a3
			ax.plot(x,y,lw=1,c="k")
			ax.format(
				xlim=(geo.W,geo.E),xlabel=L"Longitude / $\degree$",
				ylim=(geo.S,geo.N),ylabel=L"Latitude / $\degree$"
			)
		end
		
		f3.colorbar(c3,locator=300:305,length=0.6,row=(1,2))
		f3.colorbar(c3_2,row=3)
		f3.savefig(
			plotsdir("04d-DTP_IPW-sfcclimateavg-$(monthabbr(imo)).png"),
			transparent=false,dpi=400
		)
	end
end

# ╔═╡ e80aa57d-bf5a-41b7-b966-42474e8a4866
md"Month: $(@bind mmo PlutoUI.Slider(1:12,default=1, show_value=true))"

# ╔═╡ 91da6dc6-4007-402c-ac8e-b8348c14d963
load(plotsdir("04d-DTP_IPW-sfcclimateavg-$(monthabbr(mmo)).png"))

# ╔═╡ a77aa352-df99-47ee-b1de-ddef9978dd5f
md"
We see that the differences in sea surface and skin temperature are smallest around island regions which have the greatest precipitation. This is interesting because for the warmpool regions where there is always high rainfall the sea surface temperature is around 0.2 K higher than skin temperature.
"

# ╔═╡ 9f85bd6d-794c-4347-b20a-2374a6452551
md"
### D. Comparison against Domain Mean SST
"

# ╔═╡ 9d3fb2fe-16a3-4107-b4f0-b4a37491b798
ggrd = RegionGrid(geo,lsd_TRP.lon,lsd_TRP.lat)

# ╔═╡ 8cb9d654-c6e4-44a3-816b-9155ac5666dc
begin
	for imo = 1 : 12
		pplt.close()
		f4,a4 = pplt.subplots(nrows=2,aspect=asp,axwidth=asp)
		
		sstii = extractGrid(sst[:,:,imo],ggrd)
		sstm  = mean(sstii[.!isnan.(sstii)])
		c4_1 = a4[1].contourf(
			lsd.lon,lsd.lat,sstii',
			levels=300:0.5:305,extend="both"
		)
		c4_2 = a4[2].contourf(
			lsd.lon,lsd.lat,sstii'.-sstm,
			levels=vcat(-5:-1,-0.5,0.5,1:5),extend="both",cmap="RdBu_r"
		)
	
		a4[1].format(urtitle="SST / K")
		a4[2].format(urtitle=L"SST - $\mu$(SST) / K")
	
		for ax in a4
			ax.plot(x,y,lw=1,c="k")
			ax.format(
				xlim=(geo.W,geo.E),xlabel=L"Longitude / $\degree$",
				ylim=(geo.S,geo.N),ylabel=L"Latitude / $\degree$"
			)
		end
		
		a4[1].colorbar(c4_1,locator=300:305)
		a4[2].colorbar(c4_2,locator=[-5,-3,-1,1,3,5])
		f4.savefig(
			plotsdir("04d-DTP_IPW-sstdelta-$(monthabbr(imo)).png"),
			transparent=false,dpi=400
		)
	end
end

# ╔═╡ 1225457c-9ac2-4d78-aa7d-0fde23427b6e
load(plotsdir("04d-DTP_IPW-sstdelta-$(monthabbr(mmo)).png"))

# ╔═╡ b1f4cde1-345c-4315-966b-e542dc708c0a
md"
### E. Some additional Areas of Interest based on D

Based on an analysis of D, we come up with certain areas of interest for a more in depth analysis.
"

# ╔═╡ 37511844-2a99-4180-884a-8ec1e5eb0b30
begin
	if isGeoRegion("DTP_IPW_1",throw=false)
		removeGeoRegion("DTP_IPW_1")
	end
	geo1 = RectRegion(
		"DTP_IPW_1","DTP_IPW","IndoPacific Warmpool SubRegion 1", [7.5,5,110,107.5]
	)
end

# ╔═╡ ae856c1e-84c8-479b-9964-12788ff2b8f8
begin
	if isGeoRegion("DTP_IPW_2",throw=false)
		removeGeoRegion("DTP_IPW_2")
	end
	geo2 = RectRegion(
		"DTP_IPW_2","DTP_IPW","IndoPacific Warmpool SubRegion 1", [4.5,2,122,119.5]
	)
end

# ╔═╡ 609b663c-c736-46e6-892a-e3abadf4c053
begin
	if isGeoRegion("DTP_IPW_3",throw=false)
		removeGeoRegion("DTP_IPW_3")
	end
	geo3 = RectRegion(
		"DTP_IPW_3","DTP_IPW","IndoPacific Warmpool SubRegion 3", [1.25,-1.25,152.5,150]
	)
end

# ╔═╡ 40c4534d-8a0f-4ee3-8f85-15779b1b92d0
begin
	if isGeoRegion("DTP_IPW_4",throw=false)
		removeGeoRegion("DTP_IPW_4")
	end
	geo4 = RectRegion(
		"DTP_IPW_4","DTP_IPW","IndoPacific Warmpool SubRegion 4", [-6,-8.5,137.5,135]
	)
end

# ╔═╡ 51181632-1b49-48a4-95b1-1f207ea7ca92
begin
	sln1,slt1 = coordGeoRegion(GeoRegion("DTP_IPW_1"))
	sln2,slt2 = coordGeoRegion(GeoRegion("DTP_IPW_2"))
	sln3,slt3 = coordGeoRegion(GeoRegion("DTP_IPW_3"))
	sln4,slt4 = coordGeoRegion(GeoRegion("DTP_IPW_4"))
	md"Loading coordinates for subregions ..."
end

# ╔═╡ b17fca74-f18a-45b3-8660-fbb2b113d34a
begin
	pplt.close()
	f5,a5 = pplt.subplots(aspect=asp,axwidth=asp)
	
	c5 = a5[1].contourf(
		lsd.lon,lsd.lat,lsd.lsm',levels=lvl,
		cmap="delta",extend="both"
	)
	a5[1].plot(x,y,lw=1,c="k")
	a5[1].plot(sln1,slt1,lw=1,c="r")
	a5[1].plot(sln2,slt2,lw=1,c="g")
	a5[1].plot(sln3,slt3,lw=1,c="b")
	a5[1].plot(sln4,slt4,lw=1,c="y")
	a5[1].format(
		xlim=(geo.W-1,geo.E+1),xlabel=L"Longitude / $\degree$",
		ylim=(geo.S-1,geo.N+1),ylabel=L"Latitude / $\degree$"
	)
	a5[1].colorbar(c5)
	
	f5.savefig(plotsdir("04d-IPWsubregions.png"),transparent=false,dpi=400)
	load(plotsdir("04d-IPWsubregions.png"))
end

# ╔═╡ Cell order:
# ╟─1be6a70c-4f5e-11ed-15c0-2db8088f3b9d
# ╟─b82d16fd-c9a9-4e3f-ac00-3f189a5eb249
# ╟─8ca6782b-b5f8-4908-8c40-419f976b6ff2
# ╟─50e38500-5574-4aa8-9a7f-8aca1bcc7a29
# ╟─749d7a0e-ba5c-4fde-9323-0714783d3252
# ╟─1d7a6a61-4c3d-40aa-bdbb-8f84ed44ac1d
# ╟─4a161cef-b374-41f1-b44a-10c77abbcd75
# ╟─e6c078e1-8c52-4398-98f1-5230d4b37c55
# ╟─920e88d9-5f07-405a-86f7-121db3c430cb
# ╟─91145e0f-a44e-4fa3-b442-bb27cf982cff
# ╟─1e665edb-a4b8-4d0b-a348-aa0fef8a354e
# ╟─351c941e-7f1a-4e22-ab41-2d56258de20c
# ╟─68728fde-74da-484c-afd6-fba573021b05
# ╟─534bf1be-b842-42cc-a872-335fb4bf0e15
# ╟─f3df206a-49d1-4125-b291-8c99a5621f72
# ╟─a6bfffd0-8b2d-4c76-af25-986e53b20dbe
# ╟─e80aa57d-bf5a-41b7-b966-42474e8a4866
# ╟─91da6dc6-4007-402c-ac8e-b8348c14d963
# ╟─a77aa352-df99-47ee-b1de-ddef9978dd5f
# ╟─9f85bd6d-794c-4347-b20a-2374a6452551
# ╟─9d3fb2fe-16a3-4107-b4f0-b4a37491b798
# ╟─8cb9d654-c6e4-44a3-816b-9155ac5666dc
# ╟─1225457c-9ac2-4d78-aa7d-0fde23427b6e
# ╟─b1f4cde1-345c-4315-966b-e542dc708c0a
# ╟─37511844-2a99-4180-884a-8ec1e5eb0b30
# ╟─ae856c1e-84c8-479b-9964-12788ff2b8f8
# ╟─609b663c-c736-46e6-892a-e3abadf4c053
# ╠═40c4534d-8a0f-4ee3-8f85-15779b1b92d0
# ╟─51181632-1b49-48a4-95b1-1f207ea7ca92
# ╟─b17fca74-f18a-45b3-8660-fbb2b113d34a
