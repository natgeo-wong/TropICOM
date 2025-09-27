### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 8d72de06-476c-4bcf-99b3-a9b469fac93d
begin
	@quickactivate "TroPrecLS"
	using DSP
	using ERA5Reanalysis
	using NCDatasets
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ ebbbea99-9f3d-402a-ba91-8f66914d8478
md"
# 04f. Point Signals in Surface Temperatures

Following from notebook `04e`, we now seek to decompose the daily timeseries signals for sea surface temperatures over the deep-tropical (<10$\degree$) regions of the Maritime Continent and Indo-Pacific warmpool for the individual points (i.e., not the domain mean). We have downloaded 1.0º x 1.0º resolution data within the region, as this corresponds to domains roughly 100 km in size.
"

# ╔═╡ 0d258e07-3bec-4dcc-8d45-f866858a3f9e
md"
### A. Predefining and Loading ERA5 Data
"

# ╔═╡ 01907cf6-48cf-4405-a2c5-d6e29c081878
e5ds = ERA5Daily(start=Date(1979,1),stop=Date(2021,12),path=datadir())

# ╔═╡ 5ce73ef5-e423-462e-9836-6bb134602cc6
evar_skt = SingleVariable("skt")

# ╔═╡ f949e7ab-190b-4e7e-a541-b8beb6faf323
evar_sst = SingleVariable("sst")

# ╔═╡ cfa51d00-1d23-41df-bb63-907e21e76054
egeo = ERA5Region(GeoRegion("DTP_IPW"),gres=1.0)

# ╔═╡ d33bd045-d504-481c-a0e6-c107085e3788
lsd = getLandSea(egeo,path=datadir("emask"))

# ╔═╡ 3f62e505-1032-4eaf-a20b-0244ff234244
begin
	sstdy = zeros(length(lsd.lon),length(lsd.lat),Dates.value(e5ds.stop-e5ds.start)+1)
	for idt in e5ds.start : Month(1) : e5ds.stop
		ibeg = Dates.value(idt-e5ds.start) + 1
		iend = Dates.value(idt+Month(1)-e5ds.start)
		ndy  = daysinmonth(idt)
		ids = read(e5ds,evar_sst,egeo,idt)
		sstdy[:,:,ibeg:iend] = nomissing(ids[evar_sst.varID][:],NaN)
		close(ids)
	end
end

# ╔═╡ 66a2a61d-978f-4bee-9e9c-79ccc582a902
md"
### B. Basic Exploration of Time-Series Data
"

# ╔═╡ fc3dbb3b-d105-4c5a-a6ab-57c180d7234c
md"We see here that the ocean skin temperature and the sea surface temperature are very close indeed, both in the characteristics of their variability and in their mean values. Since the `DTP_IPW` region is dominated by ocean grid points (since it comprises of islands and the warmpool), the domain mean is of course more close to the ocean timeseries compared to the land timeseries.

We see that the land-averaged timeseries is both somewhat phase-shifted from the ocean-averaged and the sea surface temperature timeseries, and that the variability is markedly different even on timescales much longer than daily scales."

# ╔═╡ f2799fba-35a5-48f4-b303-8d89fb44841e
begin
	pplt.close(); f1,a1 = pplt.subplots(aspect=3)
	
	a1[1].plot(e5ds.start:Day(1):e5ds.stop,sstdy[30,11,:])
	a1[1].format(xlim=(Date(1980),Date(2021)))
	
	f1.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ 23b400e2-4e93-42cc-b1dc-7f941d774e95
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=1.5,axwidth=2)

	for ilat = 1 : length(lsd.lat), ilon = 1 : length(lsd.lon)
		pdg = periodogram(sstdy[ilon,ilat,:].-mean(sstdy[ilon,ilat,:]),fs=365)
		a2[1].plot(pdg.freq,pdg.power.*pdg.freq,lw=1,c="k")
	end
	a2[1].format(
		xscale="log",xlim=(0.1,180),yscale="log",ylim=(0.1,500),
		xlocator=[0.5,1,2,12,20,28],ylabel="Power",xlabel="Cycles per Year"
	)
	
	f2.savefig("test2.png",transparent=false,dpi=150)
	load("test2.png")
end

# ╔═╡ bfce18f3-b099-4b79-98cc-634c612142f1
md"
### C. Fourier Decomposition of the Time-Series
"

# ╔═╡ 19cdd04a-33a0-4c38-a272-d248f1f87913
md"
We see that the strength of the signal peaks at several different frequencies:
* Annually
* Semi-annually (twice a year)
* 12 times a year (roughly every month)
* 20 times a year (roughly once every 18 days)
* 28 times a year (roughly once every 12 days)

This is similar to the domain-averaged SST. 
"

# ╔═╡ 98151041-24d8-41ce-9e10-104b1f24eb4f
begin
	pplt.close(); f3,a3 = pplt.subplots(aspect=4,axwidth=5)

	sstii = sstdy[30,11,:]
	responsetype = Bandpass(10,30; fs=365)
	designmethod = Butterworth(1)
	sst_itr = filtfilt(digitalfilter(responsetype, designmethod), sstii)

	responsetype = Bandpass(1.8,2.2; fs=365)
	designmethod = Butterworth(1)
	sst_san = filtfilt(digitalfilter(responsetype, designmethod), sstii)

	responsetype = Bandpass(2.8,3.2; fs=365)
	designmethod = Butterworth(1)
	sst_san = sst_san .+ filtfilt(digitalfilter(responsetype, designmethod), sstii)

	responsetype = Bandpass(0.6,1.3; fs=365)
	designmethod = Butterworth(1)
	sst_ann = filtfilt(digitalfilter(responsetype, designmethod), sstii)

	responsetype = Lowpass(0.5; fs=365)
	designmethod = Butterworth(1)
	sst_trd = filtfilt(digitalfilter(responsetype, designmethod), sstii)

	t = e5ds.start:Day(1):e5ds.stop
	msst = mean(sstii) * ones(length(t))
	# a3[1].plot(t,sstdy[30,11,:].-sst_trd.+msst)
	a3[1].plot(t,sst_san.+mean(sstii),zorder=3,c="red")
	# a3[1].plot(t,sst_ann.+mean(sstii),zorder=4,c="green")
	# a3[1].plot(t,sst_itr.+mean(sstii),zorder=5,c="purple")
	a3[1].plot(t,msst,c="black",zorder=5)
	a3[1].format(xlim=(Date(1979),Date(2022)))

	msst = [1,1]*mean(sstii)
	p1 = a3[1].panel("r",space=1,width=0.5)
	p1.area([-1,1],msst.-0.50,msst.+0.50,alpha=0.4,zorder=1)
	p1.area([-1,1],msst.-0.30,msst.+0.30,alpha=0.4,zorder=3,c="red")
	p1.area([-1,1],msst.-0.12,msst.+0.12,alpha=0.4,zorder=4,c="green")
	a3[1].area([Date(1979),Date(2021,12,31)],msst.-0.30,msst.+0.30,alpha=0.4,zorder=6,c="black")
	p1.area([-1,1],msst.-0.30,msst.+0.30,alpha=0.4,zorder=6,c="black")
	p1.plot([-1,1],msst,c="black",zorder=5)
	
	f3.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ eef9a4b7-7ef6-4255-a026-7af63a5db081
md"But we see her from a sample point that the amplitude of the intraseasonal cycles are roughly of the same order of magnitude as the semiannual and annual cycles."

# ╔═╡ f5784e23-00cb-4d4b-af6f-c08aee2fe2d0
begin
	pplt.close(); f4,a4 = pplt.subplots(nrows=2,ncols=2,aspect=2,axwidth=2.5,sharex=0)
	
	itt = (0 : 0.001 : 86) * pi
	mann = 0.5*sin.(itt.-pi/2)
	msan = - 0.30*sin.(itt.*2) .+ 0.10*sin.(itt.*3 .-pi/3)
	mitr = 0.15*sin.(itt.*12).+0.15*sin.(itt.*28)
	fitc = mean(sstii) .+ mann .+ msan .+ mitr
	sstm = mean(sstii) * ones(length(t))
	
	# a4[1].plot((1:15706)./365,sstii.-sst_trd)
	a4[3].plot((1:15706)./365,sst_ann)
	a4[2].plot((1:15706),sst_san,label="ERA5 Data",legend="r",legend_kw=Dict("frame"=>false))
	a4[4].plot((1:15706),sst_itr)
	a4[1].plot(itt./(86*pi)*15706/365,mann.+msan.+mitr.+0.5,c="k",alpha=0.3)
	a4[3].plot(itt./(86*pi)*15706/365,mann,c="k",alpha=0.3)
	a4[2].plot(itt./(86*pi)*15706,msan,c="k",alpha=0.3)
	a4[4].plot(itt./(86*pi)*15706,mitr,c="k",alpha=0.3,label="Modelled Variability",legend="r",legend_kw=Dict("frame"=>false))
	for ax in a4
		ax.format(ylabel=L"$\Delta$ SST")
	end
	a4[1].format(urtitle="(a) Total Variability",xlim=(0,2),ylim=(-1.5,1.5))
	a4[3].format(urtitle="(b) Annual Variability",xlim=(0,20),ylim=(-0.5,0.5),xlabel="Years")
	a4[2].format(urtitle="(c) Semi-Annual Variability",xlim=(0,365*2),xlocator=0:100:730,xminorticks=0:10:365)
	a4[4].format(urtitle="(d) Intraseasonal Variability",xlabel="Days",xlim=(365*39,365*40),xlocator=0:100:730)
	
	f4.savefig("test.png",transparent=false,dpi=400)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─ebbbea99-9f3d-402a-ba91-8f66914d8478
# ╟─c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
# ╟─8d72de06-476c-4bcf-99b3-a9b469fac93d
# ╟─0d258e07-3bec-4dcc-8d45-f866858a3f9e
# ╟─01907cf6-48cf-4405-a2c5-d6e29c081878
# ╟─5ce73ef5-e423-462e-9836-6bb134602cc6
# ╟─f949e7ab-190b-4e7e-a541-b8beb6faf323
# ╟─cfa51d00-1d23-41df-bb63-907e21e76054
# ╟─d33bd045-d504-481c-a0e6-c107085e3788
# ╟─3f62e505-1032-4eaf-a20b-0244ff234244
# ╟─66a2a61d-978f-4bee-9e9c-79ccc582a902
# ╟─fc3dbb3b-d105-4c5a-a6ab-57c180d7234c
# ╠═f2799fba-35a5-48f4-b303-8d89fb44841e
# ╠═23b400e2-4e93-42cc-b1dc-7f941d774e95
# ╟─bfce18f3-b099-4b79-98cc-634c612142f1
# ╟─19cdd04a-33a0-4c38-a272-d248f1f87913
# ╠═98151041-24d8-41ce-9e10-104b1f24eb4f
# ╟─eef9a4b7-7ef6-4255-a026-7af63a5db081
# ╠═f5784e23-00cb-4d4b-af6f-c08aee2fe2d0
