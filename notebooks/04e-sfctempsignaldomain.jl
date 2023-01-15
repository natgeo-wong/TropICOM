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
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ ebbbea99-9f3d-402a-ba91-8f66914d8478
md"
# 04e. Signals in Surface Temperatures

Following from notebook `04d`, we now seek to decompose the domain-mean, ocean-mean and land-mean timeseries signals for both skin and sea surface temperatures over the deep-tropical (<10$\degree$) regions of the Maritime Continent and Indo-Pacific warmpool.
"

# ╔═╡ 43882c17-36c4-4cab-974e-e518cce6f348
function hr2daily(data)

	return dropdims(mean(reshape(data,24,:),dims=1),dims=1)
	
end

# ╔═╡ 0d258e07-3bec-4dcc-8d45-f866858a3f9e
md"
### A. Predefining and Loading ERA5 Data
"

# ╔═╡ 01907cf6-48cf-4405-a2c5-d6e29c081878
e5ds = ERA5Hourly(start=Date(1979,1),stop=Date(2021,12),path=datadir())

# ╔═╡ 5ce73ef5-e423-462e-9836-6bb134602cc6
evar_skt = SingleVariable("skt")

# ╔═╡ f949e7ab-190b-4e7e-a541-b8beb6faf323
evar_sst = SingleVariable("sst")

# ╔═╡ cfa51d00-1d23-41df-bb63-907e21e76054
egeo = ERA5Region(GeoRegion("DTP_IPW"))

# ╔═╡ 3f62e505-1032-4eaf-a20b-0244ff234244
begin
	ds,_  = read(e5ds,evar_sst,egeo,timeseries=true)
	t     = Date(e5ds.start) : Day(1) : (Date(e5ds.stop))
	sst_d = hr2daily(nomissing(ds["sst_domain"][:],NaN))
	sst_o = hr2daily(nomissing(ds["sst_ocean"][:],NaN))
	sst_l = hr2daily(nomissing(ds["sst_land"][:],NaN))
	close(ds)
	
	ds,_  = read(e5ds,evar_skt,egeo,timeseries=true)
	skt_d = hr2daily(nomissing(ds["skt_domain"][:],NaN))
	skt_o = hr2daily(nomissing(ds["skt_ocean"][:],NaN))
	skt_l = hr2daily(nomissing(ds["skt_land"][:],NaN))
	skt_l_raw = nomissing(ds["skt_land"][:],NaN)
	close(ds)
	md"Loading Sea Surface and Skin Temperature data ..."
end

# ╔═╡ 66a2a61d-978f-4bee-9e9c-79ccc582a902
md"
### B. Basic Exploration of Time-Series Data
"

# ╔═╡ c78b99e9-67f8-4a81-9ff3-72066b62fdac
begin
	pplt.close()
	f1,a1 = pplt.subplots([2,2,1],aspect=1,axwidth=2,sharey=0,sharex=0)
	
	a1[2].plot(t,sst_d,lw=1,c="k",legend="b",label="SST")
	a1[2].plot(t,skt_o,lw=0.5,legend="b",label="Ocean")
	a1[2].plot(t,skt_d,lw=0.5,legend="b",label="Domain")
	a1[2].plot(t,skt_l,lw=0.5,legend="b",label="Land",legend_kw=Dict("frame"=>false,"ncol"=>4))
	a1[2].format(xlabel="Date",ylabel="Temperature / K")

	a1[1].scatter(sst_d.-mean(sst_d),skt_o.-mean(skt_o),s=1,alpha=0.01)
	a1[1].scatter(sst_d.-mean(sst_d),skt_d.-mean(skt_d),s=1,alpha=0.01)
	a1[1].scatter(sst_d.-mean(sst_d),skt_l.-mean(skt_l),s=1,alpha=0.01)
	a1[1].format(
		xlim=(-1.5,1.5),xlocator=-1.5:0.5:1.5,ylabel=L"$\Delta$ Skin Temperature / K",
		ylim=(-1.5,1.5),ylocator=-1.5:0.5:1.5,xlabel=L"$\Delta$ SST / K",
	)
	
	f1.savefig(plotsdir("04e-sfctemp-timeseries.png"),transparent=false,dpi=400)
	load(plotsdir("04e-sfctemp-timeseries.png"))
end

# ╔═╡ fc3dbb3b-d105-4c5a-a6ab-57c180d7234c
md"We see here that the ocean skin temperature and the sea surface temperature are very close indeed, both in the characteristics of their variability and in their mean values. Since the `DTP_IPW` region is dominated by ocean grid points (since it comprises of islands and the warmpool), the domain mean is of course more close to the ocean timeseries compared to the land timeseries.

We see that the land-averaged timeseries is both somewhat phase-shifted from the ocean-averaged and the sea surface temperature timeseries, and that the variability is markedly different even on timescales much longer than daily scales."

# ╔═╡ bfce18f3-b099-4b79-98cc-634c612142f1
md"
### C. Fourier Decomposition of the Time-Series
"

# ╔═╡ 70eeff1a-8f23-4965-9923-55902b82fbec
begin
	pdg_sst = periodogram(sst_o.-mean(sst_o),fs=365)
	pdg_sko = periodogram(skt_o.-mean(skt_o),fs=365)
	pdg_skd = periodogram(skt_d.-mean(skt_d),fs=365)
	pdg_skl = periodogram(skt_l.-mean(skt_l),fs=365)
	md"Fourier Decompositions of the TimeSeries Data ..."
end

# ╔═╡ 7855c102-54f8-4d9a-8a5f-6c51c0a8dd8a
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=3,aspect=1.5,axwidth=2)
	
	a2[1].plot(pdg_skd.freq,pdg_skd.power,lw=1,alpha=0.8)
	a2[2].plot(pdg_sko.freq,pdg_sko.power,lw=1,alpha=0.8)
	a2[3].plot(pdg_skl.freq,pdg_skl.power,lw=1,alpha=0.8)

	a2[1].format(urtitle="(a) Domain")
	a2[2].format(urtitle="(b) Ocean")
	a2[3].format(urtitle="(c) Land")

	for ax in a2
		ax.plot(pdg_sst.freq,pdg_sst.power,lw=3,c="k",zorder=0)
		ax.format(
			xscale="log",xlim=(0.2,50),ylim=(0.00002,5),yscale="log",
			ylocator=10. .^(-4:0),yformatter="log",
			ylabel=L"Power (K$^2$ / f)",
			xlabel="f / Cycles per Year",
		)
	end
	
	f2.savefig(plotsdir("04e-sfctemp-fourier.png"),transparent=false,dpi=150)
	load(plotsdir("04e-sfctemp-fourier.png"))
end

# ╔═╡ 19cdd04a-33a0-4c38-a272-d248f1f87913
md"
We see that the strength of the signal peaks at several different frequencies:
* Annually
* Semi-annually (twice a year)
* 12 times a year (roughly every month)
* 20 times a year (roughly once every 18 days)
* 28 times a year (roughly once every 12 days)

We also see that the semi-annual peak is much stronger over land than it is over the ocean, but nonetheless this peak is noticeable across all the variables we are analyzing.

Let us now decompose specifically, the sea surface temperature signal, because it seems that the sea surface temperature variability is overall the main driver of surface temperature variability, even over land, at relatively short timescales.
"

# ╔═╡ e4beb500-f3c9-4e97-a5aa-74cbe97a0085
begin
	pplt.close(); f3,a3 = pplt.subplots(aspect=4,axwidth=5)

	responsetype = Bandpass(10,30; fs=365)
	designmethod = Butterworth(1)
	sst_itr = filtfilt(digitalfilter(responsetype, designmethod), sst_o)

	responsetype = Bandpass(1.5,2.5; fs=365)
	designmethod = Butterworth(1)
	sst_san = filtfilt(digitalfilter(responsetype, designmethod), sst_o)

	responsetype = Bandpass(0.9,1.1; fs=365)
	designmethod = Butterworth(1)
	sst_ann = filtfilt(digitalfilter(responsetype, designmethod), sst_o)

	responsetype = Lowpass(0.5; fs=365)
	designmethod = Butterworth(1)
	sst_trd = filtfilt(digitalfilter(responsetype, designmethod), sst_o)

	msst = mean(sst_o) * ones(length(t))
	a3[1].plot(t,sst_o.-sst_trd.+msst)
	a3[1].plot(t,sst_san.+mean(sst_o),zorder=3,c="red")
	a3[1].plot(t,sst_ann.+mean(sst_o),zorder=4,c="green")
	a3[1].plot(t,sst_itr.+mean(sst_o),zorder=5,c="purple")
	a3[1].plot(t,msst,c="black",zorder=5)
	a3[1].format(xlim=(Date(1979),Date(1995)))

	msst = [1,1]*mean(sst_o)
	p1 = a3[1].panel("r",space=1,width=0.5)
	p1.area([-1,1],msst.-0.50,msst.+0.50,alpha=0.4,zorder=1)
	p1.area([-1,1],msst.-0.30,msst.+0.30,alpha=0.4,zorder=3,c="red")
	p1.area([-1,1],msst.-0.12,msst.+0.12,alpha=0.4,zorder=4,c="green")
	p1.area([-1,1],msst.-0.09,msst.+0.09,alpha=0.4,zorder=6,c="black")
	p1.plot([-1,1],msst,c="black",zorder=5)
	
	f3.savefig(plotsdir("04e-sfctemp-decomposition.png"),transparent=false,dpi=400)
	load(plotsdir("04e-sfctemp-decomposition.png"))
end

# ╔═╡ 789ab32b-8957-4cca-8f16-2a30aa572127
begin
	pplt.close(); f4,a4 = pplt.subplots(nrows=2,ncols=2,aspect=2,axwidth=2.5,sharex=0)
	
	itt = (0 : 0.001 : 86) * pi
	mann = 0.12*sin.(itt.-pi/6)
	msan = - 0.30*sin.(itt.*2)
	mitr = 0.04*sin.(itt.*12).+0.04*sin.(itt.*20).+0.04*sin.(itt.*28)
	fitc = mean(sst_o) .+ mann .+ msan .+ mitr
	sstm = mean(sst_o) * ones(length(t))
	
	a4[1].plot((1:15706)./365,sst_o.-sst_trd)
	a4[3].plot((1:15706)./365,sst_ann)
	a4[2].plot((1:15706),sst_san,label="ERA5 Data",legend="r",legend_kw=Dict("frame"=>false))
	a4[4].plot((1:15706),sst_itr)
	a4[1].plot(itt./(86*pi)*15706/365,mann.+msan.+mitr,c="k",alpha=0.3)
	a4[3].plot(itt./(86*pi)*15706/365,mann,c="k",alpha=0.3)
	a4[2].plot(itt./(86*pi)*15706,msan,c="k",alpha=0.3)
	a4[4].plot(itt./(86*pi)*15706,mitr,c="k",alpha=0.3,label="Modelled Variability",legend="r",legend_kw=Dict("frame"=>false))
	for ax in a4
		ax.format(ylabel=L"$\Delta$ SST")
	end
	a4[1].format(urtitle="(a) Total Variability",xlim=(0,20),ylim=(-1,1))
	a4[3].format(urtitle="(b) Annual Variability",xlim=(0,20),ylim=(-0.25,0.25),xlabel="Years")
	a4[2].format(urtitle="(c) Semi-Annual Variability",xlim=(0,365),xlocator=0:50:365,xminorticks=0:10:365)
	a4[4].format(urtitle="(d) Intraseasonal Variability",xlabel="Days",xlim=(0,365),xlocator=0:50:365)
	
	f4.savefig(plotsdir("04e-sfctemp-model.png"),transparent=false,dpi=400)
	load(plotsdir("04e-sfctemp-model.png"))
end

# ╔═╡ eef9a4b7-7ef6-4255-a026-7af63a5db081
md"It is of course, noticeably harder to model the intraseasonal variability. However, we believe that our modelled intraseasonal variability is able to capture the **_nature_** of the variability of sea surface temperature on an intraseasonal timescale, in that the rate of variability is similar to the actual data, assuming of course that sea surface temperature variability on the annual and semi-annual timescales are not relevant.

We treat the annual and semi-annual cycles as less important when the WTG schemes are being considered because these cycles are of very long timescales compared to the range of timescales over which the WTG approximation acts. Even if the Rayleigh damping rate is 1 or even 0.5 days, the intraseasonal fluctuation of sea surface temperature and the atmospheric profile will still matter because the atmospheric system has inertia."

# ╔═╡ Cell order:
# ╟─ebbbea99-9f3d-402a-ba91-8f66914d8478
# ╟─c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
# ╟─8d72de06-476c-4bcf-99b3-a9b469fac93d
# ╠═43882c17-36c4-4cab-974e-e518cce6f348
# ╟─0d258e07-3bec-4dcc-8d45-f866858a3f9e
# ╟─01907cf6-48cf-4405-a2c5-d6e29c081878
# ╟─5ce73ef5-e423-462e-9836-6bb134602cc6
# ╟─f949e7ab-190b-4e7e-a541-b8beb6faf323
# ╟─cfa51d00-1d23-41df-bb63-907e21e76054
# ╟─3f62e505-1032-4eaf-a20b-0244ff234244
# ╟─66a2a61d-978f-4bee-9e9c-79ccc582a902
# ╟─c78b99e9-67f8-4a81-9ff3-72066b62fdac
# ╟─fc3dbb3b-d105-4c5a-a6ab-57c180d7234c
# ╟─bfce18f3-b099-4b79-98cc-634c612142f1
# ╟─70eeff1a-8f23-4965-9923-55902b82fbec
# ╟─7855c102-54f8-4d9a-8a5f-6c51c0a8dd8a
# ╟─19cdd04a-33a0-4c38-a272-d248f1f87913
# ╟─e4beb500-f3c9-4e97-a5aa-74cbe97a0085
# ╟─789ab32b-8957-4cca-8f16-2a30aa572127
# ╟─eef9a4b7-7ef6-4255-a026-7af63a5db081
