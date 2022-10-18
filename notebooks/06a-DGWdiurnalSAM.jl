### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ ba6d524f-b050-4a2c-9772-536c740e3b27
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 00921025-5d26-4356-ae23-7208719d5634
begin
	@quickactivate "TroPrecLS"
	using DSP
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))

	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 87ad146c-175d-11ed-274a-cb219a4294d1
md"
# 06a. Some Diurnal Stats for SAM
"

# ╔═╡ 55606e82-91e1-4dbe-b1a2-a4b0454d7c79
md"
### A. Loading and Plotting some Basic Variable Statistics
"

# ╔═╡ fd21e1a7-99ca-45a0-b596-b19a05801b18
begin
	slbvec = [
		"Slab00d02","Slab00d05","Slab00d10","Slab00d20",
		"Slab00d50","Slab01d00","Slab02d00","Slab05d00",
		"Slab10d00","Slab20d00","Slab50d00",
	]
	dmpvec = [
		"damping01d0","damping01d4","damping02d0","damping03d2",
		"damping05d0","damping07d1","damping10d0",
	]
	slbplt = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50]
	dmpplt = [1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10]
	nslb = length(slbvec)
	ndmp = length(dmpvec)
	md"Defining mixed-layer slab depth and WTG damping configurations ..."
end

# ╔═╡ 4ac1c68c-01cd-4d55-96fb-ea7dd85bda39
diurnal(data) = dropdims(mean(reshape(data,48,:),dims=2),dims=2)

# ╔═╡ 1c03b953-f2ee-4552-a884-136cf418b5e2
begin
	prcp = zeros(nslb,ndmp)
	pdrn = zeros(nslb,ndmp)
	sst = zeros(nslb,ndmp)
	sdrn = zeros(nslb,ndmp)
	olr = zeros(nslb,ndmp)
	odrn = zeros(nslb,ndmp)
	omax = zeros(nslb,ndmp)
	omin = zeros(nslb,ndmp)
	for idmp = 1 : ndmp, islb = 1 : nslb
	
		prcp[islb,idmp] = mean(diurnal(retrievevar("PREC",slbvec[islb],dmpvec[idmp])))
		sst[islb,idmp] = mean(diurnal(retrievevar("SST",slbvec[islb],dmpvec[idmp])))
		olr[islb,idmp] = mean(diurnal(retrievevar("LWNT",slbvec[islb],dmpvec[idmp])))
		pdrn[islb,idmp] = (maximum(diurnal(retrievevar("PREC",slbvec[islb],dmpvec[idmp]))) -minimum(diurnal(retrievevar("PREC",slbvec[islb],dmpvec[idmp])))) / 2
		sdrn[islb,idmp] = (maximum(diurnal(retrievevar("SST",slbvec[islb],dmpvec[idmp]))) -minimum(diurnal(retrievevar("SST",slbvec[islb],dmpvec[idmp])))) / 2
		odrn[islb,idmp] = (maximum(diurnal(retrievevar("LWNT",slbvec[islb],dmpvec[idmp]))) -minimum(diurnal(retrievevar("LWNT",slbvec[islb],dmpvec[idmp])))) / 2
		omax[islb,idmp] = maximum(diurnal(retrievevar("LWNT",slbvec[islb],dmpvec[idmp])))
		omin[islb,idmp] = minimum(diurnal(retrievevar("LWNT",slbvec[islb],dmpvec[idmp])))
		
	end
end

# ╔═╡ cd6b4d35-f70c-49d4-9f7f-934ae4859cac
begin
	pplt.close(); f1,a1 = pplt.subplots(ncols=2,nrows=2,axwidth=1.5)
	
	c1_1 = a1[1].pcolormesh(
		slbplt,dmpplt,prcp',
		cmap="blues",extend="both",levels=9:0.5:19
	)
	a1[1].colorbar(c1_1,locator=9:19,label=L"mm day$^{-1}$")
	a1[1].format(ltitle="(a) Mean Rain Rate")
	
	c1_2 = a1[2].pcolormesh(
		slbplt,dmpplt,sst',extend="both",levels=295.5:0.5:301
	)
	a1[2].colorbar(c1_2,locator=295:305,label="K")
	a1[2].format(ltitle="(c) Mean SST")

	c1_5 = a1[3].pcolormesh(
		slbplt,dmpplt,pdrn',
		cmap="blues",extend="both",levels=5:2.5:50
	)
	a1[3].colorbar(c1_5,locator=0:5:50,label=L"mm day$^{-1}$")
	a1[3].format(ltitle="(b) Diurnal Amplitude")
	
	c1_6 = a1[4].pcolormesh(
		slbplt,dmpplt,sdrn',
		extend="both",levels=10. .^(-1.5:0.25:1)
	)
	a1[4].colorbar(c1_6,locator=10. .^(-1:1),label="K")
	a1[4].format(ltitle="(d) Diurnal Amplitude")

	for ax in a1
		ax.format(
			xscale="log",xlim=(0.02,50),xlabel="Mixed-Layer Depth / m",
			yscale="log",ylim=(1,10),ylabel=L"$a_m$ / day$^{-1}$",
			ylocator=[1,2,5,10]
		)
	end
	
	f1.savefig(plotsdir("06a-diurnalamp_slabmeans.png"),transparent=false,dpi=300)
	load(plotsdir("06a-diurnalamp_slabmeans.png"))
end

# ╔═╡ 97ee3f27-3bfc-42c7-bc41-e9e0c7acbb82
begin
	pplt.close(); folr,aolr = pplt.subplots(ncols=2,nrows=2,axwidth=1.5)

	c1_3 = aolr[1].pcolormesh(
		slbplt,dmpplt,olr',cmap="greys",extend="both",levels=180:2.5:220
	)
	aolr[1].format(ltitle="(e) Mean OLR")
	aolr[1].colorbar(c1_3,label=L"W m$^{-2}$",locator=180:10:220)
	
	c1_4 = aolr[2].pcolormesh(
		slbplt,dmpplt,omax',cmap="greys",extend="both",levels=200:2:250
	)
	aolr[2].format(ltitle="(g) Maximum OLR")
	aolr[2].colorbar(c1_4,label=L"W m$^{-2}$",locator=200:10:250,minorticks=200:2:250)

	c1_7 = aolr[3].pcolormesh(
		slbplt,dmpplt,odrn',cmap="greys_r",extend="both",levels=10:5:80
	)
	aolr[3].format(ltitle="(f) Diurnal Amplitude")
	aolr[3].colorbar(c1_7,label=L"W m$^{-2}$",locator=15:15:75)
	
	c1_8 = aolr[4].pcolormesh(
		slbplt,dmpplt,omin',cmap="greys",extend="both",levels=100:5:200
	)
	aolr[4].format(ltitle="(h) Minimum OLR")
	aolr[4].colorbar(c1_8,label=L"W m$^{-2}$",locator=100:20:200)

	for ax in aolr
		ax.format(
			xscale="log",xlim=(0.02,50),xlabel="Mixed-Layer Depth / m",
			yscale="log",ylim=(1,10),ylabel=L"$a_m$ / day$^{-1}$",
			ylocator=[1,2,5,10]
		)
	end
	
	folr.savefig(plotsdir("06a-diurnalamp_olr.png"),transparent=false,dpi=300)
	load(plotsdir("06a-diurnalamp_olr.png"))
end

# ╔═╡ 0ca15667-c5f3-4e8c-8c04-2242457ccc73
md"
### B. Timeseries evolution of Model State
"

# ╔═╡ c16de05a-cc0f-4bcc-a1a1-98846b5b6b21
expname = "Slab00d02"

# ╔═╡ 0e8581c6-98b7-4e8d-aea4-0c6089680864
config = "damping01d0"

# ╔═╡ ab4fc2eb-7d24-4a22-8956-1fad02327dbb
begin
	_,p,ts = retrievedims(expname,config)
	prcpts = retrievevar("PREC",expname,config)
	sstts  = retrievevar("SST",expname,config)
	wwtgts = retrievevar("WWTG",expname,config)
	tbias  = retrievevar("TBIAS",expname,config)
	tabsts = retrievevar("TABS",expname,config)
	cldts  = retrievevar("CLD",expname,config) * 100
	qvts   = retrievevar("QV",expname,config) / 1000
	rhts   = calcrh(qvts,tabsts,p*100)
	md"Loading timeseries data ..."
end

# ╔═╡ f6c650d0-ebec-4f8f-8f54-aeff61bb261e
begin
	pplt.close(); fts,ats = pplt.subplots(ncols=2,nrows=3,aspect=1.5,sharey=0)

	expstr = replace(expname,"Slab"=>"")
	expstr = replace(expstr,"d"=>".")
	expstr = parse(Float64,expstr)
	expstr = @sprintf("%4.2f",expstr)

	constr = replace(config,"damping"=>"")
	constr = replace(constr,"d"=>".")
	constr = parse(Float64,constr)
	constr = @sprintf("%3.1f",constr)

	tlvl = vcat(-1*10. .^(1:-0.25:-0.5),10. .^(-0.5:0.25:1)) / 10
	clvl = vcat(5,10:10:90,95)
	
	ats[1].pcolormesh(ts,p,wwtgts,levels=tlvl)
	cts_1 = ats[2].pcolormesh(ts,p,tbias,levels=tlvl,extend="both")
	ats[3].pcolormesh(ts,p,cldts,levels=clvl,extend="both",cmap="Blues_r")
	cts_2 = ats[4].pcolormesh(ts,p,rhts,levels=clvl,extend="both",cmap="Blues")
	ats[5].plot(ts,prcpts)
	ats[6].plot(ts,sstts)

	for ii in 1:4
		ats[ii].format(
			xlim=(280,300),ylim=(1000,25),yscale="log",
			ylabel="Pressure / hPa",
		)
	end
	ats[1].format(ultitle=L"(a) w$_{wtg}$ / 10$^{-1}$ m s$^{-1}$")
	ats[2].format(ultitle=L"(b) T - T$_{obs}$ / K")
	ats[3].format(ultitle="(c) Cloud Fraction / %")
	ats[4].format(ultitle="(d) Relative Humidity / %")
	ats[5].format(ylim=(0,150),ylabel=L"(e) Rain Rate / mm day$^{-1}$")
	ats[6].format(ylabel="(f) SST / K")

	for ax in ats
		ax.format(
			xlabel="Time / Day",
			suptitle="MLD = $(expstr) m | " *
			L"$a_m$" * " = $(constr) " * L"day$^{-1}$"
		)
	end

	ats[2].colorbar(cts_1,locator=[-1,-sqrt(0.1),-0.1,0,0.1,sqrt(0.1),1])
	ats[4].colorbar(cts_2)

	mkpath(plotsdir("06a-timeseries"))
	
	fts.savefig(
		plotsdir("06a-timeseries","06a-timeseries-$(expname)-$(config).png"),
		transparent=false,dpi=300
	)
	load(plotsdir("06a-timeseries","06a-timeseries-$(expname)-$(config).png"))
end

# ╔═╡ ee709bcd-6872-44c7-9c8d-bf44b9c39c30
md"
### C. Signal in Precipitation Data Timeseries
"

# ╔═╡ d3165f91-391b-4d55-9405-d31ac6f35248
begin
	signalpower = zeros(7201,nslb,ndmp)
	signalfreq  = zeros(7201,nslb,ndmp)
	for idmp = 1 : ndmp, islb = 1 : nslb
		prcp = retrievevar("PREC",slbvec[islb],dmpvec[idmp])
		pdg = periodogram(prcp,fs=48)
		signalpower[:,islb,idmp] .= pdg.power
		signalfreq[:,islb,idmp]  .= pdg.freq
	end
end

# ╔═╡ f7857f2e-1a15-45e4-b460-ff2bbfc21ce9
frqbin = 10. .^(-1.01:0.02:1.01)

# ╔═╡ a3569603-9c03-426e-b3a2-f44f90bc65fe
frqplt = 10. .^(-1:0.02:1)

# ╔═╡ 514886a7-9512-4f32-b26d-63c27f010bc8
begin
	powerbin = zeros(length(frqplt),nslb,ndmp)
	for idmp = 1 : ndmp, islb = 1 : nslb, ibin in 1 : length(frqplt)
		ind = (signalfreq[:,islb,idmp] .>= frqbin[ibin]) .& 
				(signalfreq[:,islb,idmp] .<= frqbin[ibin+1])
		powerbin[ibin,islb,idmp] = mean(signalpower[ind,islb,idmp])
	end
end

# ╔═╡ 1708f378-89b2-4d93-bc22-776bcf93ae98
begin
	pplt.close(); f2,a2 = pplt.subplots(
		[[1,2],[3,4],[5,6],[7,0]],aspect=2.5,axwidth=2.5
	)

	c = a2[1].pcolormesh(
		1 ./frqplt,slbplt,log10.(powerbin[:,:,1]'),
		cmap="RdBu_r",levels=1.5:0.1:4,extend="both"
	)
	for idmp = 1 : ndmp
		dmpii = replace(dmpvec[idmp],"damping"=>"")
		dmpii = replace(dmpii,"d"=>".")
		dmpii = parse(Float64,dmpii)
		a2[idmp].pcolormesh(
			1 ./frqplt,slbplt,log10.(powerbin[:,:,idmp]'),
			cmap="RdBu_r",levels=1.5:0.1:4,extend="both"
		)
		a2[idmp].format(ultitle=L"$a_m$" * " = $(@sprintf("%4.1f",dmpii)) " * L"day$^{-1}$")
	end

	for ax in a2
		ax.format(
			yscale="log",xscale="log",xlim=(1e-1,10),ylim=(0.02,50),
			suptitle="Periodogram of Precipitation Rate",
			ylabel="Mixed Layer Depth / m",xlabel="Period / days"
		)
	end

	f2.colorbar(c,locator=1:0.5:4,length=0.4,label=L"Power Density / 10$^4$ dB")
	f2.savefig(plotsdir("06a-precipspectrum.png"),transparent=false,dpi=300)
	load(plotsdir("06a-precipspectrum.png"))
end

# ╔═╡ Cell order:
# ╟─87ad146c-175d-11ed-274a-cb219a4294d1
# ╟─ba6d524f-b050-4a2c-9772-536c740e3b27
# ╟─00921025-5d26-4356-ae23-7208719d5634
# ╟─55606e82-91e1-4dbe-b1a2-a4b0454d7c79
# ╠═fd21e1a7-99ca-45a0-b596-b19a05801b18
# ╠═4ac1c68c-01cd-4d55-96fb-ea7dd85bda39
# ╟─1c03b953-f2ee-4552-a884-136cf418b5e2
# ╟─cd6b4d35-f70c-49d4-9f7f-934ae4859cac
# ╟─97ee3f27-3bfc-42c7-bc41-e9e0c7acbb82
# ╟─0ca15667-c5f3-4e8c-8c04-2242457ccc73
# ╠═c16de05a-cc0f-4bcc-a1a1-98846b5b6b21
# ╠═0e8581c6-98b7-4e8d-aea4-0c6089680864
# ╟─ab4fc2eb-7d24-4a22-8956-1fad02327dbb
# ╟─f6c650d0-ebec-4f8f-8f54-aeff61bb261e
# ╟─ee709bcd-6872-44c7-9c8d-bf44b9c39c30
# ╟─d3165f91-391b-4d55-9405-d31ac6f35248
# ╟─f7857f2e-1a15-45e4-b460-ff2bbfc21ce9
# ╟─a3569603-9c03-426e-b3a2-f44f90bc65fe
# ╟─514886a7-9512-4f32-b26d-63c27f010bc8
# ╟─1708f378-89b2-4d93-bc22-776bcf93ae98
