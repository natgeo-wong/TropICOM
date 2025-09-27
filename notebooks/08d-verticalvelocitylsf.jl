### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

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
	using DelimitedFiles
	using ERA5Reanalysis
	using NCDatasets
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
# 08c. Vertical Profiles binned by Column Relative Humidity

In this notebooks, we bin the vertical profiles of air temperature, relative humidity and atmospheric winds by column relative humidity.
"

# ╔═╡ 455040cc-1b72-4482-a7b6-4ed0cd5125b3
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	clon = coord[:,1]; clat = coord[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ 2f1ee85d-cef1-4fd2-ba06-7c5ec1a3d4a5
md"
### A. Loading the binned 3D Profiles

We have already binned 3D profiles in notebook 08c, and now we want to load the data again.
"

# ╔═╡ 4541f9c7-60f7-43fe-8a18-8df51e13e13f
begin
	ds  = NCDataset(datadir("SAMvertprofile-control.nc"))
	bin = ds["bin"][:]
	z   = ds["z"][:]
	lvl = ds["level"][:]
	pp  = ds["pp"][:] / 100
	w   = ds["w"][:]
	close(ds)
	md"Loading the binned 3D profiles"
end

# ╔═╡ cc64c582-1271-4d18-93f0-efb32edb2a89
begin
	pplt.close(); f1,a1 = pplt.subplots(aspect=3,axwidth=5)
	
	c = a1[1].pcolormesh(bin,lvl,w',levels=vcat((-5:0.5:0)/100,0.5:0.5:5),extend="both")
	a1[1].format(yscale="log")

	f1.colorbar(c)
	f1.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ 4b35c8bc-4978-47ed-99a5-a310e5c1770e
md"
### B. Checking ERA5 Column Saturation Fraction
"

# ╔═╡ 1ae9e619-8a2f-4daa-b736-f45f2fd34fa1
e5ds = ERA5Hourly(start=Date(2021),stop=Date(2021),path=datadir())

# ╔═╡ 92555030-fba0-4a9d-9e93-c80900ffb43a
egeo = ERA5Region("DTP_IPW",gres=0.25)

# ╔═╡ 1faf29c3-6e43-4ced-9357-90fc5ca0abf3
evar = SingleVariable("csf")

# ╔═╡ 1477fc36-94a2-45ad-8dd7-490d42f6d323
lsd = getLandSea(egeo,path=datadir("emask"))

# ╔═╡ 991c62c0-1481-4f60-b0db-66371b2c3c01
begin
	eds = read(e5ds,evar,egeo,Date(2021))
	csf = eds["csf"][:]
	close(eds)
end

# ╔═╡ 8eff5641-0f6b-4f54-9758-5b5e4d4c1a0f
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=4.5,axwidth=6)
	
	c2 = a2[1].pcolormesh(lsd.lon,lsd.lat,csf[:,:,154]',cmap="drywet",levels=vcat(4:0.5:9)/10,extend="both")
	a2[1].plot(clon,clat,c="k",lw=1)
	a2[1].format(
		xlim=(minimum(lsd.lon),maximum(lsd.lon)),
		ylim=(minimum(lsd.lat),maximum(lsd.lat))
	)

	f2.colorbar(c2)
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ eb43ae4f-dae2-4243-a96e-ef985ac06abd
begin
	flsm_ds = NCDataset(datadir("flsm","flsm-$(egeo.geo.regID).nc"))
	flsm = flsm_ds["flsm"][:]
	iiocn = flsm .< 1e-5; iiocn = iiocn[:]
	close(flsm_ds)
	md"Loading filtered land-sea mask ..."
end

# ╔═╡ d221077a-3d16-4786-80a1-46c50844ec68
md"
### C. Binning the ERA5 vertical profile of vertical velocity
"

# ╔═╡ 8e48cc21-503a-4b23-9548-7b0900223d96
begin
	bin_space = 1
	csf_bin = ((-bin_space/2):bin_space:(100+bin_space/2)) ./ 100
end

# ╔═╡ d4c8309a-ac87-4df1-92e1-4e4b2ed32e06
begin
	csf_ocn = reshape(csf,length(lsd.lon)*length(lsd.lat),:)
	csf_ocn = @view csf_ocn[iiocn,:]
	md"Extracting the csf values for gridpoints over the ocean ..."
end

# ╔═╡ a3b779ef-628f-4e9f-a98b-67248b5df3d5
begin
	p1 = era5Pressures(); p1 = p1[p1.>=10]
	ωbin = zeros(length(csf_bin)-1,length(p1))
	
	for ip = 1 : length(p1)
	
		evar_w = PressureVariable("w",hPa=p1[ip])
		wds  = read(e5ds,evar_w,egeo,Date(2021))
		wair = nomissing(wds["w"][:],NaN)
		close(wds)
	
		wair = reshape(wair,length(lsd.lon)*length(lsd.lat),:)
		wair = @view wair[iiocn,:]
	
		for ibin = 1 : (length(csf_bin)-1)
			ii = (csf_ocn .< csf_bin[ibin+1]) .& (csf_ocn .> csf_bin[ibin])
			ωbin[ibin,ip] = mean(wair[ii])
		end
		
	end
	md"Binning the pressure velocity ω by csf values ..."
end

# ╔═╡ d3094db7-65cd-4aa2-b9c8-353dec7ba109
begin
	p2 = era5Pressures(); p2 = p2[p2.>=7]
	zbin = zeros(length(csf_bin)-1,length(p2)+1)
	
	for ip = 1 : length(p2)
	
		evar_z = PressureVariable("z",hPa=p2[ip])
		zds  = read(e5ds,evar_z,egeo,Date(2021))
		zair = nomissing(zds["z"][:],NaN)
		close(zds)
	
		zair = reshape(zair,length(lsd.lon)*length(lsd.lat),:)
		zair = @view zair[iiocn,:]
	
		for ibin = 1 : (length(csf_bin)-1)
			ii = (csf_ocn .< csf_bin[ibin+1]) .& (csf_ocn .> csf_bin[ibin])
			zbin[ibin,ip] = mean(zair[ii])
		end
		
	end
	md"Binning the orographic height by csf values ..."
end

# ╔═╡ 54bc5966-4372-4455-8cbe-7f711c1f7518
begin
	pplt.close(); f3,a3 = pplt.subplots(aspect=3,axwidth=4)
	
	c3 = a3[1].pcolormesh(0:bin_space:100,p1,ωbin',levels=vcat(-5:-0.5,(0:5)/100),extend="both")
	a3[1].format(ylim=(1000,10),yscale="log")
	a3[1].colorbar(c3)
	
	f3.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ c9600307-b5c6-4ba0-bafc-ee7cc7d60ce9
md"
### D. Converting Pressure Velocity to Vertical Velocity
"

# ╔═╡ 6edfa508-c69c-4ef4-b9f7-83f5cec1a272
begin
	psfc_bin = zeros(length(csf_bin)-1)
	
	evar_psfc = SingleVariable("sp")
	psfc_ds   = read(e5ds,evar_psfc,egeo,Date(2021))
	psfc = nomissing(psfc_ds["sp"][:],NaN)
	close(psfc_ds)
	
	psfc = reshape(psfc,length(lsd.lon)*length(lsd.lat),:)
	psfc = @view psfc[iiocn,:]
	
	for ibin = 1 : (length(csf_bin)-1)
		ii = (csf_ocn .< csf_bin[ibin+1]) .& (csf_ocn .> csf_bin[ibin])
		psfc_bin[ibin] = mean(psfc[ii])
	end
	
	md"Binning the surface pressure by csf values ..."
end

# ╔═╡ 5cceaf61-2099-4556-bdfc-cd190b8777c4
begin
	freq_bin = zeros(length(csf_bin)-1)
	
	for ibin = 1 : (length(csf_bin)-1)
		ii = (csf_ocn .< csf_bin[ibin+1]) .& (csf_ocn .> csf_bin[ibin])
		freq_bin[ibin] = sum(ii)
	end
	
	md"Binning the surface pressure by csf values ..."
end

# ╔═╡ 3cdf7150-4e8d-4196-862d-cbd8963f4011
begin
	dzdp_bin = zeros(length(csf_bin)-1,length(p1))
	pair_bin = zeros(length(csf_bin)-1,length(p2)+1)
	pair_bin[:,1:(end-1)] .= log10.(p2'.*100)
	pair_bin[:,end] .= log10.(psfc_bin)
	md"Creating baseline for interpolations ..."
end

# ╔═╡ 96311c7f-dc8f-468a-b138-425f73dd2247
begin
	for ibin = 1 : (length(csf_bin)-1)
	
		if !isnan(pair_bin[ibin,end])
	
			pairii = @view pair_bin[ibin,:]
			zairii = @view zbin[ibin,:]
			spl = Spline1D(pairii,zairii,k=2)
	
			for ip = 1 : length(p1)
	
				ziip = evaluate(spl,pairii[ip+1].+[-0.00001,0.00001])
				dzdp_bin[ibin,ip] = (ziip[2] - ziip[1]) / 9.81 / 0.00002 / 10. .^pairii[ip+1]
			
			end
			
		end
		
	end
	md"Converting pressure velocities to vertical velocities ..."
end

# ╔═╡ cf425cf1-c6ad-4531-8897-7f3b1b6a3bb8
begin
	pplt.close(); f4,a4 = pplt.subplots(aspect=3,axwidth=5)
	
	c4 = a4[1].contourf(
		0:bin_space:100,p1,(ωbin.*dzdp_bin)',
		levels=[-5,-3.16,-2,-1.41,-1,-0.707,-0.5,0,5,10,20,50,100,200,500]/100,
		extend="both"
	)
	a4[1].plot([0,100],[316,316])
	a4[1].format(ylim=(1000,10),yscale="log",ylabel="Pressure / hPa",xlabel="Column Relative Humidity / %")
	a4[1].colorbar(c4,label=L"w / m s$^{-1}$")

	p4 = a4[1].panel("b")
	p4.plot(0:bin_space:100,freq_bin/sum(freq_bin)*100)
	p4.format(ylabel="Density",ylim=(0,3))
	
	f4.savefig(plotsdir("08d-era5binw.png"),transparent=false,dpi=200)
	load(plotsdir("08d-era5binw.png"))
end

# ╔═╡ 418a7423-be7e-483e-9d2e-c2e98264030a
begin
	if isfile(datadir("era5vertprofile.nc"))
		rm(datadir("era5vertprofile.nc"),force=true)
	end
	bds = NCDataset(datadir("era5vertprofile.nc"),"c")

	bds.dim["bin"] = length(csf_bin) - 1
	bds.dim["levels"] = length(p1)

	ncbin = defVar(bds,"bin",Float32,("bin",),attrib=Dict(
		"units"     => "%",
		"long_name" => "column_mean_relative_humidity"
	))

	ncfrq = defVar(bds,"bin_frequency",Float32,("bin",),attrib=Dict(
		"units"     => "%",
		"long_name" => "column_mean_relative_humidity"
	))

	nclvl = defVar(bds,"level",Float32,("levels",),attrib=Dict(
		"units"     => "hPa",
		"long_name" => "pressure_level"
	))

	ncz = defVar(bds,"z",Float32,("bin","levels",),attrib=Dict(
		"units"     => "m",
		"long_name" => "height"
	))

	ncwa = defVar(bds,"w",Float32,("bin","levels"),attrib=Dict(
		"units"     => "m s**-1",
		"long_name" => "vertical_velocity"
	))

	ncbin[:] = collect(0:100)
	ncfrq[:] = freq_bin
	nclvl[:] = p1
	ncz[:]   = zbin[:,2:(end-1)]
	ncwa[:]  = ωbin.*dzdp_bin

	close(ds)
end

# ╔═╡ Cell order:
# ╟─e0fee7ac-fc9b-11ec-0741-9f810d7bcd4c
# ╟─8f2172e0-5b81-420a-8301-dbb8ed0c290a
# ╟─bb51d00a-ee75-4bd8-9744-381b87c6441b
# ╟─455040cc-1b72-4482-a7b6-4ed0cd5125b3
# ╟─2f1ee85d-cef1-4fd2-ba06-7c5ec1a3d4a5
# ╟─4541f9c7-60f7-43fe-8a18-8df51e13e13f
# ╟─cc64c582-1271-4d18-93f0-efb32edb2a89
# ╟─4b35c8bc-4978-47ed-99a5-a310e5c1770e
# ╟─1ae9e619-8a2f-4daa-b736-f45f2fd34fa1
# ╟─92555030-fba0-4a9d-9e93-c80900ffb43a
# ╟─1faf29c3-6e43-4ced-9357-90fc5ca0abf3
# ╟─1477fc36-94a2-45ad-8dd7-490d42f6d323
# ╟─991c62c0-1481-4f60-b0db-66371b2c3c01
# ╟─8eff5641-0f6b-4f54-9758-5b5e4d4c1a0f
# ╟─eb43ae4f-dae2-4243-a96e-ef985ac06abd
# ╟─d221077a-3d16-4786-80a1-46c50844ec68
# ╟─8e48cc21-503a-4b23-9548-7b0900223d96
# ╟─d4c8309a-ac87-4df1-92e1-4e4b2ed32e06
# ╟─a3b779ef-628f-4e9f-a98b-67248b5df3d5
# ╟─d3094db7-65cd-4aa2-b9c8-353dec7ba109
# ╟─54bc5966-4372-4455-8cbe-7f711c1f7518
# ╟─c9600307-b5c6-4ba0-bafc-ee7cc7d60ce9
# ╟─6edfa508-c69c-4ef4-b9f7-83f5cec1a272
# ╟─5cceaf61-2099-4556-bdfc-cd190b8777c4
# ╟─3cdf7150-4e8d-4196-862d-cbd8963f4011
# ╟─96311c7f-dc8f-468a-b138-425f73dd2247
# ╟─cf425cf1-c6ad-4531-8897-7f3b1b6a3bb8
# ╟─418a7423-be7e-483e-9d2e-c2e98264030a
