### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 27624a18-f031-11eb-36e5-97b5b4ccc70e
begin
	using Pkg; Pkg.activate()
	using DrWatson

md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ c86794a4-c4bf-4685-adcb-f079a4cc7514
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using NCDatasets
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ acb0f419-b3d3-47e5-b91d-fd17ebf4502b
md"
# 01a. Defining the Deep Tropics

We mostly define the tropics to be around 30º in latitude from the equator.  However, here we have a look at the temperature profiles, buoyancy temperature, moist static energy and (maybe) saturation moist entropy at 500 hPa in order to determine the latitudes between which the Weak Temperature Gradient hold most strongly, and then define that region as the Deep Tropics.
"

# ╔═╡ c96af41c-c02e-4e66-8c08-0512dd9c4800
md"
### A. Loading Data and Averaging ...
"

# ╔═╡ e200ea73-16a8-447f-96c4-f06f7361d8e2
begin
	tds = NCDataset(datadir("reanalysis","buoyancy","era5-GLBx0.25-b_sfc-1979.nc"))
	lon = tds["longitude"][:]; nlon = length(lon)
	lat = tds["latitude"][:];  nlat = length(lat)
	close(tds)	
	md"Finding dimensions and preallocating arrays ..."
end

# ╔═╡ f0fcf809-204f-4290-acf0-5f1a04fcefac
begin
	z    = zeros(nlon,nlat,12)
	q    = zeros(nlon,nlat,12)
	tair = zeros(nlon,nlat,12)
	ciwc = zeros(nlon,nlat,12)
	clwc = zeros(nlon,nlat,12)
	crwc = zeros(nlon,nlat,12)
	cswc = zeros(nlon,nlat,12)
	t2m  = zeros(nlon,nlat,12)
	psfc = zeros(nlon,nlat,12)
	
	for iyr = 1979 : 2019
		tfol = datadir("reanalysis","buoyancy")
		sds  = NCDataset(joinpath(tfol,"era5-GLBx0.25-b_sfc-1979.nc"))
		pds  = NCDataset(joinpath(tfol,"era5-GLBx0.25-b_air-500hPa-1979.nc"))
		z[:,:,:]    += pds["z"][:]
		q[:,:,:]    += pds["q"][:]
		tair[:,:,:] += pds["t"][:]
		ciwc[:,:,:] += pds["ciwc"][:]
		clwc[:,:,:] += pds["clwc"][:]
		crwc[:,:,:] += pds["crwc"][:]
		cswc[:,:,:] += pds["cswc"][:]
		t2m[:,:,:]  += sds["t2m"][:]
		psfc[:,:,:] += sds["sp"][:]
		close(sds); close(pds)
	end
	
	z    = dropdims(mean(z   /41,dims=3),dims=3)
	q    = dropdims(mean(q   /41,dims=3),dims=3)
	tair = dropdims(mean(tair/41,dims=3),dims=3)
	ciwc = dropdims(mean(ciwc/41,dims=3),dims=3)
	clwc = dropdims(mean(clwc/41,dims=3),dims=3)
	crwc = dropdims(mean(crwc/41,dims=3),dims=3)
	cswc = dropdims(mean(cswc/41,dims=3),dims=3)
	t2m  = dropdims(mean(t2m /41,dims=3),dims=3)
	psfc = dropdims(mean(psfc/41,dims=3),dims=3)
	
	md"Loading and averaging data ..."
end

# ╔═╡ 98e126a2-87cb-4e33-b5e0-1094991ceaf3
begin
	crd = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = crd[:,1]; y = crd[:,2];
md"Loading coastlines ..."
end

# ╔═╡ e8a9291f-eda2-43ad-8c81-b9f56ffbaed3
md"
### B. Calculating some Metrics ...

We use the following metrics, to see if they can be used as a guide to see where the WTG approximation holds most strongly, and therefore define the Deep Tropics (DTP) GeoRegion that we will use in the rest of this project.
"

# ╔═╡ 4a5625fa-706b-4b3f-a607-9ee5b8271818
md"
#### i. Temperature @ 500 hPa

As a control, we investigate the temperature at 500 hPa and see if the tropics are in WTG balance in this region ...
"

# ╔═╡ 8c76ee66-c9d4-4c61-9102-0342525cf54b
begin
	pplt.close(); fta,ata = pplt.subplots(aspect=2*4/pi,axwidth=5)
	
	cta = ata[1].contourf(lon,lat,tair',extend="both",levels=260:0.5:270)
	ata[1].plot(x,y,lw=0.5,c="k")
	ata[1].plot([0,360],[10,10],lw=1,c="k",linestyle="--")
	ata[1].plot([0,360],[-10,-10],lw=1,c="k",linestyle="--")
	ata[1].format(
		xlim=(0,360),xlocator=0:60:360,
		ylim=(-90,90),yscale="sine",ylocator=-60:15:60,
		ylabel=L"Latitude / $\degree$",suptitle="Air Temperature at 500 hPa",
		ultitle="(a)"
	)
	
	tairzon = mean(tair,dims=1)
	slat1 =  sind.(lat)
	slat2 = (slat1[2:end].+slat1[1:(end-1)])/2
	slat3 =  slat1[2:end].-slat1[1:(end-1)]
	tairzong = (tairzon[2:end].-tairzon[1:(end-1)]) ./ slat3
	
	p1_1 = ata[1].panel("r",space=0.1,width="4em")
	p1_1.plot(dropdims(tairzon,dims=1),lat)
	p1_1.plot(tairzong/5 .+265,asind.(slat2),linestyle=":")
	p1_1.plot([230,275],[10,10],lw=1,c="k",linestyle="--")
	p1_1.plot([230,275],[-10,-10],lw=1,c="k",linestyle="--")
	p1_1.plot([265,265],[90,-90],lw=1,c="k",linestyle="--")
	p1_1.format(xlim=(255,275),ultitle="(b)")
	
	tairdiff = tair .- tairzon
	tairdiff = tairdiff[:,(lat.<10) .& (lat.>-10)]
	tairdiff = dropdims(mean(tairdiff,dims=2),dims=2)
	
	p1_2 = ata[1].panel("b",space=0.1,width="4em")
	p1_2.plot(lon,tairdiff,c="k")
	p1_2.plot(lon,0.5*sind.(lon.-105),c="k",lw=0.5,linestyle="--")
	p1_2.format(xlabel=L"Longitude / $\degree$",ultitle="(c)")
	
	ata[1].colorbar(cta,loc="r",label="K")
	
	fta.savefig(plotsdir("01a-buoyancy-t_air-500hPa.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("01a-buoyancy-t_air-500hPa.png"))
end

# ╔═╡ 38d5e1f1-03fd-4863-b363-d4ad0440d6d2
md"
#### ii. Virtual Temperature @ 500 hPa

In the System of Atmospheric Modelling, Peter Blossey used virtual temperature $T_v$ as a proxy to the buoyancy of the atmosphere in order to calculate the gravity-wave adjustment.  The formula for virtual temperature is given as:

$$T_v = T\cdot\frac{(1+q/\varepsilon)}{1+q}$$

where $\varepsilon=0.622$ is the ratio of the gas constants of air and water vapor, and $q$ is the specific humidity of water vapour.
"

# ╔═╡ 161bd699-0395-4c4e-a144-8a3d91e087bc
begin
	vt  = tair .* (1 .+ q/0.622) ./ (1 .+ q)
	md"Calculating the Virtual Temperature at 500 hPa ..."
end

# ╔═╡ 417d50e0-2487-42da-9bb1-e66f4fb4e766
begin
	pplt.close(); ftv,atv = pplt.subplots(aspect=2*4/pi,axwidth=5)
	
	ctv = atv[1].contourf(lon,lat,vt',extend="both",levels=260:0.5:270)
	atv[1].plot(x,y,lw=0.5,c="k")
	atv[1].plot([0,360],[10,10],lw=1,c="k",linestyle="--")
	atv[1].plot([0,360],[-10,-10],lw=1,c="k",linestyle="--")
	atv[1].format(
		xlim=(0,360),xlocator=0:60:360,
		ylim=(-90,90),yscale="sine",ylocator=-60:15:60,
		ylabel=L"Latitude / $\degree$",suptitle="Virtual Temperature at 500 hPa",
		ultitle="(a)"
	)
	
	vtzon = mean(vt,dims=1)
	vtzong = (vtzon[2:end].-vtzon[1:(end-1)]) ./ slat3
	
	p2_1 = atv[1].panel("r",space=0.1,width="4em")
	p2_1.plot(dropdims(vtzon,dims=1),lat)
	p2_1.plot(vtzong/5 .+265,asind.(slat2),linestyle=":")
	p2_1.plot([230,275],[10,10],lw=1,c="k",linestyle="--")
	p2_1.plot([230,275],[-10,-10],lw=1,c="k",linestyle="--")
	p2_1.plot([265,265],[90,-90],lw=1,c="k",linestyle="--")
	p2_1.format(xlim=(255,275),ultitle="(b)")
	
	vtdiff = vt .- vtzon
	vtdiff = vtdiff[:,(lat.<10) .& (lat.>-10)]
	vtdiff = dropdims(mean(vtdiff,dims=2),dims=2)
	
	p2_2 = atv[1].panel("b",space=0.1,width="4em")
	p2_2.plot(lon,vtdiff,c="k")
	p2_2.plot(lon,0.5*sind.(lon.-90),c="k",lw=0.5,linestyle="--")
	p2_2.format(xlabel=L"Longitude / $\degree$",ultitle="(c)")
	
	atv[1].colorbar(ctv,loc="r",label="K")
	
	ftv.savefig(plotsdir("01a-buoyancy-vt-500hPa.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("01a-buoyancy-vt-500hPa.png"))
end

# ╔═╡ 58f78fc0-481f-42b1-85f6-61cb1fa19b91
md"
#### iii. Moist Static Energy @ 500 hPa

We also tried using Moist Static Energy at the 500 hPa pressure level.

$$s_e = c_pT + \Phi + L_vq$$

where $c_p$ is the heat capacity of dry air, $\Phi$ is the geopotential, $L_v$ is the latent heat of vapourization, and $q$ is the specific humidity of water vapour.
"

# ╔═╡ 14099821-b803-48c7-82ba-58510b1abc5a
begin
	mse = 1005.7 * tair .+ z .+ 2.260e6 * q; mse = mse / 1005.7
	md"Calculating the Moist Static Energy ..."
end

# ╔═╡ 5dd6a0f4-479d-4122-ba9c-f9435cd74d41
begin
	pplt.close(); fmse,amse = pplt.subplots(aspect=2*4/pi,axwidth=5)
	
	cmse = amse[1].contourf(lon,lat,mse',extend="both",levels=305:330)
	amse[1].plot(x,y,lw=0.5,c="k")
	amse[1].plot([0,360],[10,10],lw=1,c="k",linestyle="--")
	amse[1].plot([0,360],[-10,-10],lw=1,c="k",linestyle="--")
	amse[1].format(
		xlim=(0,360),xlocator=0:60:360,
		ylim=(-90,90),yscale="sine",ylocator=-60:15:60,
		ylabel=L"Latitude / $\degree$",suptitle="Moist Static Energy at 500 hPa / K",
		ultitle="(a)"
	)
	
	msezon = mean(mse,dims=1)
	msezong = (msezon[2:end].-msezon[1:(end-1)]) ./ slat3
	
	p3_1 = amse[1].panel("r",space=0.1,width="4em")
	p3_1.plot(dropdims(msezon,dims=1),lat)
	p3_1.plot(msezong/2 .+325,asind.(slat2),linestyle=":")
	p3_1.plot([270,370],[10,10],lw=1,c="k",linestyle="--")
	p3_1.plot([270,370],[-10,-10],lw=1,c="k",linestyle="--")
	p3_1.plot([325,325],[90,-90],lw=1,c="k",linestyle="--")
	p3_1.format(xlim=(290,360),ultitle="(b)")
	
	msediff = mse .- msezon
	msediff = msediff[:,(lat.<10) .& (lat.>-10)]
	msediff = dropdims(mean(msediff,dims=2),dims=2)
	
	p3_2 = amse[1].panel("b",space=0.1,width="4em")
	p3_2.plot(lon,msediff,c="k")
	p3_2.plot(lon,1*sind.(lon.-120).+1.5*sind.(3*lon),c="k",lw=0.5,linestyle="--")
	p3_2.format(xlabel=L"Longitude / $\degree$",ultitle="(c)",ylim=(-3,3))
	
	amse[1].colorbar(cmse,loc="r",locator=305:5:330)
	
	fmse.savefig(plotsdir("01a-buoyancy-mse-500hPa.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("01a-buoyancy-mse-500hPa.png"))
end

# ╔═╡ 11f6a8cd-7f6e-4604-a06b-3b66d394e0d0
md"
### C. Conclusion

It can thus be seen that although we would like to extend the DTP region to 15º or even 20º in latitude, the WTG approximation only really holds up to around 10º off the equator.  This of course, is also dependent on region, with the WTG approximation possibly holding over a larger range of latitudes over the Indo-Pacific warmpool compared to other regions.

We also see that MSE doesn't necessarily give us as sharp a boundary of where the WTG approximation holds as well as virtual temperature or even air temperature.  The MSE seems to be affected much more by the SST temperature distribution, while virtual temperature and air temperature at 500 hPa seem to mostly have smoothed out those variations.
"

# ╔═╡ Cell order:
# ╟─acb0f419-b3d3-47e5-b91d-fd17ebf4502b
# ╟─27624a18-f031-11eb-36e5-97b5b4ccc70e
# ╟─c86794a4-c4bf-4685-adcb-f079a4cc7514
# ╟─c96af41c-c02e-4e66-8c08-0512dd9c4800
# ╟─e200ea73-16a8-447f-96c4-f06f7361d8e2
# ╟─f0fcf809-204f-4290-acf0-5f1a04fcefac
# ╟─98e126a2-87cb-4e33-b5e0-1094991ceaf3
# ╟─e8a9291f-eda2-43ad-8c81-b9f56ffbaed3
# ╟─4a5625fa-706b-4b3f-a607-9ee5b8271818
# ╟─8c76ee66-c9d4-4c61-9102-0342525cf54b
# ╟─38d5e1f1-03fd-4863-b363-d4ad0440d6d2
# ╟─161bd699-0395-4c4e-a144-8a3d91e087bc
# ╟─417d50e0-2487-42da-9bb1-e66f4fb4e766
# ╟─58f78fc0-481f-42b1-85f6-61cb1fa19b91
# ╠═14099821-b803-48c7-82ba-58510b1abc5a
# ╟─5dd6a0f4-479d-4122-ba9c-f9435cd74d41
# ╟─11f6a8cd-7f6e-4604-a06b-3b66d394e0d0
