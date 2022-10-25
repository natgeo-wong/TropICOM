### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ a675084e-f638-11eb-2930-b70006b9445c
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 4e08778c-cd46-4a79-9b6a-96ca6447e975
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using ERA5Reanalysis
	using ImageFiltering
	using NCDatasets

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ faa41f45-5a58-4bcb-a48c-dcd17e8bad3a
md"
# 01c. Land-Sea Mask Filtering
"

# ╔═╡ 75ea588a-64df-4b31-9dae-f99f425ee55a
md"
### A. Loading Land-Sea Mask Data
"

# ╔═╡ c557eeb4-fab4-4751-af4b-bea2115e85a9
lsd = getLandSea(ERA5Region(GeoRegion("GLB"),gres=0.25),path=datadir("emask"))

# ╔═╡ dc6929b0-c95d-4c71-9a4c-cb9fa4b01e34
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]; y = coast[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ a9bae1a3-c994-4b36-b9f8-740aa14a3929
begin
	pplt.close(); f1,a1 = pplt.subplots(aspect=1,axwidth=3)
	
	c = a1[1].pcolormesh(
		lsd.lon,lsd.lat,lsd.lsm',
		levels=(0:10)/10,extend="both",cmap="Delta"
	)
	# a1[1].plot(x,y,c="k",lw=0.5)
	a1[1].plot([95.5,106,105.5,95,95.5],[6.5,-6,-6.5,6,6.5],c="r",lw=1)
	a1[1].plot([95,105.5,105,94.5,95],[6,-6.5,-7,5.5,6],c="r",lw=1)
	a1[1].plot([94.5,105,104.5,94,94.5],[5.5,-7,-7.5,5,5.5],c="r",lw=1)
	a1[1].plot([94,104.5,104,93.5,94],[5,-7.5,-8,4.5,5],c="r",lw=1)
	a1[1].plot([93.5,104,103.5,93,93.5],[4.5,-8,-8.5,4,4.5],c="r",lw=1)
	a1[1].plot([93,103.5,103,92.5,93],[4,-8.5,-9,3.5,4],c="r",lw=1)
	a1[1].format(xlim=(90,110),ylim=(-10,10))
	a1[1].colorbar(c,loc="r")
	
	f1.savefig("testlsm.png",transparent=false,dpi=200)
	PNGFiles.load("testlsm.png")
end

# ╔═╡ 95ffc474-8416-451c-ba9a-4452dc582626
md"
### B. Smoothing the Land-Sea Mask?
"

# ╔═╡ a66e0d88-7f48-4d7d-b716-3a424d93e90b
function filterlsm(olsm;iterations=1,smooth=1)
	
	it = 0
	nlsm = deepcopy(olsm)
	while it < iterations
		nlsm = log10.(imfilter(10. .^nlsm, Kernel.gaussian(smooth),"circular"));
		nlsm = (nlsm.+olsm)/2
		it += 1
	end
	
	nlsm[nlsm.<0] .= minimum(nlsm[nlsm.>0])
	
	return nlsm
	
end

# ╔═╡ d3864cc7-284b-4278-9499-0b6f4d274608
begin
	nlsm1d0 = filterlsm(lsd.lsm,iterations=10,smooth=1)
	nlsm1d5 = filterlsm(lsd.lsm,iterations=10,smooth=1.5)
	nlsm2d0 = filterlsm(lsd.lsm,iterations=10,smooth=2)
	md"Performing gaussian filtering/smoothing on land-sea mask ..."
end

# ╔═╡ 05c13f7e-58f8-4136-830d-fb20a9dc0e6f
begin
	geo = GeoRegion("DTP_IPW")
end

# ╔═╡ 8471606c-fead-450a-80f0-fd4cf308c3e5
begin
	N,S,E,W = geo.N,geo.S,geo.E,geo.W
	ggrd = RegionGrid(geo,lsd.lon,lsd.lat)
	rlsm = extractGrid(lsd.lsm,ggrd)
	flsm1d0 = extractGrid(nlsm1d0,ggrd)
	flsm1d5 = extractGrid(nlsm1d5,ggrd)
	flsm2d0 = extractGrid(nlsm2d0,ggrd)
	md"Extracting information for region ..."
end

# ╔═╡ 51d55158-32df-4135-ab12-186e295ce499
begin
	asp = (E-W+2)/(N-S+2)
	if asp > 3
		freg,areg = pplt.subplots(nrows=4,axwidth=asp*1.2,aspect=asp)
	else
		freg,areg = pplt.subplots(nrows=2,ncols=2,axwidth=asp*1.5,aspect=asp)
	end
	
	lvls = vcat(10. .^(-5:-1),0.2,0.5,0.9,0.95,0.97,0.98,0.99,0.999)

	creg = areg[1].pcolormesh(
		ggrd.lon,ggrd.lat,rlsm',
		levels=lvls,cmap="delta",extend="both"
	)
	areg[1].format(urtitle="Raw")
	
	areg[2].pcolormesh(
		ggrd.lon,ggrd.lat,flsm1d0',
		levels=lvls,cmap="delta",extend="both"
	)
	areg[2].format(urtitle="Smooth = 1")
	
	areg[3].pcolormesh(
		ggrd.lon,ggrd.lat,flsm1d5',
		levels=lvls,cmap="delta",extend="both"
	)
	areg[3].format(urtitle="Smooth = 1.5")
	
	areg[4].pcolormesh(
		ggrd.lon,ggrd.lat,flsm2d0',
		levels=lvls,cmap="delta",extend="both"
	)
	areg[4].format(urtitle="Smooth = 2")

	for ax in areg
		ax.plot(x,y,c="k",lw=0.5)
		ax.format(
			xlim=(ggrd.lon[1].-1,ggrd.lon[end].+1),ylim=(S-1,N+1),
			ylabel=L"Latitude / $\degree$",xlabel=L"Longitude / $\degree$",
			suptitle="Land-Sea Mask",grid=true
		)
	end

	if asp > 3
		freg.colorbar(creg,loc="r",length=0.4)
	else
		freg.colorbar(creg,loc="r",length=0.8)
	end
	freg.savefig(plotsdir("01c-flsm_$(geo.regID).png"),transparent=false,dpi=200)
	load(plotsdir("01c-flsm_$(geo.regID).png"))
end

# ╔═╡ 5c08a214-baaf-4f05-b505-67f49d02b4ef
begin
	fnc = datadir("flsm","flsm-$(geo.regID).nc")
	if !isdir(datadir("flsm")); mkpath(datadir("flsm")) end
	if isfile(fnc)
		rm(fnc,force=true)
	end

	ds = NCDataset(fnc,"c")
	ds.dim["longitude"] = length(ggrd.lon)
	ds.dim["latitude"]  = length(ggrd.lat)

	nclon = defVar(ds,"longitude",Float64,("longitude",),attrib = Dict(
		"units"     => "degrees_east",
		"long_name" => "longitude",
	))

	nclat = defVar(ds,"latitude",Float64,("latitude",),attrib = Dict(
		"units"     => "degrees_north",
		"long_name" => "latitude",
	))

	ncvar = defVar(ds,"flsm",Float64,("longitude","latitude",),attrib = Dict(
		"long_name"     => "land_sea_mask_filtered",
		"full_name"     => "Filtered Land-Sea Mask",
		"units"         => "0-1",
	))

	nclon[:] = ggrd.lon
	nclat[:] = ggrd.lat
	ncvar[:] = flsm1d0

	close(ds)
	md"Saving Filtered Land-Sea Mask ..."
end

# ╔═╡ Cell order:
# ╟─faa41f45-5a58-4bcb-a48c-dcd17e8bad3a
# ╟─a675084e-f638-11eb-2930-b70006b9445c
# ╟─4e08778c-cd46-4a79-9b6a-96ca6447e975
# ╟─75ea588a-64df-4b31-9dae-f99f425ee55a
# ╟─c557eeb4-fab4-4751-af4b-bea2115e85a9
# ╟─dc6929b0-c95d-4c71-9a4c-cb9fa4b01e34
# ╟─a9bae1a3-c994-4b36-b9f8-740aa14a3929
# ╟─95ffc474-8416-451c-ba9a-4452dc582626
# ╠═a66e0d88-7f48-4d7d-b716-3a424d93e90b
# ╠═d3864cc7-284b-4278-9499-0b6f4d274608
# ╟─05c13f7e-58f8-4136-830d-fb20a9dc0e6f
# ╟─8471606c-fead-450a-80f0-fd4cf308c3e5
# ╟─51d55158-32df-4135-ab12-186e295ce499
# ╟─5c08a214-baaf-4f05-b505-67f49d02b4ef
