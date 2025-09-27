### A Pluto.jl notebook ###
# v0.19.11

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
	using DelimitedFiles
	using ERA5Reanalysis
	using Interpolations
	using NCDatasets
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))

md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 90fffbc8-524d-11eb-232a-1bada28d5505
md"
# 3c. Total Cloud Cover Data

The first variable that we want to explore among the ERA5 reanalysis variables is Skin Temperature.  More details about the retrieval of reanalysis data can be found in notebook `03-reanalysis.jl`.
"

# ╔═╡ d82366b0-53b1-11eb-26c1-ff1bb6ccb027
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ c1fafcba-530f-11eb-1cc2-67d10a3fa606
md"
### A. Modelling the Diurnal Cycle of Total Cloud Cover

We wish to find the following characteristics of the diurnal cycle of total cloud cover:
* The mean $\mu$ total cloud cover
* The amplitude $A$ of the diurnal cycle (max-min)/2 of total cloud cover
* The hour $\theta$ at which total cloud cover is a maximum
"

# ╔═╡ 3565af3c-5311-11eb-34c4-2d228b05b17c
md"
ERA5 reanalysis is defined according to UTC.  Therefore, for every grid point, we must also correct to the local time before fitting the data to the model.
"

# ╔═╡ fd344560-94d3-11eb-2b79-05c905c5953f
longitude2timeshift(longitude::Real) = longitude / 180 * 12

# ╔═╡ a6a688ca-53ab-11eb-2776-b5380ffb26c1
function eradiurnal2model(data,longitude)

	nlon,nlat,nt = size(data)
	θ = zeros(nlon,nlat)
	A = zeros(nlon,nlat)
	μ = zeros(nlon,nlat)

	idat = zeros(nt+1)
	it   = 0:0.1:24; it = it[1:(end-1)]; nit = length(it)
	ts   = longitude2timeshift.(longitude)
	var  = zeros(nit)

	for ilat = 1 : nlat, ilon = 1 : nlon

		tl = (0:24) .+ ts[ilon]

        idat[1:24] = data[ilon,ilat,:]
		idat[end]  = data[ilon,ilat,1]

		itp = interpolate(idat,BSpline(Cubic(Periodic(OnGrid()))))
        stp = scale(itp,tl)
        etp = extrapolate(stp,Periodic())
		var[:] = etp[it]

		μ[ilon,ilat] = mean(@view data[ilon,ilat,:])
		A[ilon,ilat] = (maximum(idat) .- minimum(idat))/2
		θ[ilon,ilat] = argmax(var) / nit * 24

    end

	return μ,A,θ

end

# ╔═╡ c5fd453b-27b9-4a9b-9411-fee451920977
md"
### B. Defining GeoRegions
"

# ╔═╡ 5f0fb501-c9b8-44c3-8cb8-876c57f4813d
egeo = ERA5Region(GeoRegion("TRP"))

# ╔═╡ c117a827-ec0b-42f2-8ae3-404751583712
lsd = getLandSea(egeo,path=datadir("emask"))

# ╔═╡ aa05317e-530b-11eb-2ec1-93aff65659dd
md"
### C. Retrieving ERA5 Total Cloud Cover Data
"

# ╔═╡ 103f85e8-530c-11eb-047d-a537aa60075d
function retrieveera()

    ds  = NCDataset(datadir("compiled","era5hr-TRPx0.25-tcc-compiled.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	var = ds["tcc"][:] * 1
	close(ds)

	μ,A,θ = eradiurnal2model(var,lon)

    return lon,lat,μ,A,θ

end

# ╔═╡ e8141e20-53af-11eb-1a23-81d34293c5eb
begin
	lon,lat,μ,A,θ = retrieveera()
md"Modelling diurnal cycle of surface temperature"
end

# ╔═╡ bb59b8d6-53b1-11eb-3631-87ef61219c4c
begin
	pplt.close(); f,axs = pplt.subplots(nrows=3,axwidth=6,aspect=6)

	c = axs[1].contourf(lon,lat,μ',levels=0.1:0.1:0.9,cmap="Blues",extend="both")
	axs[1].plot(x,y,c="k",lw=0.5)
	axs[1].format(urtitle=L"$\mu$")
	axs[1].colorbar(c,loc="r")

	c = axs[2].contourf(lon,lat,A',levels=10. .^(-2:0.2:0),extend="both")
	axs[2].plot(x,y,c="k",lw=0.5)
	axs[2].format(urtitle="A")
	axs[2].colorbar(c,loc="r",ticks=[0.01,0.1,1])

	c = axs[3].pcolormesh(lon,lat,θ',cmap="romaO",levels=0:24,extend="both")
	axs[3].plot(x,y,c="k",lw=0.5)
	axs[3].format(urtitle=L"$\theta$ / Hour of Day")
	axs[3].colorbar(c,loc="r",locator=0:6:24,minorticks=0:24)

	for ax in axs
		ax.format(
			xlim=(0,360),ylim=(-30,30),xlocator=0:60:360,
			xlabel=L"Longitude / $\degree$",
			ylabel=L"Latitude / $\degree$",
		)
	end

	f.savefig(plotsdir("03c-tccspatial_TRP.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("03c-tccspatial_TRP.png"))
end

# ╔═╡ 5c0e5bae-554e-11eb-3f83-a364ae0a2485
md"
Test
"

# ╔═╡ 68cfc46c-5755-11eb-1702-373942539652
md"
### D. Regional Analysis

We can get quick snapshots of the results for different GeoRegions specified in this project.
"

# ╔═╡ 52b39ff8-9426-11eb-2a86-43f7da15f62e
begin
	geo = GeoRegion("SEA")
	md"Defining region coordinates ..."
end

# ╔═╡ ea7f0956-575b-11eb-3e3f-a1ba3e08b771
begin
	N,S,E,W = geo.N,geo.S,geo.E,geo.W
	ggrd = RegionGrid(geo,lon,lat)
	rμ = extractGrid(μ,ggrd)
	rA = extractGrid(A,ggrd)
	rθ = extractGrid(θ,ggrd)
	md"Extracting information for region ..."
end

# ╔═╡ 5714c13c-575c-11eb-06d4-838b4e8dbcd7
begin
	asp = (E-W+2)/(N-S+2)
	pplt.close()
	if asp > 1.5
		freg,areg = pplt.subplots(nrows=3,axwidth=asp*1.2,aspect=asp)
	else
		freg,areg = pplt.subplots(ncols=3,axwidth=2,aspect=asp)
	end

	creg = areg[1].contourf(
		ggrd.lon,ggrd.lat,rμ',
		levels=(1:9)/10,cmap="Blues",extend="both"
	)
	areg[1].plot(x,y,c="k",lw=0.5)
	areg[1].format(urtitle=L"$\mu$")
	areg[1].colorbar(creg,loc="r")

	creg = areg[2].contourf(
		ggrd.lon,ggrd.lat,rA',
		levels=10. .^(-2:0.2:0),extend="both"
	)
	areg[2].plot(x,y,c="k",lw=0.5)
	areg[2].format(urtitle="A")
	areg[2].colorbar(creg,loc="r",ticks=[0.01,0.1,1,10,100])

	creg = areg[3].pcolormesh(ggrd.lon,ggrd.lat,rθ',cmap="romaO",levels=0:24)
	areg[3].plot(x,y,c="k",lw=0.5)
	areg[3].format(urtitle=L"$\theta$ / Hour of Day")
	areg[3].colorbar(creg,loc="r",locator=0:3:24,minorticks=0:24)

	for ax in areg
		ax.format(
			xlim=(ggrd.lon[1].-1,ggrd.lon[end].+1),
			xlabel=L"Longitude / $\degree$",
			ylim=(S-1,N+1),ylabel=L"Latitude / $\degree$",
			grid=true
		)
	end

	freg.savefig(plotsdir("03c-tccspatial_$(geo.regID).png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("03c-tccspatial_$(geo.regID).png"))
end

# ╔═╡ c4792bf2-5552-11eb-3b52-997f59fd42f3
md"
### E. Binning of Precipitation Data

We now wish to bin the modelled precipitation data in order to determine the relationship between precipitation rate and the diurnal cycle over land and sea.  Is it different?
"

# ╔═╡ 3a207827-e423-4c22-b594-28c3b8360461
begin
	lsc = pplt.Colors("Delta_r",15)
	md"Colours for different regions ..."
end

# ╔═╡ f752b054-57c1-11eb-117c-ed52464aa25f
md"
#### i. Mean Total Cloud Cover $\mu$
"

# ╔═╡ 4b289fa8-57b9-11eb-0923-116c3d9444bb
begin
	lbins = collect(0:0.005:1); lpbin = (lbins[2:end].+lbins[1:(end-1)])/2
	lbin_DTP,lavg_DTP = bindatasfclnd(GeoRegion("DTP"),lbins,μ,lsd)
	lbin_SEA,lavg_SEA = bindatasfclnd(GeoRegion("SEA"),lbins,μ,lsd)
	lbin_CRB,lavg_CRB = bindatasfclnd(GeoRegion("CRB"),lbins,μ,lsd)
	lbin_TRA,lavg_TRA = bindatasfclnd(GeoRegion("TRA"),lbins,μ,lsd)
	lbin_AMZ,lavg_AMZ = bindatasfclnd(GeoRegion("AMZ"),lbins,μ,lsd)

	sbins = collect(0:0.005:1); spbin = (sbins[2:end].+sbins[1:(end-1)])/2
	sbin_DTP,savg_DTP = bindatasfcsea(GeoRegion("DTP"),sbins,μ,lsd)
	sbin_SEA,savg_SEA = bindatasfcsea(GeoRegion("SEA"),sbins,μ,lsd)
	sbin_CRB,savg_CRB = bindatasfcsea(GeoRegion("CRB"),sbins,μ,lsd)
	sbin_EPO,savg_EPO = bindatasfcsea(GeoRegion("AR6_EPO"),sbins,μ,lsd)
	sbin_EIO,savg_EIO = bindatasfcsea(GeoRegion("AR6_EIO"),sbins,μ,lsd)
	sbin_EAO,savg_EAO = bindatasfcsea(GeoRegion("AR6_EAO"),sbins,μ,lsd)

	md"Binning mean skin temperatures for different tropical regions ..."
end

# ╔═╡ e7ff7ec8-57b9-11eb-0115-abbe4aa9a1a9
begin
	pplt.close(); fbin,abin = pplt.subplots(ncols=2,aspect=2,axwidth=3);

	lgd = Dict("ncol"=>1,"frame"=>false)
	abin[1].plot(lpbin,lbin_DTP,c="k",lw=0.5);
	abin[1].plot(lpbin,lbin_CRB,c=lsc[10],lw=0.5)
	abin[1].plot(lpbin,lbin_SEA,c=lsc[5],lw=0.5)
	abin[1].plot(lpbin,lbin_AMZ,c=lsc[4],lw=0.5)
	abin[1].plot(lpbin,lbin_TRA,c=lsc[3],lw=0.5)

	abin[1].plot([1,1]*lavg_DTP,[0.1,50],c="k")
	abin[1].plot([1,1]*lavg_CRB,[0.1,50],c=lsc[10])
	abin[1].plot([1,1]*lavg_SEA,[0.1,50],c=lsc[5])
	abin[1].plot([1,1]*lavg_AMZ,[0.1,50],c=lsc[4])
	abin[1].plot([1,1]*lavg_TRA,[0.1,50],c=lsc[3])

	abin[2].plot(spbin,sbin_DTP,c="k",lw=0.5);
	abin[2].plot(spbin,sbin_EPO,c=lsc[13],lw=0.5);
	abin[2].plot(spbin,sbin_EIO,c=lsc[12],lw=0.5);
	abin[2].plot(spbin,sbin_EAO,c=lsc[11],lw=0.5);
	abin[2].plot(spbin,sbin_CRB,c=lsc[10],lw=0.5);
	abin[2].plot(spbin,sbin_SEA,c=lsc[5],lw=0.5);

	abin[2].plot([1,1]*savg_DTP,[0.1,50],c="k",label="DTP",legend="r",legend_kw=lgd)
	abin[2].plot([1,1]*savg_EPO,[0.1,50],c=lsc[13],label="AR6_EPO",legend="r")
	abin[2].plot([1,1]*savg_EIO,[0.1,50],c=lsc[12],label="AR6_EIO",legend="r")
	abin[2].plot([1,1]*savg_EAO,[0.1,50],c=lsc[11],label="AR6_EAO",legend="r")
	abin[2].plot([1,1]*savg_CRB,[0.1,50],c=lsc[10],label="CRB",legend="r")
	abin[2].plot([1,1]*savg_SEA,[0.1,50],c=lsc[5],label="SEA",legend="r")
	abin[2].plot([1,1]*NaN,[0.1,50],c=lsc[4],label="AMZ",legend="r")
	abin[2].plot([1,1]*NaN,[0.1,50],c=lsc[3],label="TRA",legend="r")

	abin[1].format(
		xlim=(minimum(lbins),maximum(lbins)),
		ylim=(0,25),#yscale="log",
		ylabel="Density",
		ltitle="(a) Land"
	)

	abin[2].format(
		xlim=(minimum(sbins),maximum(sbins)),
		# ylim=(0.1,30),yscale="log",
		xlabel="Total Cloud Cover Fraction",
		ltitle="(b) Ocean"
	)

	fbin.savefig(plotsdir("03c-tccmean.png"),transparent=false,dpi=200)
	load(plotsdir("03c-tccmean.png"))
end

# ╔═╡ 0fbb0b46-57c2-11eb-365a-a73a2ebda8e4
md"
#### ii. Diurnal Amplitude $A$ vs $\mu$
"

# ╔═╡ 252508a8-57c2-11eb-08b5-8fa673b1ac8a
begin
	lvec = collect(0:0.001:0.2); lAbin = (lvec[2:end].+lvec[1:(end-1)])/2
	lAbin_DTP,lAavg_DTP = bindatasfclnd(GeoRegion("DTP"),lvec,A,lsd)
	lAbin_SEA,lAavg_SEA = bindatasfclnd(GeoRegion("SEA"),lvec,A,lsd)
	lAbin_CRB,lAavg_CRB = bindatasfclnd(GeoRegion("CRB"),lvec,A,lsd)
	lAbin_TRA,lAavg_TRA = bindatasfclnd(GeoRegion("TRA"),lvec,A,lsd)
	lAbin_AMZ,lAavg_AMZ = bindatasfclnd(GeoRegion("AMZ"),lvec,A,lsd)

	svec = collect(0:0.001:0.2); sAbin = (svec[2:end].+svec[1:(end-1)])/2
	sAbin_DTP,sAavg_DTP = bindatasfcsea(GeoRegion("DTP"),svec,A,lsd)
	sAbin_SEA,sAavg_SEA = bindatasfcsea(GeoRegion("SEA"),svec,A,lsd)
	sAbin_CRB,sAavg_CRB = bindatasfcsea(GeoRegion("CRB"),svec,A,lsd)
	sAbin_EPO,sAavg_EPO = bindatasfcsea(GeoRegion("AR6_EPO"),svec,A,lsd)
	sAbin_EIO,sAavg_EIO = bindatasfcsea(GeoRegion("AR6_EIO"),svec,A,lsd)
	sAbin_EAO,sAavg_EAO = bindatasfcsea(GeoRegion("AR6_EAO"),svec,A,lsd)

	md"Binning amplitude of the diurnal cycle for skin temperatures in different tropical regions ..."
end

# ╔═╡ 5f58ae9c-57c2-11eb-1f04-2ddbaf2b4f1b
begin
	pplt.close(); fA,aA = pplt.subplots(ncols=2,aspect=2,axwidth=3);

	aA[1].plot(lAbin,lAbin_DTP,c="k",lw=0.5)
	aA[1].plot(lAbin,lAbin_SEA,c=lsc[10],lw=0.5)
	aA[1].plot(lAbin,lAbin_CRB,c=lsc[5],lw=0.5)
	aA[1].plot(lAbin,lAbin_AMZ,c=lsc[4],lw=0.5)
	aA[1].plot(lAbin,lAbin_TRA,c=lsc[3],lw=0.5)

	aA[1].plot([1,1]*lAavg_DTP,[0.1,50],c="k")
	aA[1].plot([1,1]*lAavg_SEA,[0.1,50],c=lsc[10])
	aA[1].plot([1,1]*lAavg_CRB,[0.1,50],c=lsc[5])
	aA[1].plot([1,1]*lAavg_AMZ,[0.1,50],c=lsc[4])
	aA[1].plot([1,1]*lAavg_TRA,[0.1,50],c=lsc[3])

	aA[2].plot(sAbin,sAbin_DTP,c="k",lw=0.5);
	aA[2].plot(sAbin,sAbin_EPO,c=lsc[13],lw=0.5);
	aA[2].plot(sAbin,sAbin_EIO,c=lsc[12],lw=0.5);
	aA[2].plot(sAbin,sAbin_EAO,c=lsc[11],lw=0.5);
	aA[2].plot(sAbin,sAbin_CRB,c=lsc[10],lw=0.5);
	aA[2].plot(sAbin,sAbin_SEA,c=lsc[5],lw=0.5);

	aA[2].plot([1,1]*sAavg_DTP,[0.1,50],c="k",label="DTP",legend="r",legend_kw=lgd)
	aA[2].plot([1,1]*sAavg_EPO,[0.1,50],c=lsc[13],label="AR6_EPO",legend="r")
	aA[2].plot([1,1]*sAavg_EIO,[0.1,50],c=lsc[12],label="AR6_EIO",legend="r")
	aA[2].plot([1,1]*sAavg_EAO,[0.1,50],c=lsc[11],label="AR6_EAO",legend="r")
	aA[2].plot([1,1]*sAavg_CRB,[0.1,50],c=lsc[10],label="CRB",legend="r")
	aA[2].plot([1,1]*sAavg_SEA,[0.1,50],c=lsc[5],label="SEA",legend="r")
	aA[2].plot([1,1]*NaN,[0.1,50],c=lsc[4],label="TRA",legend="r")
	aA[2].plot([1,1]*NaN,[0.1,50],c=lsc[3],label="AMZ",legend="r")

	aA[1].format(
		xlim=(minimum(lvec),maximum(lvec)),
		ylim=(0,10),#yscale="log",
		ylabel="Density",
		ltitle="(a) Land"
	)

	aA[2].format(
		xlim=(minimum(svec),maximum(svec)),
		ylim=(0.1,15),#yscale="log",
		xlabel="A",
		ltitle="(b) Ocean"
	)

	fA.savefig(plotsdir("03c-tccdiurnalamplitude.png"),transparent=false,dpi=200)
	load(plotsdir("03c-tccdiurnalamplitude.png"))
end

# ╔═╡ 1432fa12-57c7-11eb-0606-7be0389e8fb3
md"
#### iii. Phase Shift $\theta$
"

# ╔═╡ 76627730-57c7-11eb-2037-3f608e085a04
begin
	θvec = collect(0:0.5:24); pθbin = (θvec[2:end].+θvec[1:(end-1)])/24*pi
	pθbin = vcat(pθbin,pθbin[1]+2*pi)

	lθbin_DTP,lθavg_DTP = bindatasfclnd(GeoRegion("DTP"),θvec,θ,lsd)
	lθbin_SEA,lθavg_SEA = bindatasfclnd(GeoRegion("SEA"),θvec,θ,lsd)
	lθbin_CRB,lθavg_CRB = bindatasfclnd(GeoRegion("CRB"),θvec,θ,lsd)
	lθbin_TRA,lθavg_TRA = bindatasfclnd(GeoRegion("TRA"),θvec,θ,lsd)
	lθbin_AMZ,lθavg_AMZ = bindatasfclnd(GeoRegion("AMZ"),θvec,θ,lsd)

	sθbin_DTP,sθavg_DTP = bindatasfcsea(GeoRegion("DTP"),θvec,θ,lsd)
	sθbin_EPO,sθavg_EPO = bindatasfcsea(GeoRegion("AR6_EPO"),θvec,θ,lsd)
	sθbin_EIO,sθavg_EIO = bindatasfcsea(GeoRegion("AR6_EIO"),θvec,θ,lsd)
	sθbin_EAO,sθavg_EAO = bindatasfcsea(GeoRegion("AR6_EAO"),θvec,θ,lsd)
	sθbin_CRB,sθavg_CRB = bindatasfcsea(GeoRegion("CRB"),θvec,θ,lsd)
	sθbin_SEA,sθavg_SEA = bindatasfcsea(GeoRegion("SEA"),θvec,θ,lsd)

	md"Binning hour of maximum skin temperature for different tropical regions ..."
end

# ╔═╡ 8d739d0a-57c7-11eb-16b6-736f595e329e
begin
	pplt.close(); fθ,aθ = pplt.subplots(ncols=2,proj="polar");

	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_DTP,lθbin_DTP[1])),c="k")
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_CRB,lθbin_CRB[1])),c=lsc[10])
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_SEA,lθbin_SEA[1])),c=lsc[5])
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_TRA,lθbin_TRA[1])),c=lsc[4])
	aθ[1].plot(pθbin,sqrt.(vcat(lθbin_AMZ,lθbin_AMZ[1])),c=lsc[3])

	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_DTP,sθbin_DTP[1])),c="k",label="DTP",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_EPO,sθbin_EPO[1])),c=lsc[13],label="AR6_EPO",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_EIO,sθbin_EIO[1])),c=lsc[12],label="AR6_EIO",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_EAO,sθbin_EAO[1])),c=lsc[11],label="AR6_EAO",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_CRB,sθbin_CRB[1])),c=lsc[10],label="CRB",legend="r");
	aθ[2].plot(pθbin,sqrt.(vcat(sθbin_SEA,sθbin_SEA[1])),c=lsc[5],label="SEA",legend="r");
	aθ[2].plot(pθbin,pθbin*NaN,c=lsc[4],label="TRA",legend="r",legend_kw=lgd)
	aθ[2].plot(pθbin,pθbin*NaN,c=lsc[3],label="AMZ",legend="r")

	aθ[1].format(theta0="N",thetaformatter="tau",ltitle="(a) Land")
	aθ[2].format(theta0="N",thetaformatter="tau",ltitle="(b) Ocean")
	aθ[1].format(suptitle=L"$\theta$ / Fraction of Day")

	fθ.savefig(plotsdir("03c-tccdiurnalphase.png"),transparent=false,dpi=200)
	load(plotsdir("03c-tccdiurnalphase.png"))
end

# ╔═╡ Cell order:
# ╟─90fffbc8-524d-11eb-232a-1bada28d5505
# ╟─6f00b8fc-530c-11eb-2242-99d8544f6e14
# ╟─8f30c56c-530c-11eb-2782-33f3c4ed9e89
# ╟─d82366b0-53b1-11eb-26c1-ff1bb6ccb027
# ╟─c1fafcba-530f-11eb-1cc2-67d10a3fa606
# ╟─3565af3c-5311-11eb-34c4-2d228b05b17c
# ╠═fd344560-94d3-11eb-2b79-05c905c5953f
# ╠═a6a688ca-53ab-11eb-2776-b5380ffb26c1
# ╟─c5fd453b-27b9-4a9b-9411-fee451920977
# ╟─5f0fb501-c9b8-44c3-8cb8-876c57f4813d
# ╟─c117a827-ec0b-42f2-8ae3-404751583712
# ╟─aa05317e-530b-11eb-2ec1-93aff65659dd
# ╠═103f85e8-530c-11eb-047d-a537aa60075d
# ╟─e8141e20-53af-11eb-1a23-81d34293c5eb
# ╟─bb59b8d6-53b1-11eb-3631-87ef61219c4c
# ╟─5c0e5bae-554e-11eb-3f83-a364ae0a2485
# ╟─68cfc46c-5755-11eb-1702-373942539652
# ╠═52b39ff8-9426-11eb-2a86-43f7da15f62e
# ╟─ea7f0956-575b-11eb-3e3f-a1ba3e08b771
# ╟─5714c13c-575c-11eb-06d4-838b4e8dbcd7
# ╟─c4792bf2-5552-11eb-3b52-997f59fd42f3
# ╟─3a207827-e423-4c22-b594-28c3b8360461
# ╟─f752b054-57c1-11eb-117c-ed52464aa25f
# ╟─4b289fa8-57b9-11eb-0923-116c3d9444bb
# ╟─e7ff7ec8-57b9-11eb-0115-abbe4aa9a1a9
# ╟─0fbb0b46-57c2-11eb-365a-a73a2ebda8e4
# ╟─252508a8-57c2-11eb-08b5-8fa673b1ac8a
# ╟─5f58ae9c-57c2-11eb-1f04-2ddbaf2b4f1b
# ╟─1432fa12-57c7-11eb-0606-7be0389e8fb3
# ╟─76627730-57c7-11eb-2037-3f608e085a04
# ╟─8d739d0a-57c7-11eb-16b6-736f595e329e
