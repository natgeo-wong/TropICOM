### A Pluto.jl notebook ###
# v0.19.13

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
# 3a. Skin Temperature Data

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
### A. Modelling the Diurnal Cycle of Skin Temperature

We wish to find the following characteristics of the diurnal cycle of skin temperature:
* The mean $\mu$ skin temperature
* The amplitude $A$ of the diurnal cycle (max-min)/2 of skin temperature
* The hour $\theta$ at which skin temperature is a maximum
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
	θ = zeros(nlon,nlat,12)
	A = zeros(nlon,nlat,12)
	μ = zeros(nlon,nlat,12)

	idat = zeros(nt+1)
	it   = 0:0.1:24; it = it[1:(end-1)]; nit = length(it)
	ts   = longitude2timeshift.(longitude)
	var  = zeros(nit)

	for imo = 1 : 12, ilat = 1 : nlat, ilon = 1 : nlon

		tl = (0:24) .+ ts[ilon]

        idat[1:24] = @view data[ilon,ilat,:,imo]
		idat[end]  = data[ilon,ilat,1,imo]

		itp = interpolate(idat,BSpline(Cubic(Periodic(OnGrid()))))
        stp = scale(itp,tl)
        etp = extrapolate(stp,Periodic())
		var[:] = etp[it]

		μ[ilon,ilat,imo] = mean(@view data[ilon,ilat,:,imo])
		A[ilon,ilat,imo] = (maximum(idat) .- minimum(idat))/2
		θ[ilon,ilat,imo] = argmax(var) / nit * 24

    end

	return μ,A,θ

end

# ╔═╡ 7702f0f7-899f-4333-aba1-7bb55b136e30
ttt   = longitude2timeshift(180)

# ╔═╡ 409521c5-a78e-424d-81ae-21e8df7a581c
md"
### B. Defining GeoRegions
"

# ╔═╡ 37d2e058-22dc-4bd3-b992-b1037d2b8ada
egeo = ERA5Region(GeoRegion("TRP"))

# ╔═╡ 1fadf4ca-5755-11eb-1ece-a99313019785
lsd = getLandSea(egeo,path=datadir("emask"))

# ╔═╡ aa05317e-530b-11eb-2ec1-93aff65659dd
md"
### C. Retrieving ERA5 Skin Temperature Data
"

# ╔═╡ 103f85e8-530c-11eb-047d-a537aa60075d
function retrieveera()

    ds  = NCDataset(datadir("compiled","era5mh-TRPx0.25-sst-compiled.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	sst = nomissing(ds["sst"][:],NaN)
	close(ds)

	μ,A,θ = eradiurnal2model(sst,lon)

    return lon,lat,μ,A,θ

end

# ╔═╡ e8141e20-53af-11eb-1a23-81d34293c5eb
begin
	lon,lat,μ,A,θ = retrieveera()
	md"Modelling diurnal cycle of sea surface temperature"
end

# ╔═╡ bb59b8d6-53b1-11eb-3631-87ef61219c4c
begin
	for imo = 1 : 12
		mstr = uppercase(monthabbr(imo))
		pplt.close(); f,axs = pplt.subplots(nrows=2,axwidth=6,aspect=6)
	
		c = axs[1].contourf(lon,lat,μ[:,:,imo]',levels=295:305,extend="both")
		axs[1].plot(x,y,c="k",lw=0.5)
		axs[1].format(urtitle=L"$\mu$ / K")
		axs[1].colorbar(c,loc="r")
	
		c = axs[2].contourf(lon,lat,A[:,:,imo]'*100,levels=10. .^(-1:0.2:1),extend="both")
		axs[2].plot(x,y,c="k",lw=0.5)
		axs[2].format(urtitle=L"A / 10$^{-2}$ K")
		axs[2].colorbar(c,loc="r",ticks=[0.1,1,10,100])
	
		for ax in axs
			ax.format(
				xlim=(0,360),ylim=(-30,30),xlocator=0:60:360,
				xlabel=L"Longitude / $\degree$",
				ylabel=L"Latitude / $\degree$",
			)
		end
	
		f.savefig(plotsdir("04c-sstspatial_TRP-$mstr.png"),transparent=false,dpi=200)
	end
end

# ╔═╡ c5ed743c-4138-4cb3-906d-7b44d5b0e1cf
mo = 10

# ╔═╡ 0d16e9b1-acc8-4d85-a481-da26729f8bfb
load(plotsdir("04c-sstspatial_TRP-$(uppercase(monthabbr(mo))).png"))

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
	geo = GeoRegion("DTP_IPW")
	md"Defining Regional GeoRegion ..."
end

# ╔═╡ ea7f0956-575b-11eb-3e3f-a1ba3e08b771
begin
	N,S,E,W = geo.N,geo.S,geo.E,geo.W
	ggrd = RegionGrid(geo,lon,lat)
	nlon = length(ggrd.lon)
	nlat = length(ggrd.lat)
	rμ = zeros(nlon,nlat,12)
	rA = zeros(nlon,nlat,12)
	rθ = zeros(nlon,nlat,12)

	for imo = 1 : 12
		iμ = @view μ[:,:,imo]; irμ = @view rμ[:,:,imo]; extractGrid!(irμ,iμ,ggrd)
		iA = @view A[:,:,imo]; irA = @view rA[:,:,imo]; extractGrid!(irA,iA,ggrd)
		iθ = @view θ[:,:,imo]; irθ = @view rθ[:,:,imo]; extractGrid!(irθ,iθ,ggrd)
	end
	md"Extracting information for region ..."
end

# ╔═╡ 5714c13c-575c-11eb-06d4-838b4e8dbcd7
begin
	for imo = 1 : 12
		asp = (E-W+2)/(N-S+2)
		pplt.close()
		if asp > 1.5
			freg,areg = pplt.subplots(nrows=2,axwidth=asp*1.2,aspect=asp)
		else
			freg,areg = pplt.subplots(ncols=2,axwidth=2,aspect=asp)
		end
	
		creg = areg[1].contourf(
			ggrd.lon,ggrd.lat,rμ[:,:,imo]',
			levels=295:305,extend="both"
		)
		areg[1].plot(x,y,c="k",lw=0.5)
		areg[1].format(rtitle=L"$\mu$ / K")
		areg[1].colorbar(creg,loc="r")
	
		creg = areg[2].contourf(
			ggrd.lon,ggrd.lat,rA[:,:,imo]'*100,
			levels=10. .^(-1:0.2:1),extend="both"
		)
		areg[2].plot(x,y,c="k",lw=0.5)
		areg[2].format(rtitle=L"A / 10$^{-2}$ K")
		areg[2].colorbar(creg,loc="r",ticks=[0.1,1,10,100])
	
		for ax in areg
			ax.format(
				xlim=(ggrd.lon[1].-1,ggrd.lon[end].+1),
				xlabel=L"Longitude / $\degree$",
				ylim=(S-1,N+1),ylabel=L"Latitude / $\degree$",
				grid=true
			)
		end
	
		freg.savefig(
			plotsdir("04c-sstspatial_$(geo.regID)-$(uppercase(monthabbr(imo))).png"),
			transparent=false,dpi=200
		)
	end
end

# ╔═╡ 0fdad268-07e0-4a29-88d6-15b66a264722
load(plotsdir("04c-sstspatial_$(geo.regID)-$(uppercase(monthabbr(mo))).png"))

# ╔═╡ Cell order:
# ╟─90fffbc8-524d-11eb-232a-1bada28d5505
# ╟─6f00b8fc-530c-11eb-2242-99d8544f6e14
# ╟─8f30c56c-530c-11eb-2782-33f3c4ed9e89
# ╟─d82366b0-53b1-11eb-26c1-ff1bb6ccb027
# ╟─c1fafcba-530f-11eb-1cc2-67d10a3fa606
# ╟─3565af3c-5311-11eb-34c4-2d228b05b17c
# ╠═fd344560-94d3-11eb-2b79-05c905c5953f
# ╠═a6a688ca-53ab-11eb-2776-b5380ffb26c1
# ╠═7702f0f7-899f-4333-aba1-7bb55b136e30
# ╟─409521c5-a78e-424d-81ae-21e8df7a581c
# ╟─37d2e058-22dc-4bd3-b992-b1037d2b8ada
# ╟─1fadf4ca-5755-11eb-1ece-a99313019785
# ╟─aa05317e-530b-11eb-2ec1-93aff65659dd
# ╠═103f85e8-530c-11eb-047d-a537aa60075d
# ╟─e8141e20-53af-11eb-1a23-81d34293c5eb
# ╟─bb59b8d6-53b1-11eb-3631-87ef61219c4c
# ╠═c5ed743c-4138-4cb3-906d-7b44d5b0e1cf
# ╟─0d16e9b1-acc8-4d85-a481-da26729f8bfb
# ╟─5c0e5bae-554e-11eb-3f83-a364ae0a2485
# ╟─68cfc46c-5755-11eb-1702-373942539652
# ╠═52b39ff8-9426-11eb-2a86-43f7da15f62e
# ╟─ea7f0956-575b-11eb-3e3f-a1ba3e08b771
# ╟─5714c13c-575c-11eb-06d4-838b4e8dbcd7
# ╟─0fdad268-07e0-4a29-88d6-15b66a264722
