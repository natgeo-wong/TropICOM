### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 6f00b8fc-530c-11eb-2242-99d8544f6e14
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 8f30c56c-530c-11eb-2782-33f3c4ed9e89
begin
	@quickactivate "TroPrecLS"
	using ClimateSatellite
	using Dates
	using DelimitedFiles
	using GeoRegions
	using Interpolations
	using LsqFit
	using NCDatasets
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 90fffbc8-524d-11eb-232a-1bada28d5505
md"
# 2. Investigation of Precipitation Data

We first investigate GPM IMERGv6 precipitation data from 2001-2019 over the domains listed in the notebook `01-domains`.  The data was retrieved from the website of NASA's Precipitation Measurement Mission (PMM), and is found in 0.1º x 0.1º horizontal resolution every half-hour.

The data was downloaded using the `ClimateSatellite.jl` package, which downloads the data in monthly batches in the date range specified.  More information can be found in the package documentation, and instructions on how to use the package to download the data can be found there as well.
"

# ╔═╡ c1fafcba-530f-11eb-1cc2-67d10a3fa606
md"
### A. Modelling the Diurnal Cycle of Precipitation

We approximate that the diurnal cycle of rainfall to a cosine wave:

$$P(t) = \mu + A\cos\left((t-\theta)\frac{\pi}{12}\right)$$

where $\mu$ is the mean rainfall rate, $A$ is the amplitude of the diurnal cycle of precipitation, and $\theta$ is the hour at which precipitation peaks.  We define this model of precipitation that we will fit the GPM data to using the package `LsqFit`.

Later, as a sanity check, we will compare the $\mu$ calculated using this model against the $\mu_r$ that is found by simply temporal-averaging of the raw data
"

# ╔═╡ a0d8a08c-530f-11eb-3603-9309dcca331e
diurnalcycle(time,params) = params[1] .+ params[2] * cos.((time .- params[3]) * pi/12)

# ╔═╡ 3565af3c-5311-11eb-34c4-2d228b05b17c
md"
GPS IMERGv6 output time is defined according to UTC.  Therefore, for every grid point, we must also correct to the local time before fitting the data to the model.
"

# ╔═╡ 577a59d8-5311-11eb-2e2e-53bbeff12648
longitude2timeshift(longitude::Real) = longitude / 180 * 12

# ╔═╡ a6a688ca-53ab-11eb-2776-b5380ffb26c1
function gpm2model(data,longitude)
	
	nlon,nlat,nt = size(data)
	θ = zeros(nlon,nlat)
	A = zeros(nlon,nlat)
	μ = zeros(nlon,nlat)
	
	p0 = [0,0.5,0]; t = ((1:nt).-0.5)/2
	ts = longitude2timeshift.(longitude)
	
	for ilat = 1 : nlat, ilon = 1 : nlon

        tl = t .+ ts[ilon]
        fit = curve_fit(diurnalcycle,tl,(@view data[ilon,ilat,:]),p0)
		
        if fit.param[2] < 0
			  A[ilon,ilat] = fit.param[2] * -1
              θ[ilon,ilat] = mod(fit.param[3]+12,24)
		else; A[ilon,ilat] = fit.param[2]
			  θ[ilon,ilat] = mod(fit.param[3],24)
        end
		
		μ[ilon,ilat] = fit.param[1]

    end
	
	return μ,A,θ
	
end

# ╔═╡ aa05317e-530b-11eb-2ec1-93aff65659dd
md"
### B. Retrieving GPM Precipitation Data
"

# ╔═╡ bb90be66-554c-11eb-05de-a5db553ad4b1
md"
We fit the model to GPM IMERGv6 data from 2001 to 2018 for each individual spatial point (i.e. at 0.1º resolution).  Our results are plotted below:
"

# ╔═╡ a116023a-53ad-11eb-25e0-d5dfa8338a1b
function savediurnalphase(regID,rlon,rlat,μt,μ,A,θ)
	
	fnc = datadir("gpmimerg-prcp_rate-$regID.nc")
	if isfile(fnc); rm(fnc) end
	ds = NCDataset(fnc,"c",attrib = Dict("Conventions"=>"CF-1.6"));
	
	ds.dim["longitude"] = length(rlon)
    ds.dim["latitude"]  = length(rlat)
	
	nclongitude = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclatitude = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))
	
	ncμ = defVar(ds,"modelmean",Float32,("longitude","latitude",),attrib = Dict(
        "units" => "kg m**-2 s**-1"
    ))
	
	ncA = defVar(ds,"modelamplitude",Float32,("longitude","latitude",),attrib = Dict(
        "units" => "kg m**-2 s**-1"
    ))
	
	ncθ = defVar(ds,"modelphase",Float32,("longitude","latitude",),attrib = Dict(
        "units" => "hr"
    ))
	
	ncμt = defVar(ds,"rawmean",Float32,("longitude","latitude",),attrib = Dict(
        "units" => "kg m**-2 s**-1"
    ))
	
	nclongitude[:] = rlon
	nclatitude[:]  = rlat
	ncμ[:]  = μ
	ncA[:]  = A
	ncθ[:]  = θ
	ncμt[:] = μt
	
	close(ds)
	
end

# ╔═╡ 103f85e8-530c-11eb-047d-a537aa60075d
function retrievegpm(
    sroot::AbstractString;
    regID::AbstractString="GLB",
	timeID::Vector
)

    glon,glat = gpmlonlat()
	rlon,rlat,_ = gregiongridvec(regID,glon,glat)
    nlon = length(rlon)
	nlat = length(rlat)

    gpm  = zeros(nlon,nlat,48)

    nt = timeID[2] - timeID[1] + 1
    for yr = timeID[1] : timeID[2]

        rds,rvar = clisatanaread(
            "gpmimerg","prcp_rate","domain_yearly_mean_hourly",
            Date(yr),regID,path=sroot
        );
        gpm += rvar[:]*3600; close(rds);

    end
	
	gpm .= gpm / nt
	μ,A,θ = gpm2model(gpm,rlon)
	μt = dropdims(mean(gpm,dims=3),dims=3)
	
	savediurnalphase(regID,rlon,rlat,μt,μ,A,θ)

    return rlon,rlat,μt,μ,A,θ

end

# ╔═╡ 49d13e5c-53af-11eb-29ca-c994a7acd377
sroot = "/n/kuangdss01/lab/"; regID = "TRP"

# ╔═╡ e8141e20-53af-11eb-1a23-81d34293c5eb
if isfile(datadir("gpmimerg-prcp_rate-$regID.nc"))
	ds = NCDataset(datadir("gpmimerg-prcp_rate-$regID.nc"))
	lon = ds["longitude"][:]
	lat = ds["latitude"][:]
	μ  = ds["modelmean"][:]
	A  = ds["modelamplitude"][:]
	θ  = ds["modelphase"][:]
	μt = ds["rawmean"][:]
	close(ds)
else lon,lat,μt,μ,A,θ = retrievegpm(sroot,regID=regID,timeID=[2001,2018])
end

# ╔═╡ d82366b0-53b1-11eb-26c1-ff1bb6ccb027
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
md"Loading coastlines ..."
end

# ╔═╡ bb59b8d6-53b1-11eb-3631-87ef61219c4c
begin
	pplt.close(); f,axs = pplt.subplots(nrows=3,axwidth=5,aspect=6)
	
	c = axs[1].contourf(
		lon,lat,μ',
		cmap="Blues",cmap_kw=Dict("left"=>0.1),
		levels=vcat((1:5)/10),extend="max"
	)
	axs[1].plot(x,y,c="k",lw=0.5)
	axs[1].format(rtitle="Yearly Mean Precipitation / mm")
	axs[1].colorbar(c,loc="r")
	
	c = axs[2].contourf(
		lon,lat,A',
		cmap="Blues",cmap_kw=Dict("left"=>0.1),
		levels=vcat(0.05,0.1,0.2,0.5,1),extend="max"
	)
	axs[2].plot(x,y,c="k",lw=0.5)
	axs[2].format(rtitle=L"A / mm hr$^{-1}$")
	axs[2].colorbar(c,loc="r")
	
	c = axs[3].pcolormesh(lon,lat,θ',cmap="romaO",levels=0:0.5:24)
	axs[3].plot(x,y,c="k",lw=0.5)
	axs[3].format(rtitle=L"$\theta$ / Hour of Day")
	axs[3].colorbar(c,loc="r")
	
	for ax in axs
		ax.format(xlim=(0,360),ylim=(-30,30),xlocator=0:60:360)
	end
	
	f.savefig(plotsdir("precipmodel_TRP.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("precipmodel_TRP.png"))
end

# ╔═╡ 5c0e5bae-554e-11eb-3f83-a364ae0a2485
md"
The precipitation band is concentrated along a narrow band in the tropical oceans that represents the ITCZ and SPCZ, except around the Indo-Pacific warmpool, where we see that regions of high precipitation extend farther in latitude about the equator.

We see that the presence of land also serves to extend the distribution of precipitation over latitude.  Furthermore, the amplitude $A$ of the diurnal cycle is enhanced over land as compared to over the ocean for similar values of mean precipitation rate.

Lastly, we also see that the peak of precipitation rainfall is drastically different over land compared to the ocean.  Over land, the precipitation rate generally peaks from 1800 hours to midnight, while over the tropical ocean the rainfall peaks in the early morning.  We also see that over the tropical ocean the peak is relatively uniform, whereas the phase $\theta$ becomes much noisier in the extratropical regions, perhaps because $A$ is so much smaller that it becomes harder to distinguish $A$ from the background noise.
"

# ╔═╡ a96bfb80-5554-11eb-1fab-21f167010eea
md"We also compared $\mu_t$ and $\mu$, or the raw mean precipitation rate against the model-derived mean precipitation rate, and found that there was no difference between the two."

# ╔═╡ caa8956a-5554-11eb-250e-931b02fae447
begin
	pplt.close(); f1,axs1 = pplt.subplots(axwidth=5,aspect=6)
	
	c1 = axs1[1].contourf(
		lon,lat,(μ-μt)',cmap="RdBu_r",
		levels=vcat(-5:-1,1:5)/100,extend="both"
	)
	axs1[1].plot(x,y,c="k",lw=0.5)
	axs1[1].format(rtitle=L"$\mu - \mu_t$ / mm hr$^{-1}$")
	axs1[1].colorbar(c1,loc="b")
	
	for ax in axs1
		ax.format(xlim=(0,360),ylim=(-30,30),xlocator=0:60:360)
	end
	
	f1.savefig(plotsdir("muvmut.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("muvmut.png"))
end

# ╔═╡ 68cfc46c-5755-11eb-1702-373942539652
md"
### C. Regional Analysis

Let's get a quick snapshot of the results by region.
"

# ╔═╡ ea7f0956-575b-11eb-3e3f-a1ba3e08b771
begin
	rlon,rlat,rinfo = gregiongridvec("SEA",lon,lat)
	if maximum(rlon) > 360; rlon .= rlon .- 360 end
	rμ = regionextractgrid(μ,rinfo)
	rA = regionextractgrid(A,rinfo)
	rθ = regionextractgrid(θ,rinfo)
end

# ╔═╡ 5714c13c-575c-11eb-06d4-838b4e8dbcd7
begin
	asp = (maximum(rlon)-minimum(rlon))/(maximum(rlat)-minimum(rlat))
	pplt.close(); freg,areg = pplt.subplots(nrows=3,axwidth=4,aspect=asp)
	
	creg = areg[1].contourf(
		rlon,rlat,rμ',
		cmap="Blues",cmap_kw=Dict("left"=>0.1),
		levels=vcat((1:7)/10),extend="max"
	)
	areg[1].plot(x,y,c="k",lw=0.5)
	areg[1].format(rtitle="Yearly Mean Precipitation / mm")
	areg[1].colorbar(creg,loc="r")
	
	creg = areg[2].contourf(
		rlon,rlat,rA',
		cmap="Blues",cmap_kw=Dict("left"=>0.1),
		levels=vcat(0.05,0.1,0.2,0.3,0.4,0.5,1),extend="max"
	)
	areg[2].plot(x,y,c="k",lw=0.5)
	areg[2].format(rtitle=L"A / mm hr$^{-1}$")
	areg[2].colorbar(creg,loc="r")
	
	creg = areg[3].pcolormesh(rlon,rlat,rθ',cmap="romaO",levels=0:0.5:24)
	areg[3].plot(x,y,c="k",lw=0.5)
	areg[3].format(rtitle=L"$\theta$ / Hour of Day")
	areg[3].colorbar(creg,loc="r")
	
	for ax in areg
		ax.format(
			xlim=(minimum(rlon),maximum(rlon)),
			ylim=(minimum(rlat),maximum(rlat))
		)
	end
	
	freg.savefig(plotsdir("precipmodel_reg.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("precipmodel_reg.png"))
end

# ╔═╡ c4792bf2-5552-11eb-3b52-997f59fd42f3
md"
### D. Binning of Precipitation Data

We now wish to bin the modelled precipitation data in order to determine the relationship between precipitation rate and the diurnal cycle over land and sea.  Is it different?
"

# ╔═╡ 1fadf4ca-5755-11eb-1ece-a99313019785


# ╔═╡ Cell order:
# ╟─90fffbc8-524d-11eb-232a-1bada28d5505
# ╟─6f00b8fc-530c-11eb-2242-99d8544f6e14
# ╟─8f30c56c-530c-11eb-2782-33f3c4ed9e89
# ╟─c1fafcba-530f-11eb-1cc2-67d10a3fa606
# ╠═a0d8a08c-530f-11eb-3603-9309dcca331e
# ╟─3565af3c-5311-11eb-34c4-2d228b05b17c
# ╠═577a59d8-5311-11eb-2e2e-53bbeff12648
# ╠═a6a688ca-53ab-11eb-2776-b5380ffb26c1
# ╟─aa05317e-530b-11eb-2ec1-93aff65659dd
# ╟─bb90be66-554c-11eb-05de-a5db553ad4b1
# ╟─103f85e8-530c-11eb-047d-a537aa60075d
# ╟─a116023a-53ad-11eb-25e0-d5dfa8338a1b
# ╠═49d13e5c-53af-11eb-29ca-c994a7acd377
# ╠═e8141e20-53af-11eb-1a23-81d34293c5eb
# ╠═d82366b0-53b1-11eb-26c1-ff1bb6ccb027
# ╟─bb59b8d6-53b1-11eb-3631-87ef61219c4c
# ╟─5c0e5bae-554e-11eb-3f83-a364ae0a2485
# ╟─a96bfb80-5554-11eb-1fab-21f167010eea
# ╟─caa8956a-5554-11eb-250e-931b02fae447
# ╟─68cfc46c-5755-11eb-1702-373942539652
# ╠═ea7f0956-575b-11eb-3e3f-a1ba3e08b771
# ╟─5714c13c-575c-11eb-06d4-838b4e8dbcd7
# ╟─c4792bf2-5552-11eb-3b52-997f59fd42f3
# ╠═1fadf4ca-5755-11eb-1ece-a99313019785
