### A Pluto.jl notebook ###
# v0.12.18

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

# ╔═╡ aa05317e-530b-11eb-2ec1-93aff65659dd
md"
### A. Defining Analysis and Plotting functions

We first begin by defining the function to compile raw GPM IMERGv6 precipitation data (time-mean and diurnal cycle).
"

# ╔═╡ 103f85e8-530c-11eb-047d-a537aa60075d
function retrievegpm(
    sroot::AbstractString;
    regID::AbstractString="GLB",
	timeID::Vector
)

    glon,glat = gpmlonlat()
	rlon,rlat,_ = gregiongridvec(regID,lon,lat)
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

    return rlon,rlat,gpm

end

# ╔═╡ c1fafcba-530f-11eb-1cc2-67d10a3fa606
md"We approximate that the diurnal cycle of rainfall to a cosine wave:

$$P(t) = \mu + A\cos\left((t-\theta)\frac{\pi}{12}\right)$$

where $\mu$ is the mean rainfall rate, $A$ is the amplitude of the diurnal cycle of precipitation, and $\theta$ is the hour at which precipitation peaks.  We define this model of precipitation that we will fit the GPM data to using the package `LsqFit`.

Later, as a sanity check, we will compare the $\mu$ calculated using this model against the $\mu_r$ that is found by simply temporal-averaging of the raw data
"

# ╔═╡ a0d8a08c-530f-11eb-3603-9309dcca331e
diurnalcycle(time,params) = params[1] + params[2] * cos.((time .- params[3]) * pi/12)

# ╔═╡ 3565af3c-5311-11eb-34c4-2d228b05b17c
md"
GPS IMERGv6 output time is defined according to UTC.  Therefore, for every grid point, we must also correct to the local time before fitting the data to the model.
"

# ╔═╡ 577a59d8-5311-11eb-2e2e-53bbeff12648
longitude2timeshift(longitude::Real) = longitude / 180 * 12

# ╔═╡ Cell order:
# ╟─90fffbc8-524d-11eb-232a-1bada28d5505
# ╟─6f00b8fc-530c-11eb-2242-99d8544f6e14
# ╠═8f30c56c-530c-11eb-2782-33f3c4ed9e89
# ╟─aa05317e-530b-11eb-2ec1-93aff65659dd
# ╠═103f85e8-530c-11eb-047d-a537aa60075d
# ╟─c1fafcba-530f-11eb-1cc2-67d10a3fa606
# ╠═a0d8a08c-530f-11eb-3603-9309dcca331e
# ╟─3565af3c-5311-11eb-34c4-2d228b05b17c
# ╠═577a59d8-5311-11eb-2e2e-53bbeff12648
