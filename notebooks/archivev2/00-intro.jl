### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ d20942a4-51f2-11eb-0977-33f2082de974
md"
# TroPrecLS: Tropical Precipitation over Land and Sea

This project investigates the nature of precipitation dynamics in the tropics over land and sea, and if there are any differences in these dynamics to be found in different geographical areas.

I have created several different notebooks to help better illustrate my thought process during the course of this project.  These notebooks will be mostly used for data analysis and plotting.  For computationally heavy, memory intensive, and lengthy work, scripts will be provided and referenced as need be, so that the user can submit these jobs to a backend/cluster to obtain the data themselves.

The general outline of this project can be split into two parts:
A. Analysis of Observational and Reanalysis datasets
B. Corroborating our results from part (A) with model runs in SAM

### A. Observations and Reanalysis
1. Defining domains within the tropical region to be explored

2. Download and exploratory analysis of satellte rainfall data in the tropical domain and subdomains
    * Analysis of mean, diurnal cycle in domains

3. Download and exploratory analysis of reanalysis ERA5 data in the tropical domain and subdomains
    * Analysis of mean, diurnal cycle in domains
    * Analysis of vertical profiles for pressure-level data

4. Analysis of the seasonal and intraseasonal variability of the satellite and reanalysis data in the tropical domain and subdomains
    * What is the difference between land and sea
    * Is there a difference in this landmass over different land areas

5. Investigating the relationship between column saturation fraction (r) and precipitation rate (P) in reanalysis data
    * What is the difference between land and sea
    * Is there a difference in this landmass over different land areas

6. Investigating the relationship between column saturation fraction (r) and other climatological variables
    * Surface variables such as surface temperature
    * Pressure variables and vertical profiles of vertical wind, air temperature, relative humidity, specific humidity, etc.

### B. Corroborating Observations with SAM Runs
7. Control runs with diurnal cycle
    * Do some control runs with 128x128 km at 2km resolution
    * Get statistics for vertical wind, humidity, etc. and diurnal cycle
    * Run for different damping strengths
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─d20942a4-51f2-11eb-0977-33f2082de974
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
