### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 417ee688-5ade-11eb-2e95-91a301119e88
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 46faa412-5ade-11eb-3c37-23a7e59037a0
begin
	@quickactivate "TroPrecLS"
	using Statistics
	using PlutoUI
	using Printf
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	include(srcdir("samsnd.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 9dd4cd7e-5adb-11eb-2735-a7a4a2bb23b1
md"
# 6a. Finding an Equilibrium RCE State

This notebook will be used to help find an equilibrium RCE state that roughly returns the initial sounding profile used to initialize the model.

SAM reads in \"snd\" (sounding profiles of humidity and potential temperature) files to initialize the model upon spinup.  We run models in RCE mode over 400 day periods and obtain the profile of absolute temperature for the last 50 days.  We then compare it to the `TABSOBS` that was fed into the model, and find the temperature difference.  If the difference is too large, then we extract the vertical sounding profile, override the old \"snd\" file and rerun the model.
"

# ╔═╡ 86f0a60a-5ae5-11eb-0b1a-935e8703b842
md"
### A. Extracting Temperature Data

First, we extract the temperature and observed temperature (`TABS` and `TABSOBS` respectively) data from the statistical output.  We have modified the SAM code such that `TABSOBS` is invariant in height and time.  Thus, we plot the difference between the two variables for the last 50 days of the RCE model run, after a 350 day spinup.
"

# ╔═╡ ac8b9d4c-5ade-11eb-06f4-33bff063bbde
config = "3P"

# ╔═╡ 32acc944-624b-11eb-14df-556c55bd89e4
ndys = 200

# ╔═╡ ab0909d0-5ade-11eb-2622-0fee6450004b
begin
	z,p,t = retrievedims("Control","$(config)RCE")
	v2D = retrievevar("PREC","Control","$(config)RCE")
	tem = retrievevar("TABS","Control","$(config)RCE")
	tob = retrievevar("TABSOBS","Control","$(config)RCE")
md"Loading data from $(config)RCE run ..."
end

# ╔═╡ e82cd648-5ade-11eb-2739-ad93c0b6d3af
begin
	tdiff = dropdims(mean(tem[:,(end-ndys+1):end],dims=2),dims=2).-tob[:,1]
	tdts  = tem .- tob[:,1]
md"Calculating difference between `TABS` and `TABSOBS` ..."
end

# ╔═╡ 0c8cc7ec-5ae3-11eb-1c52-0fecbb575c43
begin
	trms = sqrt(mean(tdiff.^2))
	trms = @sprintf("%.3f",trms)
	
md"Assuming that the tropopause is at 100 hPa, the root-mean-square of the temperature difference between the model temperature `TABS` and the observed temperature `TABSOBS` is $(trms) K.  The profile of the temperature difference is shown below:"
end

# ╔═╡ 978d0442-5b33-11eb-39e5-cdea9c6efa4c
begin
	
	pplt.close(); fts,ats = pplt.subplots(
		[1,2,2,2],aspect=1/3,axwidth=0.5,
		sharex=0
	)
	
	ats[1].plot(tdiff,p,lw=0.5)
	ats[1].scatter(tdiff,p,s=3)
	ats[1].format(
		xlim=(-0.15,0.15),ylim=(1010,15),yscale="log",
		xlabel=L"T - T$_{OBS}$ / K",ylabel="Pressure / hPa",
	)
	
	c = ats[2].contourf(
		t.-80,p,tdts,
		cmap="RdBu_r",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=vcat(-10,-7.07,-5,-3.16,-2,-1.41,-1,-0.5,0.5,1,1.41,2,3.16,5,7.07,10)/10
	)
	ats[2].format(
		xlim=(800,1000),
		ylim=(1010,15),yscale="log",
		xlabel="time / days",ylabel="Pressure / hPa",
		urtitle=L"T$_{RMS}$" * " = $(trms) K"
	)
	
	fts.colorbar(c,loc="r",width=0.2)
	fts.savefig("rcetdts-$(config).png",transparent=false,dpi=200)
	load("rcetdts-$(config).png")
	
end

# ╔═╡ 39868b46-5ae5-11eb-0bf5-274642855e12
md"
### B. Creating Sounding Profile

Regardless, we create the sounding profile based on the last 50 days of RCE data.  We then have to decide for ourselves whether we want to rerun the RCE model based on this new sounding data, or to move on to using this RCE profile for WTG simulations.
"

# ╔═╡ 4fd6e272-5b32-11eb-2f60-bbd4f8b9fc12
md"Create SND file? $(@bind dosnd PlutoUI.Slider(0:1))"

# ╔═╡ 094999c8-5ae5-11eb-1526-f38c604184cb
if isone(dosnd)
	createsndmean(
		"initsnd$(config)",
		exp="Control",config="$(config)RCE",
		psfc=1009.32,ndays=50
	)
md"Creating the SND file initsnd$(config) ..."
else
md"We have decided not to create the SND file initsnd$(config) yet ..."
end

# ╔═╡ Cell order:
# ╟─9dd4cd7e-5adb-11eb-2735-a7a4a2bb23b1
# ╟─417ee688-5ade-11eb-2e95-91a301119e88
# ╟─46faa412-5ade-11eb-3c37-23a7e59037a0
# ╟─86f0a60a-5ae5-11eb-0b1a-935e8703b842
# ╠═ac8b9d4c-5ade-11eb-06f4-33bff063bbde
# ╠═32acc944-624b-11eb-14df-556c55bd89e4
# ╠═ab0909d0-5ade-11eb-2622-0fee6450004b
# ╠═e82cd648-5ade-11eb-2739-ad93c0b6d3af
# ╠═0c8cc7ec-5ae3-11eb-1c52-0fecbb575c43
# ╠═978d0442-5b33-11eb-39e5-cdea9c6efa4c
# ╟─39868b46-5ae5-11eb-0bf5-274642855e12
# ╟─4fd6e272-5b32-11eb-2f60-bbd4f8b9fc12
# ╟─094999c8-5ae5-11eb-1526-f38c604184cb
