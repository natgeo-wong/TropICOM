### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ cadb2d1c-5055-11eb-1f40-edf74c826843
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ f05f7958-5055-11eb-3b7d-4723da6d65a1
begin
	@quickactivate "TroPrecLS"
	using Statistics
	using NCDatasets
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("samsnd.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 2b1aec94-5051-11eb-04e6-379d1a4cfa3f
md"
# X. Sounding Profiles from RCE Simulations

Based on the Control RCE run `sst301d70`, we construct overall mean, and diurnal cycle, sounding files based on the last 100 days of the simulation.  These functions `createsndmean` and `createsnddiurnal` can be found and explored in `src/samsnd.jl`.
"

# ╔═╡ 4fae4bc0-5059-11eb-124b-e121399c73e5
md"
### A. Creating the Sounding Profile

We create the sounding profiles based on the `Control/sst301d70` simulation.
"

# ╔═╡ bd914aaa-5056-11eb-0070-438d7c3ff9c2
begin
	p3S,snd3S = createsndmean(
		"initsnd3S",exp="Control",config="3SRCE",
		psfc=1009.32,extract=true
	)
	p3P,snd3P = createsndmean(
		"initsnd3P",exp="Control",config="3PRCE",
		psfc=1009.32,extract=true
	)
	# p3P,snd3P = createsndmean(
	# 	"initsnd3P",exp="Control",config="3PRCE",
	# 	psfc=1009.32,ndays=50,extract=true
	# )
	length(p3P)
end

# ╔═╡ 4f4f3016-506f-11eb-1206-391a46bd6579
md"
### B. Sanity check for the Sounding Profile

Below, we plot the sounding profiles against the pressure levels that were extracted from the `OUT_STAT` files, both the mean profile (of temperature $T$ and potential temeprature $\theta$) and the diurnal cycle (of $T$):
"

# ╔═╡ 17986926-505b-11eb-172c-b93d26450f0b
begin
	pplt.close()
	f,axs = pplt.subplots(ncols=2,aspect=1/3,axwidth=0.75,sharex=0,wspace=0.1)
	
	tem3S = snd3S[:,3] .* (p3S/p3S[1]).^0.287
	tem3P = snd3P[:,3] .* (p3P/p3P[1]).^0.287
	
	lgd = Dict("ncols"=>1,"frame"=>false)
	axs[1].plot(tem3S-tem3P,p3S,lw=1)
	# axs[1].plot(snd3P[:,3] .* (p3P/p3P[1]).^0.287,p3S,lw=1)
	axs[1].format(
		#xlim=(180,320),xlabel="T / K",
		ylim=(1010,50),ylabel="Pressure / hPa"
	)

	axs[2].plot(snd3S[:,3],p3S,lw=1,label="3S",legend="r",legend_kw=lgd)
	axs[2].plot(snd3P[:,3],p3P,lw=1,label="3P",legend="r")
	axs[2].format(xlabel=L"$\theta$ / K",xlim=(270,500))
	
	f.savefig(plotsdir("soundmean.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("soundmean.png"))
end

# ╔═╡ Cell order:
# ╟─2b1aec94-5051-11eb-04e6-379d1a4cfa3f
# ╠═cadb2d1c-5055-11eb-1f40-edf74c826843
# ╠═f05f7958-5055-11eb-3b7d-4723da6d65a1
# ╟─4fae4bc0-5059-11eb-124b-e121399c73e5
# ╠═bd914aaa-5056-11eb-0070-438d7c3ff9c2
# ╟─4f4f3016-506f-11eb-1206-391a46bd6579
# ╠═17986926-505b-11eb-172c-b93d26450f0b
