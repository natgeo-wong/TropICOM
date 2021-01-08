### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ cadb2d1c-5055-11eb-1f40-edf74c826843
begin
	using DrWatson
	using Pkg
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ f05f7958-5055-11eb-3b7d-4723da6d65a1
begin
	@quickactivate "TroPrecLS"
	Pkg.instantiate()
	using Statistics
	
	using ImageShow, FileIO, ImageMagick
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
	p,msnd = createsndmean(
		"initsndmean",exp="Control",config="sst301d70",
		psfc=1009.32,extract=true
	)
	dt,p,dsnd = createsnddiurnal(
		"initsnddiurnal",exp="Control",config="sst301d70",
		psfc=1009.32,extract=true
	)
	p,dt
end

# ╔═╡ 4f4f3016-506f-11eb-1206-391a46bd6579
md"
### B. Sanity check for the Sounding Profile

Below, we plot the sounding profiles against the pressure levels that were extracted from the `OUT_STAT` files, both the mean profile (of temperature $T$ and potential temeprature $\theta$) and the diurnal cycle (of $T$):
"

# ╔═╡ 17986926-505b-11eb-172c-b93d26450f0b
begin
	pplt.close()
	f,axs = pplt.subplots([1,2,2,2,2],aspect=1/3,axwidth=1,sharex=0,wspace=0.1)
	
	axs[1].plot(msnd[:,3] .* (p/p[1]).^0.287,p,c="k",lw=1,label="T",legend="lr")
	axs[1].plot(
		msnd[:,3],p,c="k",lw=1,linestyle="-.",
		label=L"\theta",legend="lr",legend_kw=Dict("ncols"=>1,"frame"=>false)
	)
	axs[1].format(
		xlim=(180,500),xlabel="Temperature / K",
		ylim=(1010,50),ylabel="Pressure / hPa",
		title="Vertical Profile"
	)
	
	temp = dsnd[:,3,:] .* (p/p[1]).^0.287
	c = axs[2].contourf(
		dt*24,p,temp .- mean(temp,dims=2),
		cmap="RdBu_r",levels=vcat(-5:-1,1:5)/20,
		extend="both"
	)
	axs[2].format(
		xlim=(0,24),xlabel="Hour of Day",xlocator=0:4:24,
		ylim=(1010,50),ylabel="Pressure / hPa",
		title="Diurnal Cycle"
	)
	
	f.colorbar(c,loc="r")
	f.savefig("soundmean.png",transparent=false,dpi=200)
	load("soundmean.png")
end

# ╔═╡ cfcba23e-505f-11eb-22a3-5d50a50e7002


# ╔═╡ c7801c6c-5060-11eb-02bd-d319b41c1a29


# ╔═╡ Cell order:
# ╟─2b1aec94-5051-11eb-04e6-379d1a4cfa3f
# ╠═cadb2d1c-5055-11eb-1f40-edf74c826843
# ╠═f05f7958-5055-11eb-3b7d-4723da6d65a1
# ╟─4fae4bc0-5059-11eb-124b-e121399c73e5
# ╠═bd914aaa-5056-11eb-0070-438d7c3ff9c2
# ╟─4f4f3016-506f-11eb-1206-391a46bd6579
# ╠═17986926-505b-11eb-172c-b93d26450f0b
# ╠═cfcba23e-505f-11eb-22a3-5d50a50e7002
# ╠═c7801c6c-5060-11eb-02bd-d319b41c1a29
