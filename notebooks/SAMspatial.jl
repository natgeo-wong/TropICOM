### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ a78f9aae-5443-11eb-2c51-add2147140c2
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ b9a81586-5443-11eb-2f74-0563bee08169
begin
	@quickactivate "TroPrecLS"
	using NCDatasets
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 193437cc-5445-11eb-2d8d-914892b1e129
t = 67

# ╔═╡ 9c926732-5444-11eb-1edc-d5c74fa4363c
begin
	fnc = datadir("Control/3DRCE/OUT_2D/RCE_TroPrecLS-Control_64.2Dbin_1.nc")
	ds  = NCDataset(fnc)
	x = ds["x"][:]/1000; y = ds["y"][:]/1000; tt = ds["time"][t]
	prcp = ds["PW"][:,:,t]
	close(ds)
end

# ╔═╡ 33f5578a-5445-11eb-33ce-278ad2f7c250
begin
	pplt.close(); f,axs = pplt.subplots(axwidth=3)
	c = axs[1].contourf(x,y,prcp',cmap="Blues",levels=0:10)
	axs[1].format(rtitle="Time: $tt")
	f.colorbar(c,loc="r")
	f.savefig("2Dtest.png",transparent=false,dpi=200)
	PNGFiles.load("2Dtest.png")
end

# ╔═╡ Cell order:
# ╠═a78f9aae-5443-11eb-2c51-add2147140c2
# ╠═b9a81586-5443-11eb-2f74-0563bee08169
# ╠═193437cc-5445-11eb-2d8d-914892b1e129
# ╠═9c926732-5444-11eb-1edc-d5c74fa4363c
# ╠═33f5578a-5445-11eb-33ce-278ad2f7c250
