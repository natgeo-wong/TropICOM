### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ e56a0d84-1e8e-11ed-27ce-c3cdf389563d
begin
	using Pkg; Pkg.activate()
	using DrWatson
end

# ╔═╡ ad2fcb70-fc14-4b02-84a1-9244ef8e740e
begin
	@quickactivate "TroPrecLS"
	using NCDatasets
	using StatsBase
	
	using PyCall, LaTeXStrings
	using PNGFiles, ImageShow
	
	pplt = pyimport("proplot")
end

# ╔═╡ 5f4a3ed0-56af-429d-8143-41178e18f971
begin
	ds  = NCDataset(datadir("MockWalker","OUT_2D","RCE_MockWalker_64.2Dbin_1.nc"))
	t   = ds["time"][:]
	x   = ds["x"][:]
	olr = ds["LWNT"][:]
	close(ds)

	x = dropdims(mean(reshape(x,32,:),dims=1),dims=1)
	t = dropdims(mean(reshape(t,24,:),dims=1),dims=1)
	olr = dropdims(mean(reshape(olr,32,:,24,300),dims=(1,3)),dims=(1,3))
end

# ╔═╡ b471fc0f-96d5-4e70-bd16-ef91aabc5456
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=1/2)
	
	axs[1].pcolormesh(x,t,olr')
	
	fig.savefig("testmockwalker.png",transparent=false,dpi=150)
	load("testmockwalker.png")
end

# ╔═╡ 06d91eeb-8581-4ee8-baa2-dcc46e04f83d
begin
	d3D = NCDataset(datadir("MockWalker","OUT_3D","RCE_MockWalker_64.bin2D_3.nc"))
	x3D = d3D["x"][:]
	p3D = d3D["p"][:]
	qrd = d3D["QN"][:]
	close(d3D)
end

# ╔═╡ 764f57f8-ff2f-4b62-89f4-3196bc9c9936
begin
	pplt.close(); f3D,a3D = pplt.subplots(aspect=6,axwidth=6)
	
	a3D[1].pcolormesh(x3D,p3D,log10.(qrd[:,:,5]'))
	a3D[1].format(yscale="log")
	
	f3D.savefig("testmockwalker3D.png",transparent=false,dpi=150)
	load("testmockwalker3D.png")
end

# ╔═╡ Cell order:
# ╟─e56a0d84-1e8e-11ed-27ce-c3cdf389563d
# ╠═ad2fcb70-fc14-4b02-84a1-9244ef8e740e
# ╟─5f4a3ed0-56af-429d-8143-41178e18f971
# ╟─b471fc0f-96d5-4e70-bd16-ef91aabc5456
# ╠═06d91eeb-8581-4ee8-baa2-dcc46e04f83d
# ╠═764f57f8-ff2f-4b62-89f4-3196bc9c9936
