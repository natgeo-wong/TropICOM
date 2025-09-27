### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 68d60a02-3798-11ed-36a6-13483ca696b8
begin
	using Pkg; Pkg.activate()
	using DrWatson
end

# ╔═╡ c8be23ea-3768-4822-99e0-b5bc8f85186d
begin
	@quickactivate "TroPrecLS"
	using NCDatasets
	using PyCall, LaTeXStrings
	using PNGFiles, ImageShow
	
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
end

# ╔═╡ e221ad2f-ee09-43ba-943e-383182efe8d6
begin
	dsw = NCDataset(datadir(
		"Slab00d20","OUT_STAT",
		"WTG_TroPrecLS-Slab00d20-relaxscale01d4.nc"
		# "DGW_TroPrecLS-Slab00d02-damping01d4.nc"
	))
	tw  = dsw["time"][:]
	p   = dsw["p"][:]
	prc = dsw["PREC"][:]
	clw = dsw["CLD"][:]
	close(dsw)
end

# ╔═╡ a13f629e-45ec-4a78-ab04-809b0f8e37ce
begin
	dsd = NCDataset(datadir(
		"Slab00d20","OUT_STAT",
		"DGW_TroPrecLS-Slab00d20-damping02d0.nc"
	))
	td  = dsd["time"][:]
	pd  = dsd["PREC"][:]
	cld = dsd["CLD"][:]
	close(dsd)
end

# ╔═╡ 3a5f67c5-234e-4997-a9b6-5741ef0c6887
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=3,nrows=3)
	axs[1].plot(tw,prc)
	axs[1].plot(tw,pd)
	# axs[2].plot(tw,ins)
	c = axs[2].pcolormesh(tw,p,clw,extend="both")
	c = axs[3].pcolormesh(tw,p,cld,extend="both")
	axs[1].format(xlim=(40,50))
	axs[3].colorbar(c)
	fig.savefig("testwtg.png",transparent=false,dpi=200)
	load("testwtg.png")
end

# ╔═╡ Cell order:
# ╠═68d60a02-3798-11ed-36a6-13483ca696b8
# ╠═c8be23ea-3768-4822-99e0-b5bc8f85186d
# ╠═e221ad2f-ee09-43ba-943e-383182efe8d6
# ╠═a13f629e-45ec-4a78-ab04-809b0f8e37ce
# ╠═3a5f67c5-234e-4997-a9b6-5741ef0c6887
