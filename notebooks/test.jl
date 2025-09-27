### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 8d72de06-476c-4bcf-99b3-a9b469fac93d
begin
	@quickactivate "TroPrecLS"
	using ERA5Reanalysis
	using NCDatasets
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 01907cf6-48cf-4405-a2c5-d6e29c081878
begin
	ds = NCDataset(datadir(
		"Aggregation","OUT_2D",
		"RCE_TroPrecLS-Aggregation-spinup_64_0000118080.2Dbin_1.nc"
	))
	olr = ds["LWNT"][:,:,1]
	close(ds)
end

# ╔═╡ a136825d-94b9-495f-aae4-fa84e519cdef
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=8,axwidth=8,sharey=0)
	
	axs[1].pcolormesh(olr[1:512,:]')
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ 05f081b2-fe94-4b2e-ba49-7c8a1bb1b49e
begin
	pplt.close(); f2,a2 = pplt.subplots(axwidth=1.5)
	
	a2[1].plot(pwv,prcp)
	a2[1].format(ylim=(0.001,100),yscale="log",xlim=(10,50))
	
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
# ╟─8d72de06-476c-4bcf-99b3-a9b469fac93d
# ╠═01907cf6-48cf-4405-a2c5-d6e29c081878
# ╠═a136825d-94b9-495f-aae4-fa84e519cdef
# ╠═05f081b2-fe94-4b2e-ba49-7c8a1bb1b49e
