### A Pluto.jl notebook ###
# v0.19.16

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
		"IslandSize","Damping01d0","OUT_STAT",
		"DGW_TroPrecLS-IslandSize-010km.nc"
	))
	t    = ds["time"][:]
	p    = ds["p"][:]
	prcp = ds["PREC"][:] / 24
	pwv  = ds["PW"][:]
	sst  = ds["SST"][:]
	ins  = ds["SOLIN"][:]
	tdif = ds["TABS"][:] .- ds["TABSOBS"][:]
	wwtg = ds["WWTG"][:]
	close(ds)
end

# ╔═╡ a136825d-94b9-495f-aae4-fa84e519cdef
begin
	pplt.close(); fig,axs = pplt.subplots(nrows=2,aspect=2,axwidth=2,sharey=0)
	
	axs[1].plot(t,prcp)
	axs[2].plot(t,sst)

	axs[1].format(ylim=(0,100),ylabel=L"Rain / mm hr$^{-1}$")
	axs[2].format(ylim=(280,320),ylabel="SST / K")

	for ax in axs
		ax.format(xlim=(0,15),xlabel="Days")
	end
	
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
