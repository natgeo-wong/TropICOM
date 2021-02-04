### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 24abbe3a-5bb7-11eb-160b-1323efad463b
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 24601ef8-5bb7-11eb-1cd8-198dac960d3a
begin
	@quickactivate "TroPrecLS"
	using NCDatasets
	using Statistics
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 420e093c-5ba4-11eb-07da-c9a80044c8f1
md"
# 7a. Surface Energy Imbalance

So far, the sea-surface temperature has been held a constant in all our simulations. However, this means that there is an energy imbalance at the surface.  In WTG, we calculate the surface energy imbalance for an ocean of fixed SST, and then remove this imbalance in our mixed-layer experiments.

In this notebook, we perform these calculations, and what we calculate here will be put into PRM files.
"

# ╔═╡ 57f52568-5bb9-11eb-1e7f-c34b6efe0bac
begin
	ndy = 50
	z,p,t = retrievedims("DiurnalAmp","slabInfd")
	td,tstep,tshift,beg = t2d(t,ndy)
md"Loading dimensional data ..."
end

# ╔═╡ 4f40bd7e-5bb9-11eb-34f2-91e1f959c59a
begin
	lw = retrievevar("LWNS","DiurnalAmp","slabInfd")[(end-beg):end]
	sw = retrievevar("SWNS","DiurnalAmp","slabInfd")[(end-beg):end]
	sh = retrievevar("SHF","DiurnalAmp","slabInfd")[(end-beg):end]
	lh = retrievevar("LHF","DiurnalAmp","slabInfd")[(end-beg):end]
md"Loading surface energy balance data from the fixed SST run ..."
end

# ╔═╡ 1ab3ad70-5c2a-11eb-187b-75e524a4581f
seb = mean(sw .- lw .- sh .- lh)

# ╔═╡ 6e1170f4-5bb9-11eb-0d38-61befdd2ad88
mean(sw), mean(lw), mean(sh), mean(lh)

# ╔═╡ Cell order:
# ╟─420e093c-5ba4-11eb-07da-c9a80044c8f1
# ╟─24abbe3a-5bb7-11eb-160b-1323efad463b
# ╟─24601ef8-5bb7-11eb-1cd8-198dac960d3a
# ╟─57f52568-5bb9-11eb-1e7f-c34b6efe0bac
# ╟─4f40bd7e-5bb9-11eb-34f2-91e1f959c59a
# ╠═1ab3ad70-5c2a-11eb-187b-75e524a4581f
# ╠═6e1170f4-5bb9-11eb-0d38-61befdd2ad88
