### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ ba6d524f-b050-4a2c-9772-536c740e3b27
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 00921025-5d26-4356-ae23-7208719d5634
begin
	@quickactivate "TroPrecLS"
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))

	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 87ad146c-175d-11ed-274a-cb219a4294d1
md"
# 06a. Some Diurnal Stats for SAM
"

# ╔═╡ fd21e1a7-99ca-45a0-b596-b19a05801b18
begin
	slbvec = [
		"Slab00d02","Slab00d05","Slab00d10","Slab00d20",
		"Slab00d50","Slab01d00","Slab02d00","Slab05d00",
		"Slab10d00","Slab20d00","Slab50d00",
	]
	dmpvec = [
		"damping01d0","damping01d4","damping02d0","damping03d2",
		"damping05d0","damping07d1","damping10d0",
	]
	nslb = length(slbvec)
	ndmp = length(dmpvec)
end

# ╔═╡ 1c03b953-f2ee-4552-a884-136cf418b5e2
begin
	prcp = zeros(nslb,ndmp)
	for idmp = 1 : ndmp, islb = 1 : nslb
	
		prcp[islb,idmp] = mean(retrievevar("PREC",slbvec[islb],dmpvec[idmp]))
		
	end
end

# ╔═╡ Cell order:
# ╟─87ad146c-175d-11ed-274a-cb219a4294d1
# ╟─ba6d524f-b050-4a2c-9772-536c740e3b27
# ╠═00921025-5d26-4356-ae23-7208719d5634
# ╠═fd21e1a7-99ca-45a0-b596-b19a05801b18
# ╠═1c03b953-f2ee-4552-a884-136cf418b5e2
