### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ b5d1e16c-59a0-11ed-3e1a-917178c535a4
begin
	using Pkg; Pkg.activate()
	using DrWatson

md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ abf9aa3c-a2fb-4c4f-ac2c-c2034906b335
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using ERA5Reanalysis
	using NCDatasets
	using Printf
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 88ed76f6-75ab-4a82-8625-a4e4cbcb05ef
begin
	crd = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = crd[:,1]; y = crd[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ 148326c8-59dd-4921-832e-1b628ae83ce2
e5ds = ERA5Daily(start=Date(1979),stop=Date(1979),path=datadir())

# ╔═╡ 13f35777-035c-4095-838a-a9be4ddc67f0
evar = SingleVariable("sst")

# ╔═╡ b0bbcdb3-c71e-445e-816f-7a937d01e662
ereg = ERA5Region("TRP",gres=1)

# ╔═╡ 1328a162-63cf-4d0d-a8a1-8293112b3ba8
lsd = getLandSea(ereg,path=datadir("emask"))

# ╔═╡ a1cbe04b-3d5d-47d0-8cba-1ca3b42fc2c7
begin
	dsr = read(e5ds,evar,ereg,Date(1979))
	ssr = nomissing(dsr[evar.varID][:],NaN)
	close(dsr)
	ds5 = read(e5ds,evar,ereg,Date(1979),smooth=true,smoothlon=5,smoothlat=5)
	ss5 = nomissing(ds5[evar.varID][:],NaN)
	close(ds5)
	ds1 = read(e5ds,evar,ereg,Date(1979),smooth=true,smoothlon=120,smoothlat=10)
	ss1 = nomissing(ds1[evar.varID][:],NaN)
	close(ds1)
end

# ╔═╡ 9c1085cf-bedf-4329-b3ac-2c18bc9d8f1e
begin
	for it = 7 : 25
		pplt.close(); fig,axs = pplt.subplots(nrows=3,aspect=6,axwidth=5)
		axs[1].pcolormesh(
			lsd.lon,lsd.lat,
			dropdims(mean(ssr[:,:,(it-6):(it+6)],dims=3),dims=3)',
			levels=300:0.5:305,extend="both"
		)
		axs[2].pcolormesh(
			lsd.lon,lsd.lat,
			dropdims(mean(ss1[:,:,(it-6):(it+6)],dims=3),dims=3)',
			levels=300:0.5:305,extend="both"
		)
		axs[3].pcolormesh(
			lsd.lon,lsd.lat,
			dropdims(mean(ssr[:,:,(it-6):(it+6)],dims=3),dims=3)' .-
			dropdims(mean(ss1[:,:,(it-6):(it+6)],dims=3),dims=3)',
			levels=vcat(-1,-0.5:0.1:-0.1,-0.05,0.05,0.1:0.1:0.5,1),extend="both"
		)
	
		for ax in axs
			ax.plot(x,y,lw=1,c="k")
			ax.format(ylim=(-20,20),xlim=(30,270))
			# ax.format(ylim=(-30,30),xlim=(30,300))
		end
		
		fig.savefig("test-$it.png",transparent=false,dpi=150)
	end
end

# ╔═╡ 1754ec4c-961b-4fbe-9b0e-1917ca183e51
load("test-7.png")

# ╔═╡ Cell order:
# ╟─b5d1e16c-59a0-11ed-3e1a-917178c535a4
# ╟─abf9aa3c-a2fb-4c4f-ac2c-c2034906b335
# ╟─88ed76f6-75ab-4a82-8625-a4e4cbcb05ef
# ╟─148326c8-59dd-4921-832e-1b628ae83ce2
# ╠═13f35777-035c-4095-838a-a9be4ddc67f0
# ╟─b0bbcdb3-c71e-445e-816f-7a937d01e662
# ╟─1328a162-63cf-4d0d-a8a1-8293112b3ba8
# ╟─a1cbe04b-3d5d-47d0-8cba-1ca3b42fc2c7
# ╟─9c1085cf-bedf-4329-b3ac-2c18bc9d8f1e
# ╠═1754ec4c-961b-4fbe-9b0e-1917ca183e51
