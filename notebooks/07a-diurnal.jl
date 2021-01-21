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
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 420e093c-5ba4-11eb-07da-c9a80044c8f1
md"
# 7a. Amplitude of the Diurnal Cycle

In this notebook, we do some basic investigation into the relationship between slab depth and the amplitude of the diurnal cycle when the model is run in RCE state.

Relevant model parameters:
* SST = 301.7 K
* Insolation Peak = 1354.23 W m$^{-2}$
* Momentum Damping $a_m$ = 2
* Momentum Damping $a_m$ exponent = 1
"

# ╔═╡ 57f52568-5bb9-11eb-1e7f-c34b6efe0bac
function retrieveprcp(configlist,beg,tstep,tshift)
	
	nconfig = length(configlist);
	prcp_d  = zeros(tstep+2,nconfig)
	
	for iconfig = 1 : nconfig
		prcp = retrievevar("PREC","DiurnalAmp",configlist[iconfig])
		prcp_d[:,iconfig] = diurnal2D(prcp[(end-beg):end],tstep,tshift);
	end
	
	return prcp_d
	
end

# ╔═╡ 4f40bd7e-5bb9-11eb-34f2-91e1f959c59a
function plotprcpWTG(axs,axsii,td,prcp,configlist,colors)
	
	nconfig = length(configlist)
	
	for iconfig = 1 : nconfig
		config = configlist[iconfig]
		config = replace(config,"slab"=>"")
		config = replace(config,"d"=>".")
		if config != "Inf."
			  config = parse(Float64,config)
		else; config = L"\infty"
		end
		axs[axsii].plot(
			td,prcp[:,iconfig],color=colors[iconfig],
			label="$(config) m",legend="r"
		)
	end
	
end

# ╔═╡ 1ab3ad70-5c2a-11eb-187b-75e524a4581f
begin
	configs = [
		"slab00d1","slab00d2","slab00d3","slab00d5","slab00d7",
		"slab01d0","slab01d4","slab02d0","slab03d2","slab05d0",
		"slab07d1","slab10d0","slab14d1","slab20d0","slabInfd"
	]
	ncon = length(configs)
	colors = pplt.Colors("blues",ncon)
end

# ╔═╡ 6e1170f4-5bb9-11eb-0d38-61befdd2ad88
begin
	ndy = 50
	# z,p,t = retrievedims("DiurnalAmp","slabInfd")
	# td,tstep,tshift,beg = t2d(t,ndy);
	# np = length(p)
md"Loading dimensions from Control Experiment ..."
end

# ╔═╡ 1a31f88e-5c2a-11eb-0c1f-257863a43cf5
# prcp = retrieveprcp(configs,beg,tstep,tshift)

# ╔═╡ 7d401cd6-5c2e-11eb-3842-917545e546ef
begin
	
# 	pplt.close(); f,axs = pplt.subplots(axwidth=4,aspect=2)
	
# 	plotprcpWTG(axs,1,td,prcp_WTG,configs,colors)
# 	axs[1].format(
# 		ylim=(0,6),ylabel=L"Rainfall Rate / mm day$^{-1}$",
# 		xlim=(0,24),xlabel="Hour of Day",
# 	)
	
# 	f.savefig(plotsdir("DiurnalAmp.png"),transparent=false,dpi=200)
# 	PNGFiles.load(plotsdir("DiurnalAmp.png"))
	
end

# ╔═╡ Cell order:
# ╟─420e093c-5ba4-11eb-07da-c9a80044c8f1
# ╟─24abbe3a-5bb7-11eb-160b-1323efad463b
# ╟─24601ef8-5bb7-11eb-1cd8-198dac960d3a
# ╠═57f52568-5bb9-11eb-1e7f-c34b6efe0bac
# ╠═4f40bd7e-5bb9-11eb-34f2-91e1f959c59a
# ╠═1ab3ad70-5c2a-11eb-187b-75e524a4581f
# ╠═6e1170f4-5bb9-11eb-0d38-61befdd2ad88
# ╠═1a31f88e-5c2a-11eb-0c1f-257863a43cf5
# ╠═7d401cd6-5c2e-11eb-3842-917545e546ef
