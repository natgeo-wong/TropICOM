### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 247edb6d-65a6-4fb7-beab-46604b78cfe9
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 5ebccf8b-f701-436c-b654-9acb6856473d
begin
	@quickactivate "TroPrecLS"
	using Printf
	using Statistics
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 440fd8da-8f4a-11eb-108e-993e7f81183e
md"
# 9d. Binning CSF values

As an extension to notebook 09c, we do a quick binning of the CSF values for the different slabs.

Relevant model parameters:
* SST = 301.7 K
* Insolation Peak = 1354.23 W m$^{-2}$
* Momentum Damping $a_m$ = 2
* Momentum Damping $a_m$ exponent = 0
"

# ╔═╡ d8248390-1416-4322-bebe-9b913e4cdc1c
md"
### A. Loading Data for WTG Control Experiments
"

# ╔═╡ 41ed6c57-c7fc-4bbe-bcb7-08414d9f84b3
function retrieve2D(variable,configlist)
	
	nconfig = length(configlist);
	d2D  = zeros(4800,5,nconfig)
	
	for ic = 1 : nconfig
		for ii = 1 : 5
			fnc = outstatname("Control",configlist[ic],false,true,ii)
			d2D[:,ii,ic] .= retrievevar(variable,fnc)[9601:end]
		end
	end
	
	d2D = reshape(d2D,4800*5,nconfig)
	
	return d2D
	
end

# ╔═╡ eed5fd10-e1c8-448e-80c7-44eea38eca01
function retrievecsf(exp,configlist)
	
	nconfig = length(configlist);
	var  = zeros(4800,5,nconfig)
	
	for ic = 1 : nconfig
		for ii = 1 : 5
			fnc = outstatname(exp,configlist[ic],false,true,ii)
			tcw = retrievevar("PW",fnc)[9601:end] / 1000
			ta  = retrievevar("TABS",fnc)[:,9601:end]
			qv  = retrievevar("QV",fnc)[:,9601:end] / 1000
			pp  = retrievevar("p",fnc)
			rh  = calcrh(qv,ta,pp)
			_,swp = calccsf(rh,qv,pp)
			var[:,ii,ic] += tcw ./ swp
		end
	end
	
	var = reshape(var,4800*5,nconfig)
	
	return var
	
end

# ╔═╡ 4a1ad4f6-5a83-467f-823e-b3f5d92591b7
begin
	configs = [
		"slab00d05","slab00d07","slab00d10","slab00d14","slab00d20",
		"slab00d32","slab00d50","slab00d71","slab01d00","slab01d41",
		"slab02d00","slab03d16","slab05d00","slab07d07","slab10d00",
		"slab14d14","slab20d00","slab31d62","slab50d00",
	]
	ncon = length(configs)
	colors = pplt.Colors("blues",ncon+4)
	lndocn = pplt.Colors("Delta_r",ncon+4)
	md"Defining experimental configurations for WTG experiments ..."
end

# ╔═╡ 375f5fe4-cd08-4361-8976-577d17fd0ee3
begin
	prcp = retrieve2D("PREC",configs) / 24
	csf  = retrievecsf("Control",configs) * 100
md"Loading results from the WTG Slab-depth experiments ..."
end

# ╔═╡ 8d552f5e-42e5-4042-bc3c-a7600a1c8ed6
begin
	cmean = dropdims(mean(csf,dims=1),dims=1)
end

# ╔═╡ 2d4de3a6-d575-49ce-992d-e661a5cc377c
size(csf)

# ╔═╡ 9f580ce2-4dc2-4d05-ad26-0382a5e17b54
md"
### B. Binning Precipitation Rate against Column Saturation Fraction
"

# ╔═╡ 9892a59f-6a61-48a6-aeb2-8aecd66416a4
function csfvprcpbin(prcp,csf,cvec,csep)
	pvec = zeros(length(cvec)); jj = 0;
    pfrq = zeros(Int64,length(cvec))
    for cii in cvec
        pii = @view prcp[ (csf.>(cii-csep)) .& (csf.<=(cii+csep)) ]
        jj = jj + 1; pvec[jj] = mean(pii); pfrq[jj] = length(pii)
    end

    return pvec,pfrq
end

# ╔═╡ 355247fb-e98f-4a20-b7b9-3ce6f227253a
begin
	cbin = collect(70:0.2:100); cstep = (cbin[2]-cbin[1])/2; nbins = length(cbin)
	pfrq = zeros(nbins,ncon)
	pvec = zeros(nbins,ncon)
	cfrq = zeros(nbins,ncon)
	
	for icon = 1 : ncon
		pvec[:,icon],pfrq[:,icon] = csfvprcpbin(prcp[:,icon],csf[:,icon],cbin,cstep)
	end
	
	for icon = 1 : ncon, ibin = 5 : (nbins-4)
		cfrq[ibin,icon] = pfrq[ibin-4,icon] * 0.0625 + pfrq[ibin-3,icon] * 0.125  +
		                  pfrq[ibin-2,icon] * 0.25   + pfrq[ibin-1,icon] * 0.5    +
		                  pfrq[ibin+1,icon] * 0.5    + pfrq[ibin+2,icon] * 0.25   +
		                  pfrq[ibin+3,icon] * 0.125  + pfrq[ibin+4,icon] * 0.0625 +
		                  pfrq[ibin,icon]
	end
	cfrq = cfrq / 3
	md"Finding CSF frequencies ..."
end

# ╔═╡ 55babf48-524d-42cf-b814-73be9361b727
begin
	pplt.close(); ff,af = pplt.subplots(aspect=3,axwidth=6)
	
	for icon = 1 : 19
		label = configs[icon]
		label = replace(label,"slab"=>"")
		label = replace(label,"d"=>".")
		label = parse(Float64,label)
		label = @sprintf("%05.2f",label)
		label = "$(label) m"
		af[1].plot(
			cbin,cfrq[:,icon]./sum(pfrq[:,icon])*nbins,c=lndocn[icon+2],
			label=label,legend="r",legend_kw=Dict("ncol"=>2,"frame"=>false)
		)
		af[1].plot(
			ones(2)* cmean[icon],[0,7.5],c=lndocn[icon+2],
		)
	end
	af[1].format(
		xlim=(70,100),xlabel="Column Relative Humidity",
		ylim=(0,7.5),ylabel="Probability Density"
	)
	
	ff.savefig(plotsdir("samcsffreq.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("samcsffreq.png"))
end

# ╔═╡ 7b41f5c9-f5db-4167-8ea3-25fa9225cce2
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=2,axwidth=4)
	
	for icon = 1 : 19
		label = configs[icon]
		label = replace(label,"slab"=>"")
		label = replace(label,"d"=>".")
		label = parse(Float64,label)
		label = @sprintf("%05.2f",label)
		label = "$(label) m"
		a2[1].plot(
			cbin,pvec[:,icon]./pfrq[:,icon],c=lndocn[icon+2],
			label=label,legend="r",legend_kw=Dict("ncol"=>2,"frame"=>false)
		)
	end
	a2[1].format(
		xlim=(75,100),xlabel="Column Relative Humidity",
		ylim=10. .^(-5,1),ylabel=L"Precipitation Rate / mm hr$^{-1}$",yscale="log"
	)
	
	f2.savefig(plotsdir("samcsfvprcpinit.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("samcsfvprcpinit.png"))
end

# ╔═╡ Cell order:
# ╟─440fd8da-8f4a-11eb-108e-993e7f81183e
# ╟─247edb6d-65a6-4fb7-beab-46604b78cfe9
# ╟─5ebccf8b-f701-436c-b654-9acb6856473d
# ╟─d8248390-1416-4322-bebe-9b913e4cdc1c
# ╠═41ed6c57-c7fc-4bbe-bcb7-08414d9f84b3
# ╠═eed5fd10-e1c8-448e-80c7-44eea38eca01
# ╟─4a1ad4f6-5a83-467f-823e-b3f5d92591b7
# ╠═375f5fe4-cd08-4361-8976-577d17fd0ee3
# ╠═8d552f5e-42e5-4042-bc3c-a7600a1c8ed6
# ╠═2d4de3a6-d575-49ce-992d-e661a5cc377c
# ╟─9f580ce2-4dc2-4d05-ad26-0382a5e17b54
# ╟─9892a59f-6a61-48a6-aeb2-8aecd66416a4
# ╟─355247fb-e98f-4a20-b7b9-3ce6f227253a
# ╟─55babf48-524d-42cf-b814-73be9361b727
# ╠═7b41f5c9-f5db-4167-8ea3-25fa9225cce2
