### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 681658b0-5914-11eb-0d65-bbace277d145
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 6dce35fc-5914-11eb-0ce2-0d4e164e1898
begin
	@quickactivate "TroPrecLS"
	using SpecialFunctions
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 6b. Transitioning from RCE to WTG

In this notebook, we investigate and develop a way to implement the WTG forcing gradually in the System of Atmospheric Modelling.  The sudden introduction of WTG large-scale forcing often causes a model to enter a \"shocked\" state that unnaturally forces the model into a different state.  Here, we develop a method that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ c22f840c-5af6-11eb-084d-d30d8891f460
md"
### A. Implementing the WTG Approximation in Models

Blossey et al. [2009] were the first to implement to WTG approximation in SAM.  The strength of the WTG is affected by the momentum-damping coefficient $a_m$ in the below formula:

$$\frac{\partial}{\partial p}
\left( \frac{f^2+a_m^2}{a_m}\frac{\partial\omega'}{\partial p} \right)
\approx \frac{k^2R_d}{p} T_v'$$

In our CRM simulations, we assume that $f$, the coriolis parameter, is 0, within the tropics (and therefore negligible).  Thus, the formula simplifies itself to

$$\frac{\partial}{\partial p}
\left( a_m\frac{\partial\omega'}{\partial p} \right)
\approx \frac{k^2R_d}{p} T_v'$$

If we do a further simplification by assuming that $a_m$ is constant in height (which was noted in Blossey et al. [2009], though they favoured that $a_m \propto p/p_\text{ref}$), this further simplifies the equation down to

$$\frac{\partial^2\omega'}{\partial p^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v'$$

Therefore, we see that for a bigger $a_m$, the induced $\omega'$ (or the WTG-induced vertical adjustment) is smaller.  Thus, bigger $a_m$ would result in less deviation from an RCE solution.  In other words, the magnitude of the WTG forcing is proportional to $a_m^{-1}$.
"

# ╔═╡ 205c5c9a-5914-11eb-09af-4313417b0df2
wtgstrength(am::Real) = 1/am

# ╔═╡ 75ed9ab0-5af7-11eb-1cba-ad51d005bd0d
md"
### B. Implementing a smooth transition
"

# ╔═╡ e7bd22f2-5ae9-11eb-1f82-09333889a5b2
twtg_scale = 0.25

# ╔═╡ 310fd792-5914-11eb-0d81-bd1eb22fab05
begin
	pplt.close(); f,axs = pplt.subplots(nrows=3,aspect=3,axwidth=3,sharey=0);
	
	for am in 2 .^(1:9)
		tmax = 1; twtg_max = twtg_scale * tmax
		t = 0:0.001:tmax; 
		test = erf.((twtg_max/2 .-t)/(twtg_max/5)) * (10-log(2,am))/2 .+ 5 .+ log(2,am)/2
		amc = 2 .^(test)
		wtg = wtgstrength.(amc)

		axs[1].plot(t,test,c="k")
		axs[2].plot(t,amc,c="k")
		axs[3].plot(t,wtg,c="k")
	end
	
	axs[1].plot([1,1]*twtg_scale,[0,10],c="r")
	axs[2].plot([1,1]*twtg_scale,[0,1024],c="r")
	axs[3].plot([1,1]*twtg_scale,[0,1],c="r")
	
	axs[1].format(ylim=(0,10),ylabel=L"$\log_2 (a_m)$")
	axs[2].format(ylim=(0,1024),ylabel=L"a_m")
	axs[3].format(ylim=(0.001,1),ylabel=L"WTG $\propto$ $a_m^{-1}$",yscale="log",
		ylocator=[0.001,0.01,0.1,1],
		xlim=(0,1),xlabel="Nondimensionalized time",
		suptitle="twtg_scale = $(twtg_scale)")
	
	f.savefig(plotsdir("wtgcoeff.png"),transparent=false,dpi=200)
	load(plotsdir("wtgcoeff.png"))
end

# ╔═╡ 2cc46a1a-5b35-11eb-2b84-d536d6f7f7d0
md"
### C. Implementation in SAM
"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
exp = "3P"

# ╔═╡ bdfa9872-5b35-11eb-059e-ad1ac171d295
am_wtg = 1

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	_,_,t = retrievedims("$(exp)WTGamExp$(am_wtg)","damping02d00")
	t = t .- 80
	configvec = [
		"damping02d00","damping04d00","damping08d00","damping16d00",
		"damping32d00","damping64d00","damping128d0","damping256d0",
		"damping512d0"
	]
	ncon = length(configvec)
	lgd = Dict("frame"=>false,"ncols"=>1)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 49ad8588-5b62-11eb-2e9d-a357ab55bc2a
begin
	pplt.close(); fts,ats = pplt.subplots(aspect=1.5,axwidth=3)
	
	for icon in 1 : ncon
		config = configvec[icon]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		prcp = retrievevar("PREC","$(exp)WTGamExp$(am_wtg)",configvec[icon])
		prcp = dropdims(mean(reshape(prcp,24,:),dims=1),dims=1)
		ats[1].plot(
			1:200,prcp,c="Blue$(icon)",
			label=(L"$a_m =$" * " $config"),
			legend="r",legend_kw=lgd
		)
	end
	
	ats[1].format(
		xlim=(100,200),xlabel="Time / Days",
		ylim=(0.1,100),yscale="log",ylabel=L"Precipitation Rate / mm day$^{-1}$",
		suptitle="$(exp)WTGamExp$(am_wtg)"
	)
	fts.savefig(plotsdir(
		"rce2wtg-$(exp)WTGamExp$(am_wtg).png"),
		transparent=false,dpi=200
	)
	load(plotsdir("rce2wtg-$(exp)WTGamExp$(am_wtg).png"))
end

# ╔═╡ 7c9cc168-5b6a-11eb-332a-3bfcfe9453fb
expvec = ["3P"]; nexp = length(expvec)

# ╔═╡ 048e5a08-5b3b-11eb-0b13-efe1a204e66a
begin
	amvec = [0,1];   nam  = length(amvec)
	prcp = zeros(length(configvec))
	
	pplt.close(); fp,ap = pplt.subplots(aspect=1.5,axwidth=3)
	
	for iexp = 1 : nexp, iam in 1 : nam
		for icon in 1 : ncon
			prcp[icon] = mean(retrievevar(
				"PREC","$(expvec[iexp])WTGamExp$(amvec[iam])",configvec[icon]
			)[3601:4800])
		end
		ap[1].plot(
			2 .^(1:9),prcp,
			label=("$(expvec[iexp]), "* L"$a_m =$" * " $(amvec[iam])"),
			legend="r",legend_kw=lgd
		)
		ap[1].scatter(2 .^(1:9),prcp,s=10)
	end
	
	
	ap[1].format(
		xlim=(1,1000),xscale="log",xlabel=L"$a_m$",
		ylim=(0.1,100),yscale="log",ylabel=L"Precipitation Rate / mm day$^{-1}$",
	)
	fp.savefig(plotsdir(
		"rce2wtg-$(exp)WTG-prcp.png"),
		transparent=false,dpi=200
	)
	load(plotsdir("rce2wtg-$(exp)WTG-prcp.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─c22f840c-5af6-11eb-084d-d30d8891f460
# ╟─205c5c9a-5914-11eb-09af-4313417b0df2
# ╟─75ed9ab0-5af7-11eb-1cba-ad51d005bd0d
# ╠═e7bd22f2-5ae9-11eb-1f82-09333889a5b2
# ╟─310fd792-5914-11eb-0d81-bd1eb22fab05
# ╟─2cc46a1a-5b35-11eb-2b84-d536d6f7f7d0
# ╠═d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╠═bdfa9872-5b35-11eb-059e-ad1ac171d295
# ╠═a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─49ad8588-5b62-11eb-2e9d-a357ab55bc2a
# ╠═7c9cc168-5b6a-11eb-332a-3bfcfe9453fb
# ╟─048e5a08-5b3b-11eb-0b13-efe1a204e66a
