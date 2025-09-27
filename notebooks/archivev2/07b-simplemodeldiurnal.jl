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
	using DelimitedFiles
	using NCDatasets
	using PlutoUI
	using SimpleIslandModels
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ db521478-61e1-4184-960b-7cab95a48b50
md"
# 7b. The Diurnal Cycle in a Simple Model

Text
"

# ╔═╡ af993da1-6181-4c26-8531-e8537ae629d9
TableOfContents()

# ╔═╡ 12ea7ea3-dba6-403d-8853-560389f95ef8
md"
### A. A simple model of island climatology

In this notebook, we run models of island climatology with a fixed atmosphere (i.e. the atmosphere is under WTG) in representation of tropical climatology.

First, we want to calculate the temperature of the atmosphere when the surface temperature is fixed at 300 K.
"

# ╔═╡ 56ad3673-1925-4623-a8f4-16985191e7a9
sfc_init = CreateSurface(sst=300,isfixed=true,mld=0.1)

# ╔═╡ bf70a374-6e98-44cc-8356-e5d68da00e46
atm_init = CreateAtmosphere(sfc_init,εsw=0.05,εlw=0.95,Ta=200,isfixed=false)

# ╔═╡ 10dd7f4c-6d17-4a28-ada4-9c04305ea32e
model_init = CreateModel(sfc_init,atm_init)

# ╔═╡ 031da261-f7b1-4e8b-b91d-470571576c1d
begin
	_,_,Ta,Fs,_ = solveanalytic(model_init)
	md"We solve for the balanced atmospheric temperature: $Ta K"
end

# ╔═╡ 645016d2-d179-4eea-a5e4-f6ae82ec3dc2
md"Now we solve for the mixed-layer surface temperature using the calculated balanced atmospheric temperature."

# ╔═╡ cc7a7e11-6845-47ac-864a-367624b7f690
sfc = CreateSurface(sst=300,isfixed=false,mld=20)

# ╔═╡ f45828bd-acec-4a3f-867e-d66c73e2613b
atm = CreateAtmosphere(sfc,εsw=0.05,εlw=0.95,Ta=250,isfixed=true)

# ╔═╡ 479a88a0-4207-4ff4-8ad6-eba37f09e282
mld = CreateModel(sfc,atm)

# ╔═╡ 9d0acddc-6166-4517-9815-72657dad2b31
begin
	_,Ts,_,_,_ = solveanalytic(mld)
	md"We solve for the balanced surface temperature: $Ts K"
end

# ╔═╡ 293fa075-4dc1-4d23-8dd5-137ec4f99bbe
vars = run(mld,nsteps=288000,dt=300,nstats=6,Ts0=300,Ta0=Ta)

# ╔═╡ 556b4d51-d416-412e-b0c2-43d2b64a41a6
begin
	pplt.close(); f1,a1 = pplt.subplots(nrows=1,aspect=3,axwidth=3)
	
	a1[1].plot(vars.t/86400,vars.Tₛ)
	a1[1].format(xlim=(0,1000),xlabel="Time / Days",ylabel=L"$T_s$ / K")
	
	f1.savefig("test.png")
	load("test.png")
end

# ╔═╡ 5c862bd8-8823-403e-a788-c0f9fab49b18
md"
We see that the average surface temperature is about 320 K which is higher than the given surface temperature of 300 K. This is not entirely unexpected, because when the surface temperature was fixed, it was actually an energy sink of about $Fs W/m2
"

# ╔═╡ 684ab42d-9da3-4bc9-95c6-70ecb313d3b1
md"
### B. A More Thorough Exploration
"

# ╔═╡ a98281be-6b50-4700-960a-eb8caec6f2c6
begin
	depthlist = [
		0.0001,0.0001*sqrt(2),0.0002,0.0002*sqrt(2.5),0.0005,0.0005*sqrt(2),
		0.001,0.001*sqrt(2),0.002,0.002*sqrt(2.5),0.005,0.005*sqrt(2),
		0.01,0.01*sqrt(2),0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
		0.1,0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
		1.,1*sqrt(2),2.,2*sqrt(2.5),5.,5*sqrt(2),
		10.,10*sqrt(2),20.,20*sqrt(2.5),50.,50*sqrt(2),100.
	]
	
	ndepth = length(depthlist)
	md"Loading vector of mixed-layer depth in m ..."
end

# ╔═╡ bac56247-8499-40ea-a674-155e77219716
begin
	vars_vec = Vector{SimpleIslandModels.Variables}(undef,ndepth)
	for idepth in 1 : ndepth
		sfc_ii = CreateSurface(sst=Ts,isfixed=false,mld=depthlist[idepth])
		mld_ii = CreateModel(sfc_ii,atm)
		vars_vec[idepth] = run(mld_ii,nsteps=2880000,dt=30,nstats=10,Ts0=300,Ta0=Ta)
	end
	md"Running model for surfaces corresponding to different slab depths ..."
end

# ╔═╡ f5d6f889-7786-4f33-86eb-8f011896089b
begin
	Ts_vec = zeros(ndepth)
	for idepth = 1 : ndepth
		Ts_vec[idepth] = mean(vars_vec[idepth].Tₛ[(vars_vec[idepth].t/86400).>900])
	end
	md"Finding the surface temperature mean ..."
end

# ╔═╡ 084ab778-4d01-48d3-93b0-7237af3314f1
lsc = pplt.get_colors("delta_r",17)

# ╔═╡ 925bdb9a-666f-4d21-8a86-508e42cdc41c
begin
	pplt.close()
	f2,a2 = pplt.subplots([[1,2,0],[1,2,5],[3,4,5],[3,4,0]],aspect=1.5,axwidth=2,sharey=0)
	
	a2[1].plot(depthlist[[1,end]],[Ts,Ts],c="k",linestyle="--")
	a2[1].scatter(depthlist[1:2:end],Ts_vec[1:2:end],c="blue7",zorder=5)
	a2[1].text(0.001,Ts-3,L"T_s")
	a2[1].format(xscale="log",xlabel="Slab Depth / m",ylabel=L"$\mu(T_s)$ / K",ylim=(300,325))

	ii = 0
	for idepth = vcat(1,7,13:2:25,31,37)
		global ii += 1
		t = vars_vec[idepth].t[1:288] / 86400 * 24
		Ts_ii = vars_vec[idepth].Tₛ[(vars_vec[idepth].t/86400).>900]
		Ts_di = dropdims(mean(reshape(Ts_ii,288,:),dims=2),dims=2)
		Ts_di = (Ts_di .- minimum(Ts_di)) / (maximum(Ts_di)- minimum(Ts_di))
		a2[5].plot(t,Ts_di,c=lsc[ii+3],label="$(depthlist[idepth]) m",legend="r",legend_kw=Dict("ncol"=>1,"frame"=>false))
	end
	ii = 0
	for idepth = 1 : 2 : 37
		global ii += 1
		t = vars_vec[idepth].t[1:288] / 86400 * 24
		Ts_ii = vars_vec[idepth].Tₛ[(vars_vec[idepth].t/86400).>900]
		Ts_di = dropdims(mean(reshape(Ts_ii,288,:),dims=2),dims=2)
		a2[2].scatter(depthlist[idepth],maximum(Ts_di)-minimum(Ts_di),c="blue7")
		Ts_di = (Ts_di .- minimum(Ts_di)) / (maximum(Ts_di)- minimum(Ts_di))
		a2[3].scatter(depthlist[idepth],argmax(Ts_di)/12,c="blue7")
		a2[4].scatter(depthlist[idepth],argmin(Ts_di)/12,c="blue7")
	end
	a2[5].format(xlim=(0,24),xlabel="Hour of Day",xlocator=0:3:24,ylocator=0:1,yticklabels=[L"min($T_s$)",L"max($T_s$)"])
	a2[2].format(xscale="log",ylabel=L"Amp($T_s$)",yscale="log")
	a2[3].format(xscale="log",ylabel=L"$\theta$(max($T_s$))")
	a2[4].format(xscale="log",ylabel=L"$\theta$(min($T_s$))")

	for ii = 1 : 4
		a2[ii].format(xlim=(1e-4,100))
	end
	
	f2.savefig(plotsdir("07b-simplemodeldiurnal.png"),transparent=false,dpi=400)
	load(plotsdir("07b-simplemodeldiurnal.png"))
end

# ╔═╡ Cell order:
# ╟─db521478-61e1-4184-960b-7cab95a48b50
# ╟─c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
# ╟─8d72de06-476c-4bcf-99b3-a9b469fac93d
# ╟─af993da1-6181-4c26-8531-e8537ae629d9
# ╟─12ea7ea3-dba6-403d-8853-560389f95ef8
# ╟─56ad3673-1925-4623-a8f4-16985191e7a9
# ╟─bf70a374-6e98-44cc-8356-e5d68da00e46
# ╟─10dd7f4c-6d17-4a28-ada4-9c04305ea32e
# ╟─031da261-f7b1-4e8b-b91d-470571576c1d
# ╟─645016d2-d179-4eea-a5e4-f6ae82ec3dc2
# ╟─cc7a7e11-6845-47ac-864a-367624b7f690
# ╟─f45828bd-acec-4a3f-867e-d66c73e2613b
# ╟─479a88a0-4207-4ff4-8ad6-eba37f09e282
# ╟─9d0acddc-6166-4517-9815-72657dad2b31
# ╟─293fa075-4dc1-4d23-8dd5-137ec4f99bbe
# ╟─556b4d51-d416-412e-b0c2-43d2b64a41a6
# ╟─5c862bd8-8823-403e-a788-c0f9fab49b18
# ╟─684ab42d-9da3-4bc9-95c6-70ecb313d3b1
# ╟─a98281be-6b50-4700-960a-eb8caec6f2c6
# ╟─bac56247-8499-40ea-a674-155e77219716
# ╟─f5d6f889-7786-4f33-86eb-8f011896089b
# ╟─084ab778-4d01-48d3-93b0-7237af3314f1
# ╠═925bdb9a-666f-4d21-8a86-508e42cdc41c
