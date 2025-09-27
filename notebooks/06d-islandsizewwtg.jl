### A Pluto.jl notebook ###
# v0.19.25

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
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ db521478-61e1-4184-960b-7cab95a48b50
md"
# 6d. Mean Island $w_{wtg}$

Text
"

# ╔═╡ af993da1-6181-4c26-8531-e8537ae629d9
TableOfContents()

# ╔═╡ a99d31a1-4ef8-4cc9-894c-342fd5424e36
begin
	sizelist = [
		2,2*sqrt(2.5),5,5*sqrt(2),
	    10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100,
	    100*sqrt(2),200,200*sqrt(2.5),500,500*sqrt(2),1000,
	    1000*sqrt(2),2000,2000*sqrt(2.5),5000,
	]
	nsize = length(sizelist)
	md"Loading vector of Island sizes in km ..."
end

# ╔═╡ 37a716a7-2425-4f9d-960f-a0be3744a223
begin
	depthlist = [
		0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
		0.1,0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
		1.,1*sqrt(2),2.,2*sqrt(2.5),5.,5*sqrt(2),
		10.,10*sqrt(2),20.,20*sqrt(2.5),50.
	]
	ndepth = length(depthlist)
	md"Loading vector of mixed-layer depth in m ..."
end

# ╔═╡ e8499402-439d-4de2-8be2-1a43b2b506d5
begin
	blues_WTG = pplt.get_colors("roma",(nsize))
	md"Loading Colors ..."
end

# ╔═╡ 18245629-9a24-4d97-a4c3-80877b19929b
function retrieve3D(varID,sizelist,depthlist)

	nsize  = length(sizelist)
	ndepth = length(depthlist)
	z    = zeros(64,nsize,ndepth) * NaN
	p    = zeros(64,nsize,ndepth) * NaN
	d3D  = zeros(64,nsize,ndepth) * NaN

	for id in 1 : ndepth
		depthstr = @sprintf("%05.2f",depthlist[id])
        depthstr = replace(depthstr,"."=>"d")
		for is in 1 : nsize
			sizestr = @sprintf("%04d",sizelist[is])
			config  = "size$(sizestr)km-depth$(depthstr)m"

			fnc = outstatname("DGW","IslandSize3D",config)
			if isfile(fnc)
				ds  = NCDataset(fnc)
				if length(ds["time"][:]) >= 4800
					var = retrievevar_fnc(varID,fnc)
					d3D[:,is,id] .= mean(var,dims=2)
					p[:,is,id] .= ds["p"][:]
					z[:,is,id] .= ds["z"][:]
				end
				close(ds)
			end
		end
	end

	return z,p,d3D

end

# ╔═╡ b2675e23-18ae-4ca1-abb6-e39a0a5be617
function findtropopause(tair,z)

	nlvl = size(tair,1)
	dims = size(tair)[2:end]
	tair = reshape(tair,nlvl,:); nexp = size(tair,2)
	z    = reshape(z,nlvl,:)
	
	tgrad = (tair[3:end,:] .- tair[1:(end-2),:]) ./ 
	        (z[1:(end-2),:] .- z[3:end,:])
	tropp = zeros(Int,nexp)

	for iexp = 1 : nexp

		tgradii = @view tgrad[:,iexp]
		itrop = 20
		while tgradii[itrop] > 2e-3
			itrop += 1
		end
		tropp[iexp] = itrop
		
	end

	tropp = reshape(tropp,dims)

	return tropp

end

# ╔═╡ 64ab3846-2de2-4ca3-a48d-018c6bd5aac7
function baroclinicweights(wwtg,z,ktrop)
	
	nlvl  = size(wwtg,1)
	dims  = size(wwtg)[2:end]
	wwtg  = reshape(wwtg,nlvl,:); nexp = size(wwtg,2)
	z     = reshape(z,nlvl,:)
	ktrop = reshape(ktrop,:)
	nbaro = 16
	weights = zeros(nbaro,nexp)

	for iexp = 1 : nexp

		iw = @view wwtg[:,iexp]
		iz = @view z[:,iexp]
		iwgt = @view weights[:,iexp]
		nk = ktrop[nexp]
		iibl = findfirst(iz.>1000)
		if isnothing(iibl); iibl = 20 end
		ztrop = z[nk,nexp]

		for isin = 1 : nbaro

			iwgt[isin] = iw[iibl] * sin(pi*iz[iibl]*isin/ztrop) * iz[iibl]
			
		end

		for ik = (iibl+1) : nk, isin = 1 : nbaro
			iwgt[isin] += (iw[ik] * sin(pi*iz[ik]*isin/ztrop) + 
			               iw[ik-1] * sin(pi*iz[ik-1]*isin/ztrop)) *
			              (iz[ik] - iz[ik-1])
		end

		iwgt .= iwgt / sum(abs.(iwgt))
		
	end

	weights = reshape(weights,tuple(vcat(nbaro,collect(dims))...))

	return weights
	
end

# ╔═╡ 1455f008-06c8-4f79-a852-ca7d4a324fe8
md"
### A. $w_{wtg}$ Statistics
"

# ╔═╡ 5c04acab-2807-4526-858f-7601b368130b
begin
	z,p,wwtg = retrieve3D("WWTG",sizelist,depthlist)
	_,_,tair = retrieve3D("TABS",sizelist,depthlist)
	pp = zeros(64)
	for ii = 1 : 64
		ip = @view p[ii,:,:]; ip = ip[.!isnan.(ip)]
		pp[ii] = mean(ip)
	end
	md"Binning mean large-scale vertical velocity resulting"
end

# ╔═╡ ef8788d0-a5b0-4e0b-8858-89b1b28c6714
begin
	tropp = findtropopause(tair,z)
	md"Finding index at which tropopause occurs ..."
end

# ╔═╡ 18674823-e7c9-4d47-97ce-7a56827b3c8e
begin
	baromodes = baroclinicweights(wwtg,z,tropp)
	md"Decompose the vertical velocity profiles into baroclinic modes ..."
end

# ╔═╡ 8cacf996-b900-41e5-a184-2a79d785d1fa
begin
	pplt.close()
	fig,axs = pplt.subplots(nrows=2,ncols=3,aspect=1.5,axwidth=2,hspace=1,wspace=1)

	lvls = (0.01:0.01:0.1)*2
	lvls = vcat(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50) / 10
	lvls = vcat(-reverse(lvls),lvls)
	c1 = axs[1].contourf(depthlist,pp,wwtg[:,end-2,:],extend="both",levels=lvls)
	axs[2].contourf(depthlist,pp,wwtg[:,end-4,:],extend="both",levels=lvls)
	axs[3].contourf(depthlist,pp,wwtg[:,end-6,:],extend="both",levels=lvls)
	axs[4].contourf(depthlist,pp,wwtg[:,end-8,:],extend="both",levels=lvls)
	axs[5].contourf(depthlist,pp,wwtg[:,end-10,:],extend="both",levels=lvls)
	axs[6].contourf(depthlist,pp,wwtg[:,end-12,:],extend="both",levels=lvls)
	
	axs[1].plot(depthlist,pp[tropp[end,:]])
	axs[2].plot(depthlist,pp[tropp[end-1,:]])
	axs[3].plot(depthlist,pp[tropp[end-2,:]])
	axs[4].plot(depthlist,pp[tropp[end-3,:]])
	axs[5].plot(depthlist,pp[tropp[end-4,:]])
	axs[6].plot(depthlist,pp[tropp[end-5,:]])

	axs[1].format(ultitle=L"(a) $a_m\lambda^2$ = 20$\times$10$^5$ km$^2$ day$^{-1}$")
	axs[2].format(ultitle=L"(b) $a_m\lambda^2$ = 10$\times$10$^5$ km$^2$ day$^{-1}$")
	axs[3].format(ultitle=L"(c) $a_m\lambda^2$ = 5$\times$10$^5$ km$^2$ day$^{-1}$")
	axs[4].format(ultitle=L"(d) $a_m\lambda^2$ = 2$\times$10$^5$ km$^2$ day$^{-1}$")
	axs[5].format(ultitle=L"(e) $a_m\lambda^2$ = 1$\times$10$^5$ km$^2$ day$^{-1}$")
	axs[6].format(ultitle=L"(f) $a_m\lambda^2$ = 0.5$\times$10$^5$ km$^2$ day$^{-1}$")

	for ax in axs
		ax.format(
			xscale="log",xlim=(0.02,50),xlabel="Mixed Layer Depth / m",
			yscale="log",ylim=(1000,25),ylabel="Pressure / hPa"
		)
	end
	fig.colorbar(c1,length=0.8,label=L"$w_{wtg}$ / m s$^{-1}$")
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ a2c25a4a-24aa-4ac1-9dd2-29f0bb4633a0
begin
	pplt.close()
	fbaro,abaro = pplt.subplots(aspect=0.9,axwidth=1.5,ncols=3)

	for ax in abaro
		ax.format(
			xscale="log",xlim=(0.05,50),ylabel="Mixed Layer Depth / m",
			yscale="log",ylim=(0.02,50),xlabel=L"$a_m\lambda^2$ / 10$^5$ km$^2$ day$^{-1}$",
			# ylocator=[0.1,0.2,0.5,1,2,5,10,20],
			# xlocator=[10,20,50,100,200,500,1000]/100
		)
	end

	barolvls = 0.25:0.05:0.75
	barolvls = vcat(-reverse(barolvls),0,barolvls)
	cbaro = abaro[1].pcolormesh(
		sizelist/100,depthlist,baromodes[1,:,:]',cmap="RdBu_r",
		levels=barolvls,extend="both"
	)

	abaro[2].pcolormesh(
		sizelist/100,depthlist,
		(baromodes[2,:,:]./baromodes[1,:,:])',cmap="RdBu_r",
		levels=barolvls,extend="both"
	)

	abaro[3].pcolormesh(
		sizelist/100,depthlist,
		dropdims(sum(abs.(baromodes[1:2,:,:]),dims=1),dims=1)',cmap="RdBu_r",
		levels=barolvls,extend="both"
	)

	abaro[1].format(ultitle=L"(a) $c_1$ / $\sum |c_i|$")
	abaro[2].format(ultitle=L"(b) $c_2$ / $c_1$")
	abaro[3].format(ultitle=L"(c) (|$c_1$| + |$c_2$|) / $\sum c_i$")

	fbaro.colorbar(cbaro,locator=vcat(-0.75:0.1:-0.25,0.25:0.1:0.75))
	fbaro.savefig(plotsdir("06d-islandsize-wwtgbaro.png"),transparent=false,dpi=400)
	load(plotsdir("06d-islandsize-wwtgbaro.png"))
end

# ╔═╡ 3d14c84b-6114-4d6d-87f8-c116e1dcad56
md"So we see that the first two baroclinic modes of $w_{wtg}$ make up 70% of the entire vertical profile.  So, we can roughly just approximate things using the first and second baroclinic modes"

# ╔═╡ ac6ad433-01f9-48f2-9896-0b9a32571d13
md"
### B. Difference with Large-Scale Simulations ...
"

# ╔═╡ 9fbbe7b6-578b-4639-bee9-0a47dfd1cccc
begin
	ds  = NCDataset(datadir("SAMvertprofile-control.nc"))
	bin = ds["bin"][:]
	zls = ds["z"][:]
	lvl = ds["level"][:]
	ta  = ds["t"][:]
	w   = ds["w"][:]
	close(ds)
	md"Loading the binned 3D profiles"
end

# ╔═╡ f8830c27-a878-4425-8970-6034d11f24ed
tropp_ls = findtropopause(ta',zls)

# ╔═╡ b9bc61e8-3e4e-48f6-873e-95e54df0f67b


# ╔═╡ 4b1ea115-b382-4ba2-94cf-bfe54987fc48
begin
	bm1 = zeros(101)
	bm2 = zeros(101)
	for ib = 15 : 70
		bm1[ib] = -0.007 * sqrt(ib-15) / sqrt(55)
		bm2[ib] = 0.005 * sqrt(ib-15) / sqrt(55)
	end
	for ib = 70 : 85
		bm1[ib] = -0.007 + 0.007 * (ib-70)^2 / (15)^2
		bm2[ib] = 0.005 -  0.045 * (ib-70)^2 / (15)^2
	end
	for ib = 70 : 85
		bm2[ib] = 0.005 -  0.055 * (ib-70)^2 / (15)^2
	end
	for ib = 85 : 101
		bm1[ib] = 7.5 * ((ib-70)/(30))^10
		bm2[ib] = -0.05 + 5 * ((ib-85)/(15))^10
	end
	wa = zeros(64,101)
		for ibin = 1 : 99, iz = 1 : 64
		if iz < 43
			wa[iz,ibin] = bm1[ibin] * sin(pi*zls[iz]/zls[43]) .+
			        	  bm2[ibin] * sin(2pi*zls[iz]/zls[43])
		end
	end
end

# ╔═╡ c8202948-2cce-4715-8b00-7f9b6bb98994
begin
	pplt.close(); f1,a1 = pplt.subplots(aspect=3,axwidth=5,nrows=2)

	c = a1[1].pcolormesh(bin,lvl,w',levels=lvls*1,extend="both")
	a1[1].plot(bin,lvl[tropp_ls])
	a1[1].format(yscale="log")

	a1[2].pcolormesh(bin,lvl,wa,levels=lvls*1,extend="both")
	a1[2].format(yscale="log")
	p1 = a1[2].panel("b")
	p1.plot(0:100,bm1)
	p1.plot(0:100,bm2)
	p1.format(yscale="symlog",yscale_kw=Dict("linthresh"=>0.01))

	f1.colorbar(c)
	f1.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─db521478-61e1-4184-960b-7cab95a48b50
# ╟─c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
# ╟─8d72de06-476c-4bcf-99b3-a9b469fac93d
# ╟─af993da1-6181-4c26-8531-e8537ae629d9
# ╟─a99d31a1-4ef8-4cc9-894c-342fd5424e36
# ╟─37a716a7-2425-4f9d-960f-a0be3744a223
# ╟─e8499402-439d-4de2-8be2-1a43b2b506d5
# ╠═18245629-9a24-4d97-a4c3-80877b19929b
# ╠═b2675e23-18ae-4ca1-abb6-e39a0a5be617
# ╠═64ab3846-2de2-4ca3-a48d-018c6bd5aac7
# ╟─1455f008-06c8-4f79-a852-ca7d4a324fe8
# ╟─5c04acab-2807-4526-858f-7601b368130b
# ╟─ef8788d0-a5b0-4e0b-8858-89b1b28c6714
# ╟─18674823-e7c9-4d47-97ce-7a56827b3c8e
# ╟─8cacf996-b900-41e5-a184-2a79d785d1fa
# ╟─a2c25a4a-24aa-4ac1-9dd2-29f0bb4633a0
# ╟─3d14c84b-6114-4d6d-87f8-c116e1dcad56
# ╟─ac6ad433-01f9-48f2-9896-0b9a32571d13
# ╟─9fbbe7b6-578b-4639-bee9-0a47dfd1cccc
# ╟─f8830c27-a878-4425-8970-6034d11f24ed
# ╟─b9bc61e8-3e4e-48f6-873e-95e54df0f67b
# ╟─4b1ea115-b382-4ba2-94cf-bfe54987fc48
# ╟─c8202948-2cce-4715-8b00-7f9b6bb98994
