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
# 6b. Analysis of Island Sizes

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

# ╔═╡ 1455f008-06c8-4f79-a852-ca7d4a324fe8
md"
### A. Precipitation and Surface Temperature Statistics
"

# ╔═╡ bb4b054f-87ca-4076-97be-207097ceb79e
begin
	prcpmat = zeros(nsize,ndepth) * NaN
	pmaxmat = zeros(nsize,ndepth) * NaN
	pminmat = zeros(nsize,ndepth) * NaN
	pampmat = zeros(nsize,ndepth) * NaN
	ptmxmat = zeros(nsize,ndepth) * NaN
	ptmnmat = zeros(nsize,ndepth) * NaN
	tsfcmat = zeros(nsize,ndepth) * NaN
	tmaxmat = zeros(nsize,ndepth) * NaN
	tminmat = zeros(nsize,ndepth) * NaN
	tampmat = zeros(nsize,ndepth) * NaN
	ttmxmat = zeros(nsize,ndepth) * NaN
	ttmnmat = zeros(nsize,ndepth) * NaN
	md"Preallocation of climatological arrays"
end

# ╔═╡ 1e9f2a5f-3a84-43e3-9ff8-a8adce7c8077
for idepth in 1 : ndepth
	depthstr = @sprintf("%05.2f",depthlist[idepth])
	depthstr = replace(depthstr,"."=>"d")
	for isize in 1 : nsize
		sizestr = @sprintf("%04d",sizelist[isize])
		config = "size$(sizestr)km-depth$(depthstr)m"

		fnc = outstatname("DGW","IslandSize3D",config)
		if isfile(fnc)
			ds  = NCDataset(fnc)
			if length(ds["time"][:]) >= 4800
				time = ds["time"][1:48] .- floor(ds["time"][1])
				prcp = ds["PREC"][:] / 24
				ptmx = dropdims(mean(reshape(prcp,48,:),dims=2),dims=2)
				prcpmat[isize,idepth] = mean(ptmx)
				pmaxmat[isize,idepth] = maximum(ptmx)
				pminmat[isize,idepth] = minimum(ptmx)
				pampmat[isize,idepth] = maximum(ptmx) - minimum(ptmx)
				ptmxmat[isize,idepth] = time[argmax(ptmx)] * 24
				if !iszero(minimum(ptmx))
					ptmnmat[isize,idepth] = time[argmin(ptmx)] * 24
				else
					ptmnmat[isize,idepth] = NaN
				end
				tsfc = ds["SST"][:]
				ttmx = dropdims(mean(reshape(tsfc,48,:),dims=2),dims=2)
				tsfcmat[isize,idepth] = mean(ttmx) - 300
				tmaxmat[isize,idepth] = maximum(ttmx) - 300
				tminmat[isize,idepth] = minimum(ttmx) - 300
				tampmat[isize,idepth] = maximum(ttmx) - minimum(ttmx)
				ttmxmat[isize,idepth] = time[argmax(ttmx)] * 24
				ttmnmat[isize,idepth] = time[argmin(ttmx)] * 24
			end
			close(ds)
		end

	end
end

# ╔═╡ 4ab80acb-afe4-45ad-8698-1131a77c5c79
begin
	pplt.close()
	fprcp,aprcp = pplt.subplots(axwidth=1.5,nrows=2,ncols=3,wspace=[1.5,3])

	for ax in aprcp
		ax.format(
			xscale="log",xlim=(0.05,50),ylabel="Mixed Layer Depth / m",
			yscale="log",ylim=(0.02,50),xlabel=L"$a_m\lambda^2$ / 10$^5$ km$^2$ day$^{-1}$",
			# ylocator=[0.1,0.2,0.5,1,2,5,10,20],
			# xlocator=[10,20,50,100,200,500,1000]/100
		)
	end
	
	cprcp_1 = aprcp[1].pcolormesh(
		sizelist/100,depthlist,prcpmat',cmap="blues",
		levels=10. .^(-1:0.1:2),extend="both"
	)
	ctxt_1 = aprcp[1].contour(
		sizelist/100,depthlist,prcpmat',levels=10. .^(-1:2),
		c="k",linestyle="--",lw=0.5
	); aprcp[1].clabel(ctxt_1,inline=true,fontsize=10)
	cprcp_2 = aprcp[2].pcolormesh(
		sizelist/100,depthlist,pmaxmat',cmap="blues",
		levels=10. .^(-1:0.2:2),extend="both"
	)
	ctxt_2 = aprcp[2].contour(
		sizelist/100,depthlist,pmaxmat',levels=10. .^(-1:2),
		c="k",linestyle="--",lw=0.5
	); aprcp[2].clabel(ctxt_2,inline=true,fontsize=10)
	cprcp_3 = aprcp[3].pcolormesh(
		sizelist/100,depthlist,ptmxmat',levels=0:24,cmap="RdBu",
	)
	cprcp_4 = aprcp[4].pcolormesh(
		sizelist/100,depthlist,pampmat',cmap="blues",
		levels=10. .^(-1:0.2:2),extend="both"
	)
	ctxt_4 = aprcp[4].contour(
		sizelist/100,depthlist,pampmat',levels=10. .^(-1:2),
		c="k",linestyle="--",lw=0.5
	); aprcp[4].clabel(ctxt_4,inline=true,fontsize=10)
	cprcp_5 = aprcp[5].pcolormesh(
		sizelist/100,depthlist,pminmat',cmap="blues",
		levels=10. .^(-1:0.2:2),extend="both"
	)
	ctxt_5 = aprcp[5].contour(
		sizelist/100,depthlist,pminmat',levels=vcat(0,10. .^(-1:2)),
		c="k",linestyle="--",lw=0.5
	); aprcp[5].clabel(ctxt_5,inline=true,fontsize=10)
	cprcp_6 = aprcp[6].pcolormesh(
		sizelist/100,depthlist,ptmnmat',levels=0:24,cmap="RdBu"
	)

	aprcp[1].format(ltitle=L"(a) $\mu(P)$",suptitle="Rainfall Climatology")
	aprcp[2].format(ltitle=L"(c) max($P$)")
	aprcp[3].format(ltitle=L"(e) $\theta$(max($P$))")
	aprcp[4].format(ltitle=L"(b) max($P$) - min($P$)")
	aprcp[5].format(ltitle=L"(d) min($P$)")
	aprcp[6].format(ltitle=L"(f) $\theta$(min($P$))")

	fprcp.colorbar(cprcp_1,length=0.8,loc="l",label=L"Rainfall Rate $P$ / mm hr$^{-1}$",locator=10. .^(-2:2))
	fprcp.colorbar(cprcp_6,length=0.8,label="Hour of Day",locator=0:3:24,minorlocator=[])
	
	fprcp.savefig(plotsdir("06b-islandsize-precip.png"),transparent=false,dpi=400)
	load(plotsdir("06b-islandsize-precip.png"))
end

# ╔═╡ 99e126ae-c746-4e83-8fc5-c09a5dd7980d
begin
	pplt.close()
	ftsfc,atsfc = pplt.subplots(axwidth=1.5,nrows=2,ncols=3,wspace=[1.5,3])

	for ax in atsfc
		ax.format(
			xscale="log",xlim=(0.05,20),ylabel="Mixed Layer Depth / m",
			yscale="log",ylim=(0.02,50),xlabel=L"$a_m\lambda^2$ / 10$^5$ km$^2$ day$^{-1}$",
			# ylocator=[0.1,0.2,0.5,1,2,5,10,20],
			# xlocator=[10,20,50,100,200,500,1000]/100
		)
	end
	
	ctsfc_1 = atsfc[1].pcolormesh(
		sizelist/100,depthlist,tsfcmat',cmap="RdBu_r",
		levels=vcat(-20,-15,-10,-5:5,10,15,20),extend="both"
	)
	ctsfc_2 = atsfc[2].pcolormesh(
		sizelist/100,depthlist,tmaxmat',cmap="RdBu_r",
		levels=vcat(-20,-15,-10,-5:5,10,15,20),extend="both"
	)
	ctsfc_3 = atsfc[3].pcolormesh(
		sizelist/100,depthlist,ttmxmat',levels=6:0.5:18,cmap="RdBu",extend="both"
	)
	ctsfc_txt = atsfc[3].contour(
		sizelist/100,depthlist,ttmxmat',levels=6:3:18,c="k",lw=0.5,linestyle="--"
	); ; atsfc[3].clabel(ctsfc_txt,inline=true,fontsize=10)
	ctsfc_4 = atsfc[4].pcolormesh(
		sizelist/100,depthlist,tampmat',cmap="RdBu_r",
		levels=vcat(-20,-15,-10,-5:5,10,15,20),extend="both"
	)
	ctsfc_5 = atsfc[5].pcolormesh(
		sizelist/100,depthlist,tminmat',cmap="RdBu_r",
		levels=vcat(-20,-15,-10,-5:5,10,15,20),extend="both"
	)
	ctsfc_6 = atsfc[6].pcolormesh(
		sizelist/100,depthlist,ttmnmat',levels=6:0.5:18,cmap="RdBu",extend="both"
	)

	atsfc[1].format(ltitle=L"(a) $\mu(T_s) - 300$",suptitle="Surface Temperature Climatology")
	atsfc[2].format(ltitle=L"(c) max($T_s$) - 300")
	atsfc[3].format(ltitle=L"(e) $\theta$(max($T_s$))")
	atsfc[4].format(ltitle=L"(b) max($T_s$) - min($T_s$)")
	atsfc[5].format(ltitle=L"(d) min($T_s$) - 300")
	atsfc[6].format(ltitle=L"(f) $\theta$(min($P$))")

	ftsfc.colorbar(ctsfc_1,length=0.8,loc="l",label="Temperature / K")
	ftsfc.colorbar(ctsfc_6,length=0.8,label="Hour of Day",locator=0:3:24,minorlocator=[])
	
	ftsfc.savefig(plotsdir("06b-islandsize-tsurface.png"),transparent=false,dpi=400)
	load(plotsdir("06b-islandsize-tsurface.png"))
end

# ╔═╡ 03cf4cf1-9d98-4f0e-83c3-036fe0ecae2b
md"
### B. Outgoing Longwave Radiation Statistics
"

# ╔═╡ 21246adc-d19f-4eb6-8bb5-53f3384e640e
begin
	olrμmat = zeros(nsize,ndepth) * NaN
	olrAmat = zeros(nsize,ndepth) * NaN
	olmxmat = zeros(nsize,ndepth) * NaN
	olmnmat = zeros(nsize,ndepth) * NaN
	otmxmat = zeros(nsize,ndepth) * NaN
	otmnmat = zeros(nsize,ndepth) * NaN
	md"Preallocation of OLR arrays"
end

# ╔═╡ e0d48ddb-0ddc-4e25-90eb-c38f183d2d56
for idepth in 1 : ndepth
	depthstr = @sprintf("%05.2f",depthlist[idepth])
	depthstr = replace(depthstr,"."=>"d")
	for isize in 1 : nsize
		sizestr = @sprintf("%04d",sizelist[isize])

		config = "size$(sizestr)km-depth$(depthstr)m"

		fnc = outstatname("DGW","IslandSize3D",config)
		if isfile(fnc)
			ds  = NCDataset(fnc)
			if length(ds["time"][:]) >= 4800
				time = ds["time"][1:48] .- floor(ds["time"][1])
				olr  = ds["LWNTOA"][(48*50+1):end]
				olrμmat[isize,idepth] = mean(olr)
				otmx = dropdims(mean(reshape(olr,48,:),dims=2),dims=2)
				otmxmat[isize,idepth] = time[argmax(otmx)] * 24
				otmnmat[isize,idepth] = time[argmin(otmx)] * 24
				olrAmat[isize,idepth] = maximum(otmx) - minimum(otmx)
				olmxmat[isize,idepth] = maximum(otmx)
				olmnmat[isize,idepth] = minimum(otmx)
			end
			close(ds)
		end

	end
end

# ╔═╡ bb611707-a6a1-4009-bdd4-0bb1923ee528
begin
	pplt.close(); folr,aolr = pplt.subplots(axwidth=1.5,nrows=2,ncols=3,)
	
	colr_1 = aolr[1].pcolormesh(
		sizelist/100,depthlist,olrμmat',
		levels=175:5:250,extend="both",cmap="gray_r"
	)
	colr_2 = aolr[2].pcolormesh(
		sizelist/100,depthlist,olmxmat',levels=200:5:300,cmap="gray_r",extend="both"
	)
	colr_3 = aolr[3].pcolormesh(
		sizelist/100,depthlist,otmxmat',levels=0:24,cmap="RdBu",
	)
	colr_4 = aolr[4].pcolormesh(
		sizelist/100,depthlist,olrAmat',levels=20:20:200,cmap="gray",extend="both"
	)
	colr_5 = aolr[5].pcolormesh(
		sizelist/100,depthlist,olmnmat',levels=75:5:200,cmap="gray_r",extend="both"
	)
	colr_6 = aolr[6].pcolormesh(
		sizelist/100,depthlist,otmnmat',levels=0:24,cmap="RdBu",
	)
	
	aolr[1].format(ltitle=L"(a) $\mu$(OLR) / W m$^{-2}$",suptitle="OLR Climatology")
	aolr[2].format(ltitle=L"(c) max(OLR) / W m$^{-2}$")
	aolr[3].format(ltitle=L"(e) $\theta$(max(OLR))")
	aolr[4].format(ltitle="(b) max(OLR) - min(OLR)")
	aolr[5].format(ltitle=L"(d) min(OLR) / W m$^{-2}$")
	aolr[6].format(ltitle=L"(f) $\theta$(min(OLR))")

	for ax in aolr
		ax.format(
			xscale="log",xlim=(0.05,20),ylabel="Mixed Layer Depth / m",
			yscale="log",ylim=(0.05,20),xlabel="Island Radius / km",
			# ylocator=[0.1,0.2,0.5,1,2,5,10,],
			# xlocator=[10,20,50,100,200,500,1000]
		)
	end

	aolr[1].colorbar(loc="r",colr_1)
	aolr[2].colorbar(loc="r",colr_2)
	aolr[4].colorbar(loc="r",colr_4)
	aolr[5].colorbar(loc="r",colr_5)
	folr.colorbar(loc="r",colr_3,locator=0:3:24,length=0.8,minorlocator=[],label="Hour of Day")
	
	folr.savefig(plotsdir("06b-islandsize-olr.png"),transparent=false,dpi=400)
	load(plotsdir("06b-islandsize-olr.png"))
end

# ╔═╡ Cell order:
# ╟─db521478-61e1-4184-960b-7cab95a48b50
# ╟─c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
# ╟─8d72de06-476c-4bcf-99b3-a9b469fac93d
# ╟─af993da1-6181-4c26-8531-e8537ae629d9
# ╟─a99d31a1-4ef8-4cc9-894c-342fd5424e36
# ╟─37a716a7-2425-4f9d-960f-a0be3744a223
# ╟─1455f008-06c8-4f79-a852-ca7d4a324fe8
# ╠═bb4b054f-87ca-4076-97be-207097ceb79e
# ╠═1e9f2a5f-3a84-43e3-9ff8-a8adce7c8077
# ╟─4ab80acb-afe4-45ad-8698-1131a77c5c79
# ╟─99e126ae-c746-4e83-8fc5-c09a5dd7980d
# ╟─03cf4cf1-9d98-4f0e-83c3-036fe0ecae2b
# ╟─21246adc-d19f-4eb6-8bb5-53f3384e640e
# ╟─e0d48ddb-0ddc-4e25-90eb-c38f183d2d56
# ╟─bb611707-a6a1-4009-bdd4-0bb1923ee528
