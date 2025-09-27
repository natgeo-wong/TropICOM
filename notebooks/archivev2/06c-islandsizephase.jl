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
# 6c. Phases of the Diurnal Cycle

Text
"

# ╔═╡ af993da1-6181-4c26-8531-e8537ae629d9
TableOfContents()

# ╔═╡ a99d31a1-4ef8-4cc9-894c-342fd5424e36
begin
	sizelist = [5,10,20,50,100,200,500,1000,2000]
	nsize = length(sizelist)
	md"Loading vector of Island sizes in km ..."
end

# ╔═╡ 37a716a7-2425-4f9d-960f-a0be3744a223
begin
	depthlist = [0.05,0.1,0.2,0.5,1.,2.,5.,10.,20.,]
	ndepth = length(depthlist)
	md"Loading vector of mixed-layer depth in m ..."
end

# ╔═╡ e8499402-439d-4de2-8be2-1a43b2b506d5
begin
	blues_WTG = pplt.get_colors("roma",(nsize))
	md"Loading Colors ..."
end

# ╔═╡ 18245629-9a24-4d97-a4c3-80877b19929b
function retrievemaxtime(varID,sizelist,depthlist)

	nsize  = length(sizelist)
	ndepth = length(depthlist)
	d2D  = zeros(192,nsize,ndepth)
	vtmp = zeros(96)

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
					var = reshape(var,96,:); ndy = size(var,2)
					for idy = 1 : ndy
						vtmp .= var[:,idy]
						imax = argmax(vtmp)
						if vtmp[imax] > 0.25*24
							d2D[(2*imax-1):(2*imax),is,id] .+= 1
						end
					end
					d2D[:,is,id] .= d2D[:,is,id] ./ sum(d2D[:,is,id]) * 96
				end
				close(ds)
			end
		end
	end

	return d2D

end

# ╔═╡ 1455f008-06c8-4f79-a852-ca7d4a324fe8
md"
### A. Precipitation and Surface Temperature Statistics
"

# ╔═╡ 5c04acab-2807-4526-858f-7601b368130b
begin
	prcpθ = retrievemaxtime("PREC",sizelist,depthlist)
	md"Binning precipitation maximum times"
end

# ╔═╡ 8cacf996-b900-41e5-a184-2a79d785d1fa
begin
	pplt.close()
	fθ,aθ = pplt.subplots([[1,4,7],[2,5,8],[3,6,9]],proj="polar",axwidth=1.5
	);

	θvec = zeros(2,96)
	θvec[1,:] = collect(0:0.25:23.75)/24*2*pi
	θvec[2,:] = collect(0.25:0.25:24)/24*2*pi
	θvec = θvec[:]
	θvec = vcat(θvec,2*pi)

	for ii in 1 : 9
		for is in nsize : -1 : 1
	
			ivec = sqrt.(prcpθ[:,is,ii])
			ivec = vcat(ivec,ivec[1])
			if ii != 8
				aθ[ii].plot(θvec,ivec,c=blues_WTG[is],lw=1)
			else
				aθ[ii].plot(
					θvec,ivec,c=blues_WTG[is],lw=1,
					label=L"a_m\lambda^2" * " = $(@sprintf("%.1e",sizelist[is]*1000)) " * L"km$^2$ day$^{-1}$",
					legend="r",legend_kw=Dict("ncol"=>1,"frame"=>false)
				)
			end
	
		end
		aθ[ii].format(ultitle="MLD = $(depthlist[ii]) m")
	end
	

	for ax in aθ
		ax.format(
			theta0="N",thetaformatter="tau",
			rlim=(0,5),rlabelpos=135,rlocator=1:4,thetadir=-1,
			suptitle=L"$\theta$ / Fraction of Day"
		)
		ax.set_xticks((0:24)/24*2*pi)
		ax.set_xticklabels(vcat(collect(0:23),""))
	end

	fθ.savefig(plotsdir("06c-islandsize-rainphase.png"),transparent=false,dpi=400)
	load(plotsdir("06c-islandsize-rainphase.png"))
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
# ╟─1455f008-06c8-4f79-a852-ca7d4a324fe8
# ╟─5c04acab-2807-4526-858f-7601b368130b
# ╟─8cacf996-b900-41e5-a184-2a79d785d1fa
