### A Pluto.jl notebook ###
# v0.19.9

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
	using NumericalIntegration
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ 47218e37-00d8-4a0a-9f31-c3722636de8a
begin
	ds = NCDataset(datadir("Slab02d00/OUT_STAT/DGW_TroPrecLS-Slab02d00-damping01d0.nc"))
	t  = ds["time"][:]
	pp = ds["p"][:] * 100
	pw = ds["PW"][:]
	qv = ds["QV"][:] / 1000
	qo = ds["QVOBS"][:] / 1000
	ta = ds["TABS"][:]
	to = ds["TABSOBS"][:]
	pr = ds["PREC"][:]
	ts = ds["SST"][:]
	w  = ds["WWTG"][:]
	close(ds)
end

# ╔═╡ fefd739b-f0ea-44dc-98ba-10b3496243c7
rh = calcrh(qv,ta,pp)

# ╔═╡ bbe06590-36c2-4686-89d0-79da8cb3da71
qv

# ╔═╡ bd95b4d4-4258-4858-8d41-0d68ccd93437
function calcswp(RH,QV,P)

	pvec = reverse(P); nt = size(RH,2)
	QVsat = zeros(length(pvec))
	swp = zeros(nt)

    for it = 1 : nt
		QVsat .= reverse(QV[:,it]) ./ reverse(RH[:,it]/100)
		swp[it] = integrate(pvec,QVsat) / 9.81 / 1000
    end

	return swp

end

# ╔═╡ 5a938176-9a86-48f2-b682-7f37d8c26b98
sw = calcswp(rh,qv,pp) * 1000

# ╔═╡ 40306fcc-9ad5-4295-8ca9-3a447787a155
size(ta)

# ╔═╡ 7b48d3e0-f56f-41f1-bac3-f35145cd3d70
ta2 = dropdims(mean(reshape(ta[:,12001:14400],64,96,:),dims=3),dims=3)

# ╔═╡ ab6ad48d-eddf-4d81-b636-ebf645508697
begin
	pplt.close(); fig,axs = pplt.subplots(
		[[1,1,1,1,1,1,4],[2,2,2,2,2,2,0],[3,3,3,3,3,3,0]],
		aspect=3,axwidth=4)
	
	# c = axs[1].contourf(t,pp/100,rh,cmap="Blues",levels=50:10:100,extend="both")
	c = axs[1].contourf(t,pp/100,ta.-ta[:,1],cmap="RdBu_r",extend="both")
	axs[1].format(ylim=(1000,25),yscale="log")
	axs[2].plot(t,pw ./sw)
	axs[3].plot(t,pr)
	axs[4].plot(mean(ta,dims=2),pp/100)

	axs[3].format(xlim=(250,300),ylim=(0,150))

	fig.colorbar(c)
	fig.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ 0aeb9695-06c9-4e7a-9135-7e9248de93b1
minimum(pr)

# ╔═╡ f939c8a3-cf7d-4e93-9f97-bf89e5ad28fe
begin
	pr2 = dropdims(mean(reshape(pr,2,:),dims=1),dims=1)
	pr2 = pr2[1201:7200]
end

# ╔═╡ 5f301d23-d974-41d0-9e6e-78683115875b
crh = dropdims(mean(reshape(pw./sw * 100,2,:),dims=1),dims=1)

# ╔═╡ 141d4128-6f85-4295-9c84-af75cc3f73f8
crh2 = crh[1201:7200]

# ╔═╡ c054ee62-f199-4f31-a23d-995c578016e4
bcrh = vcat(74.5:100.5)

# ╔═╡ 34679b9b-8056-4e8c-be17-32f1eb50a42c
bcrh_p = 75:100

# ╔═╡ 76c30632-d07b-46f0-9dcd-9592608e6c73
bprcp = zeros(length(bcrh_p))

# ╔═╡ ebd3df23-3a73-4048-bd0d-f8c11d058823
begin
	for ii = 1 : length(bcrh_p)
		bprcp[ii] = mean(pr2[(crh2.>bcrh[ii]) .& (crh2.<bcrh[ii+1])])
	end
	pplt.close(); f2,a2 = pplt.subplots(aspect=2,axwidth=4)
	
	a2[1].scatter(bcrh_p,bprcp / 24)
	a2[1].format(yscale="log",xlim=(60,100),ylim=(0.0005,20))
	
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
# ╟─8d72de06-476c-4bcf-99b3-a9b469fac93d
# ╠═47218e37-00d8-4a0a-9f31-c3722636de8a
# ╠═fefd739b-f0ea-44dc-98ba-10b3496243c7
# ╠═bbe06590-36c2-4686-89d0-79da8cb3da71
# ╠═bd95b4d4-4258-4858-8d41-0d68ccd93437
# ╠═5a938176-9a86-48f2-b682-7f37d8c26b98
# ╠═40306fcc-9ad5-4295-8ca9-3a447787a155
# ╠═7b48d3e0-f56f-41f1-bac3-f35145cd3d70
# ╠═ab6ad48d-eddf-4d81-b636-ebf645508697
# ╠═0aeb9695-06c9-4e7a-9135-7e9248de93b1
# ╠═f939c8a3-cf7d-4e93-9f97-bf89e5ad28fe
# ╠═5f301d23-d974-41d0-9e6e-78683115875b
# ╠═141d4128-6f85-4295-9c84-af75cc3f73f8
# ╠═c054ee62-f199-4f31-a23d-995c578016e4
# ╠═34679b9b-8056-4e8c-be17-32f1eb50a42c
# ╠═76c30632-d07b-46f0-9dcd-9592608e6c73
# ╟─ebd3df23-3a73-4048-bd0d-f8c11d058823
