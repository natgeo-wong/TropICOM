using NCDatasets
using Statistics

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

include(srcdir("sam.jl"))

function plottair(
    tvec::Vector{<:Real}, height::Vector{<:Real}, data::AbstractArray{<:Real,2},
    axsii::PyObject
)

    c = axsii.contourf(
        tvec,height,data,cmap="RdBu_r",extend="both",
        norm="segmented",levels=[-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5]
    ); axsii.colorbar(c,loc="b")

    axsii.format(
        xlim=(0,24),xlocator=[0,3,6,9,12,15,18,21,24],
        title="Deviation from Mean",rtitle="Air Temperature / K"
    )

end

function plotcld(
    tvec::Vector{<:Real}, height::Vector{<:Real}, data::AbstractArray{<:Real,2},
    axsii::PyObject
)

    c = axsii.contourf(
        tvec,height,data.-mean(data,dims=2),
        #cmap="Blues",norm="segmented",levels=[0,1,2,5,10,20,50,80,90,95,99,100]
        cmap="drywet",levels=[-10,-5,-2,-1,-0.5,-0.2,0.2,0.5,1,2,5,10],extend="both"
    ); axsii.colorbar(c,loc="r")

    axsii.format(
        xlim=(0,24),xlocator=[0,3,6,9,12,15,18,21,24],xlabel="Hour of Day",
        rtitle="Cloud Cover Fraction / %"
    )

end

function plotrhum(
    tvec::Vector{<:Real}, height::Vector{<:Real}, data::AbstractArray{<:Real,2},
    axsii::PyObject
)

    c = axsii.contourf(
        tvec,height,data.-mean(data,dims=2),
        #cmap="Blues",norm="segmented",levels=[0,1,2,5,10,20,50,80,90,95,99,100]
        cmap="drywet",levels=vcat(-5:-1,1:5),extend="both"
    ); axsii.colorbar(c,loc="r")

    axsii.format(
        xlim=(0,24),xlocator=[0,3,6,9,12,15,18,21,24],xlabel="Hour of Day",
        rtitle="Relative Humidity / %"
    )

end

function plotwwtg(
    tvec::Vector{<:Real}, height::Vector{<:Real}, data::AbstractArray{<:Real,2},
    axsii::PyObject
)

    c = axsii.contourf(
        tvec,height,data,cmap="RdBu_r",extend="both",
        norm="segmented",levels=[-10,-5,-2,-1,-0.5,0.5,1,2,5,10]/100
    ); axsii.colorbar(c,loc="r")

    axsii.format(
        xlim=(0,24),xlocator=[0,3,6,9,12,15,18,21,24],xlabel="Hour of Day",
        rtitle=L"Weak Temperature Gradient / m s$^{-1}$"
    )

end

function plotqtend(
    tvec::Vector{<:Real}, height::Vector{<:Real}, data::AbstractArray{<:Real,2},
    axsii::PyObject
)

    c = axsii.contourf(
        tvec,height,data,cmap="drywet",extend="both",
        norm="segmented",levels=[-50,-20,-10,-5,-2,-1,1,2,5,10,20,50]/100
    ); axsii.colorbar(c,loc="r")

    axsii.format(
        xlim=(0,24),xlocator=[0,3,6,9,12,15,18,21,24],xlabel="Hour of Day",
        rtitle=L"Q Tendency due to W Large Scale Advect / g/kg hr$^{-1}$"
    )

end

function retrieve3Dvar(
    experiment::AbstractString, config::AbstractString, plotWTG::Bool=false
)

    rce = NCDataset(datadir(joinpath(
        experiment,config,"OUT_STAT",
        "RCE_TroPrecLS-$(experiment).nc"
    )))

    p    = rce["p"][:]; t = rce["time"][:] .- 80; nt = length(t)
    cld  = rce["CLD"][:]*100; qv = rce["QV"][:]/10;
    tair = rce["TABS"][:]; qvtd = rce["QVTEND"]/24; wwtg = 0;

    rh = zeros(size(qv,1),nt)
    for it = 1 : nt, ilvl = 1 : length(p)

        rh[ilvl,it] = qv[ilvl,it] / temp2qsat(tair[ilvl,it],p[ilvl]*100)

    end

    if plotWTG; wwtg = rce["WWTG"][:] end

    return p,t,cld,rh,tair,qvtd,wwtg

end

function diurnal(data::AbstractArray{<:Real,2},tstep::Integer,tshift::Integer)

    nz = size(data,1)
    data = dropdims(mean(reshape(data,nz,tstep,:),dims=3),dims=3);
    data = circshift(data,(0,tshift))

    return cat(dims=2,data[:,end],data,data[:,1])

end

function plot3Ddiurnal(
    experiment::AbstractString, config::AbstractString;
    plotWTG::Bool=true, days::Integer=100
)

    z,t,cld,rh,tair,qvtd,wwtg = retrieve3Dvar(experiment,config,plotWTG)

    t,tstep,tshift = t2d(t); beg = days*tstep - 1
    cld  = diurnal(cld[:,(end-beg):end],tstep,tshift);  mcld  = mean(cld,dims=2)
    rh   = diurnal(rh[:,(end-beg):end],tstep,tshift);   mrh   = mean(rh,dims=2)
    tair = diurnal(tair[:,(end-beg):end],tstep,tshift); mtair = mean(tair,dims=2)
    qvtd = diurnal(qvtd[:,(end-beg):end],tstep,tshift); mqvtd = mean(qvtd,dims=2)
    if plotWTG
        wwtg  = diurnal(wwtg[:,(end-beg):end],tstep,tshift)
        mwwtg = mean(wwtg,dims=2)
    end

    pplt.close(); arr = [[0,0,1,2,2,2,0,0],[3,4,4,4,5,6,6,6],[7,8,8,8,9,10,10,10]];
    f,axs = pplt.subplots(arr,axwidth=1,aspect=0.5,sharex=1,wspace=[0,0,0,nothing,0,0,0])

    axs[1].plot(dropdims(mean(tair,dims=2),dims=2),z,lw=1)
    axs[1].format(ylabel="Pressure / hPa",title="Daily Mean",xlim=(180,300));
    plottair(t,z,tair .- mtair,axs[2])

    axs[3].plot(dropdims(mean(cld,dims=2),dims=2),z,lw=1);
    axs[3].format(ylabel="Pressure / hPa",xlim=(0,75))
    plotcld(t,z,cld,axs[4])

    axs[5].plot(dropdims(mean(rh,dims=2),dims=2),z,lw=1);
    axs[5].format(ylabel="Pressure / hPa",xlim=(0,100))
    plotrhum(t,z,rh,axs[6])

    if plotWTG
        axs[7].plot(dropdims(mean(wwtg,dims=2),dims=2),z,lw=1);
        axs[7].format(ylabel="Pressure / hPa",xlim=(-0.1,0.1))
        plotwwtg(t,z,wwtg,axs[8])
    else
        axs[8].format(
            xlim=(0,24),xlocator=[0,3,6,9,12,15,18,21,24],xlabel="Hour of Day",
            rtitle=L"Weak Temperature Gradient / m s$^{-1}$"
        )
    end

    axs[9].plot(dropdims(mean(qvtd,dims=2),dims=2),z,lw=1);
    axs[9].format(ylabel="Pressure / hPa",xlim=(-0.4,0.4))
    plotqtend(t,z,qvtd,axs[10])

    axs[1].format(suptitle="$(uppercase(experiment)) | $(uppercase(config))");

    for ii = 1 : 5; axs[ii].format(ylim=(1010,50)) end

    if !isdir(plotsdir("SAM_STAT-3D")); mkpath(plotsdir("SAM_STAT-3D")) end
    f.savefig(plotsdir("SAM_STAT-3D/$experiment-$config.png"),transparent=false,dpi=200)

end
