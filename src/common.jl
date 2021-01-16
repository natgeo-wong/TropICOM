using Dates
using Statistics
using GeoRegions

prect(N::Real,S::Real,W::Real,E::Real) = [W,E,E,W,W],[S,S,N,N,S]

yrmo2str(date::TimeType) = Dates.format(date,dateformat"yyyymm")

function ncoffsetscale(data::AbstractArray)

    dmax = maximum(data[.!ismissing.(data)])
    dmin = minimum(data[.!ismissing.(data)])
    scale = (dmax-dmin) / 65533;
    offset = (dmax+dmin-scale) / 2;

    return scale,offset

end

function bindatasfclnd(coords,bins,var,lon,lat,lsm)

    tlon,tlat,rinfo = regiongridvec(coords,lon,lat);
    rvar = regionextractgrid(var,rinfo)
    rlsm = regionextractgrid(lsm,rinfo)
    rwgt = ones(size(rlsm)) .* cosd.(reshape(tlat,1,:))

    lvar = rvar[rlsm.>0.5];
	lvar = lvar[.!ismissing.(lvar)]; lvar = lvar[.!isnan.(lvar)]
    lbin = fit(Histogram,lvar,bins).weights;
	lbin = lbin ./ sum(lbin) * (length(bins) - 1)

    rvar = rvar .* cosd.(reshape(tlat,1,:))
    lvar = rvar[rlsm.>0.5];
	lvar = lvar[.!ismissing.(lvar)]; lvar = lvar[.!isnan.(lvar)]
    lvar = lvar / mean(rwgt[rlsm.>0.5])

    return lbin,mean(lvar)

end

function bindatasfcsea(coords,bins,var,lon,lat,lsm)

    tlon,tlat,rinfo = regiongridvec(coords,lon,lat);
    rvar = regionextractgrid(var,rinfo)
    rlsm = regionextractgrid(lsm,rinfo)
    rwgt = ones(size(rlsm)) .* cosd.(reshape(tlat,1,:))

    svar = rvar[rlsm.<0.5];
	svar = svar[.!ismissing.(svar)]; svar = svar[.!isnan.(svar)]
    sbin = fit(Histogram,svar,bins).weights;
	sbin = sbin ./ sum(sbin) * (length(bins) - 1)

    rvar = rvar .* cosd.(reshape(tlat,1,:))
    svar = rvar[rlsm.<0.5];
	svar = svar[.!ismissing.(svar)]; svar = svar[.!isnan.(svar)]
    svar = svar / mean(rwgt[rlsm.<0.5])

    return sbin,mean(svar)

end

function getmean(coords,var,lon,lat,nlvl,lsm)

    tlon,tlat,rinfo = regiongridvec(coords,lon,lat);
    rvar = regionextractgrid(var,rinfo)
    rlsm = regionextractgrid(lsm,rinfo)
    sprf = zeros(nlvl); lprf = zeros(nlvl)

    for ilvl = 1 : nlvl

        ivar = rvar[:,:,ilvl];  ivar = ivar .* cosd.(reshape(tlat,1,:))
        lavg = ivar[rlsm.>0.5]; lavg = lavg[.!ismissing.(lavg)]; lavg = lavg[.!isnan.(lavg)]
        savg = ivar[rlsm.<0.5]; savg = savg[.!ismissing.(savg)]; savg = savg[.!isnan.(savg)]
        lprf[ilvl] = mean(lavg)
        sprf[ilvl] = mean(savg)

    end

    return lprf,sprf

end
