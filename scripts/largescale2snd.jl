using DrWatson
@quickactivate "TroPrecLS"

using NCDatasets

include(srcdir("samsnd.jl"))

ds = NCDataset(datadir("SAMvertprofile.nc"))
z  = ds["z"][:]; nz = length(z)
cr = ds["bin"][:]

snddata = zeros(nz,6)

for icrh = 20 : 5 : 95

    ii = argmin(abs.(cr.-icrh))

    pp = ds["level"][:] .+ ds["pp"][ii,:] / 100
    qv = ds["q"][ii,:]
    ta = ds["t"][ii,:]
    pt = ta ./ (pp/1000).^(287/1004)

    snddata[:,1] .= z
    snddata[:,2] .= pp
    snddata[:,3] .= pt
    snddata[:,4] .= qv

    createsndmean("colrelhum$icrh",snddata;psfc=1009.32)

end

close(ds);