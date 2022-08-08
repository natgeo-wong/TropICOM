using DrWatson
@quickactivate "TroPrecLS"

using NCDatasets

include(srcdir("samsnd.jl"))

ds = NCDataset(datadir("SAMvertprofile.nc"))
z  = ds["z"][:]; nz = length(z)
cr = ds["bin"][:]

snddata = zeros(nz,7)

for icrh = 20 : 5 : 95

    ii = argmin(abs.(cr.-icrh))

    pp = ds["level"][:] .+ ds["pp"][ii,:] / 100
    wa = ds["wa"][ii,:]

    snddata[:,1] .= z
    snddata[:,2] .= pp
    snddata[:,7] .= wa

    createlsfmean("colrelhum$icrh",snddata;psfc=1009.32)

end

close(ds);