# using Dates
# using DrWatson
# using Glob
# using NCDatasets
using Statistics

function coarse1D(var;dx=1)

    return dropdims(mean(reshape(var,dx,:),dims=1),dims=1)

end

function coarse2D(var;dx=1,dy=1)

    nx,ny = size(var)
    nx /= dx; nx = Int(nx)
    ny /= dy; ny = Int(ny)
    return dropdims(mean(reshape(var,dx,nx,dy,ny),dims=(1,3)),dims=(1,3))

end

function coarse3D(var;dx=1,dy=1,dz=1)

    nx,ny,nz = size(var)
    nx /= dx; nx = Int(nx)
    ny /= dy; ny = Int(ny)
    nz /= dz; nz = Int(nz)
    return dropdims(mean(reshape(var,dx,nx,dy,ny,dz,nz),dims=(1,3,5)),dims=(1,3,5))

end