using Printf

lsfinit(nvert::Integer) = zeros(nvert,7)

function lsfprint(flsf::AbstractString,lsf::Array{<:Real,2},p::Real)

    nz = size(lsf,1)

    open(flsf,"w") do io
        @printf(io,"  z[m] p[mb] tpls[K/s] qls[kg/kg/s] uls_hor vls_hor wls[m/s]\n")
        @printf(io,"%10.2f, %10d, %10.2f\n",0.00,nz,p)
    end

    open(flsf,"a") do io
        for iz = 1 : nz
            @printf(
                io,"%16.8f, %16.8f, %16.8e, %16.8e, %16.8f, %16.8f, %16.8f\n",
                lsf[iz,1],lsf[iz,2],lsf[iz,3],
                lsf[iz,4],lsf[iz,5],lsf[iz,6],lsf[iz,7]
            )
        end
    end

    open(flsf,"a") do io
        @printf(io,"%10.2f, %10d, %10.2f\n",10000.00,nz,p)
    end

    open(flsf,"a") do io
        for iz = 1 : nz
            @printf(
                io,"%16.8f, %16.8f, %16.8e, %16.8e, %16.8f, %16.8f, %16.8f\n",
                lsf[iz,1],lsf[iz,2],lsf[iz,3],
                lsf[iz,4],lsf[iz,5],lsf[iz,6],lsf[iz,7]
            )
        end
    end

end
