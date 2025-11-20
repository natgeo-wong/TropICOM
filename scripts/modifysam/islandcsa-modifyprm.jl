using DrWatson
@quickactivate "TropICOM"
using Printf

include(srcdir("common.jl"))

dSSTv  = [0,0.02,0.05,0.1,0.2,0.5,1,2,5]
depthv = [0.02,0.1,1,10,50]
sizev  = [10,20,50,100]

oprm = rundir("modifysam","prmtemplates","IslandCSA","statistics.prm")
open(oprm,"r") do rprm
    oldstr = read(rprm,String)
    for dSST in dSSTv, depth in depthv, size in sizev, position in 384 : 96 : 1152

        nprm = prmdir(
            "IslandCSA","dSST[dSST]K","mld[depth]m",
            "size[size]km-x[position]km.prm"
        )
        nprm = replace(nprm,"[dSST]"=>@sprintf("%04.2f",dSST))
        nprm = replace(nprm,"[depth]"=>@sprintf("%05.2f",depth))
        nprm = replace(nprm,"[size]"=>@sprintf("%03d",size))
        nprm = replace(nprm,"[position]"=>@sprintf("%04d",position))

        mkpath(dirname(nprm))

        open(nprm,"w") do fprm
            newstr = replace(oldstr,"[position]"=>@sprintf("%04d",position))
            newstr = replace(newstr,"[dSST]"=>@sprintf("%.2f",dSST))
            newstr = replace(newstr,"[dSST/2]"=>@sprintf("%.2f",dSST/2))
            newstr = replace(newstr,"[depth]"=>@sprintf("%.2f",depth))
            newstr = replace(newstr,"[size]"=>@sprintf("%03d",size))
            write(fprm,newstr)
        end

    end
end

oprm = rundir("modifysam","prmtemplates","IslandCSA","spinup.prm")
open(oprm,"r") do rprm
    oldstr = read(rprm,String)
    for dSST in dSSTv, depth in depthv, size in sizev, position in 384 : 96 : 1152

        nprm = prmdir(
            "IslandCSA","dSST[dSST]K","mld[depth]m",
            "spinup-size[size]km-x[position]km.prm"
        )
        nprm = replace(nprm,"[dSST]"=>@sprintf("%04.2f",dSST))
        nprm = replace(nprm,"[depth]"=>@sprintf("%05.2f",depth))
        nprm = replace(nprm,"[size]"=>@sprintf("%03d",size))
        nprm = replace(nprm,"[position]"=>@sprintf("%04d",position))

        mkpath(dirname(nprm))

        open(nprm,"w") do fprm
            newstr = replace(oldstr,"[position]"=>@sprintf("%04d",position))
            newstr = replace(newstr,"[dSST]"=>@sprintf("%.2f",dSST))
            newstr = replace(newstr,"[dSST/2]"=>@sprintf("%.2f",dSST/2))
            newstr = replace(newstr,"[depth]"=>@sprintf("%.2f",depth))
            newstr = replace(newstr,"[size]"=>@sprintf("%03d",size))
            write(fprm,newstr)
        end

    end
end