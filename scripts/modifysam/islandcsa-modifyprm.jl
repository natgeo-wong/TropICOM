using DrWatson
@quickactivate "TropICOM"
using Printf

include(srcdir("common.jl"))

oprm = rundir("modifysam","prmtemplates","IslandCSA","position.prm")
open(oprm,"r") do rprm
    oldstr = read(rprm,String)
    for position in 384 : 24 : 1536

        posstr = @sprintf("%04d",position)
        mkpath(prmdir("IslandCSA"))
        nprm = prmdir("IslandCSA","$(posstr).prm")
        open(nprm,"w") do fprm
            newstr = replace(oldstr,"[position]"=>posstr)
            write(fprm,newstr)
        end
    end
end