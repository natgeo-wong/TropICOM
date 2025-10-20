using DrWatson
@quickactivate "TropICOM"
using Printf

include(srcdir("common.jl"))

slist = [
    5,5*sqrt(2),10,10*sqrt(2),20,20*sqrt(2.5),
    50,50*sqrt(2),100,100*sqrt(2),200,200*sqrt(2.5),
    500,500*sqrt(2),1000,1000*sqrt(2),2000
] / 1000
dlist = [
    0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
    0.1,0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
    1.,1*sqrt(2),2.,2*sqrt(2.5),5.,5*sqrt(2),
    10.,10*sqrt(2),20.,20*sqrt(2.5),50.,
]

oprm = rundir("modifysam","prmtemplates","IslandDGW","spinup.prm")
open(oprm,"r") do rprm
    oldstr = read(rprm,String)
    for depth in dlist, islandsize in slist
        expname = dampingstrprnt(islandsize)
        runname = depthstrprnt(depth)

        mkpath(prmdir("IslandDGW",expname))
        nprm = prmdir("IslandDGW",expname,"spinup-$(runname).prm")
        open(nprm,"w") do fprm
            newstr = replace(oldstr,"[expname]"=>expname)
            newstr = replace(newstr,"[runname]"=>runname)
            newstr = replace(newstr,"[depth]"=>@sprintf("%7e",depth))
            newstr = replace(newstr,"[am]"=>@sprintf("%7e",islandsize))
            write(fprm,newstr)
        end
    end
end


oprm = rundir("modifysam","prmtemplates","IslandDGW","statistics.prm")
open(oprm,"r") do rprm
    oldstr = read(rprm,String)
    for depth in dlist, islandsize in slist
        expname = dampingstrprnt(islandsize)
        runname = depthstrprnt(depth)

        mkpath(prmdir("IslandDGW",expname))
        nprm = prmdir("IslandDGW",expname,"$(runname).prm")
        open(nprm,"w") do fprm
            newstr = replace(oldstr,"[expname]"=>expname)
            newstr = replace(newstr,"[runname]"=>runname)
            newstr = replace(newstr,"[depth]"=>@sprintf("%7e",depth))
            newstr = replace(newstr,"[am]"=>@sprintf("%7e",islandsize))
            write(fprm,newstr)
        end
    end
end