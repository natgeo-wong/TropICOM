using DrWatson
@quickactivate "TropICOM"
using Printf

include(srcdir("common.jl"))

case    = "NEUTRAL"
email   = ""
doBuild = true

slist = [
    2,2*sqrt(2.5),5,5*sqrt(2),10,10*sqrt(2),20,20*sqrt(2.5),
    50,50*sqrt(2),100,100*sqrt(2),200,200*sqrt(2.5),
    500,500*sqrt(2),1000,1000*sqrt(2),2000,2000*sqrt(2.5),5000
] / 1000
dlist = [
    0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
    0.1,0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
    1.,1*sqrt(2),2.,2*sqrt(2.5),5.,5*sqrt(2),
    10.,10*sqrt(2),20.,20*sqrt(2.5),50.,
]

orun  = rundir("modifysam","runtemplates","modelrun_islandwtg.sh")
open(orun,"r") do rrun
    oldstr = read(rrun,String)
    for depth in dlist, islandsize in slist
        expname = dampingstrprnt(islandsize)
        runname = depthstrprnt(depth)

        nrun = rundir("IslandDGW",expname,"spinup-$(runname).sh")
        open(nrun,"w") do fprm
            newstr = replace(oldstr,"[email]"   => email)
            newstr = replace(newstr,"[dirname]" => projectdir())
            newstr = replace(newstr,"[project]" => "IslandDGW")
            newstr = replace(newstr,"[expname]" => expname)
            newstr = replace(newstr,"[runname]" => "spinup-$(runname)")
            newstr = replace(newstr,"[sndname]" => case)
            newstr = replace(newstr,"[lsfname]" => "noforcing")
            write(fprm,newstr)
        end

        nrun = rundir("IslandDGW",expname,"$(runname).sh")
        open(nrun,"w") do fprm
            newstr = replace(oldstr,"[email]"   => email)
            newstr = replace(newstr,"[dirname]" => projectdir())
            newstr = replace(newstr,"[project]" => "IslandDGW")
            newstr = replace(newstr,"[expname]" => expname)
            newstr = replace(newstr,"[runname]" => runname)
            newstr = replace(newstr,"[sndname]" => case)
            newstr = replace(newstr,"[lsfname]" => "noforcing")
            write(fprm,newstr)
        end
    end
end

orun  = rundir("modifysam","runtemplates","Build_islandwtg.csh")
if doBuild
    open(orun,"r") do frun
        s = read(frun,String)
        for depth in dlist, islandsize in slist
        expname = dampingstrprnt(islandsize)
            
            nrun = rundir("IslandDGW",expname,"Build.csh")
            open(nrun,"w") do wrun
                sn = replace(s ,"[datadir]" => datadir())
                sn = replace(sn,"[project]" => "IslandDGW")
                write(wrun,sn)
            end

        end
    end
end