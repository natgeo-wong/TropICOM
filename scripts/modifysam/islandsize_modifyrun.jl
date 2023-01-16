using DrWatson
@quickactivate "TroPrecLS"
using Printf

orun  = projectdir("run","islandsize_modelrun.sh")

slist = [
    10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100,
    100*sqrt(2),200,200*sqrt(2.5),500,500*sqrt(2),1000
]
dlist = [
    0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
    0.1,0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
    1.,1*sqrt(2),2.,2*sqrt(2.5),5.,5*sqrt(2),
    10.,10*sqrt(2),20.,20*sqrt(2.5),50.,
]

open(orun,"r") do rrun
    oldstr = read(rrun,String)
    for depth in dlist
        depthstr = @sprintf("%05.2f",depth)
        depthstr = replace(depthstr,"."=>"d")
        for islandsize in slist
            sizestr = @sprintf("%04d",islandsize)
            sizestr = sizestr[2:end]

            nrun = projectdir(
                "run","IslandSize","size$(sizestr)km","depth$(depthstr)m.sh"
            )
            open(nrun) do fprm
                newstr = replace(oldstr,"[email]"=>"")
                newstr = replace(newstr,"[project]"=>"TroPrecLS")
                newstr = replace(newstr,"[experiment]"=>"size$(sizestr)km")
                newstr = replace(newstr,"[config]"=>"depth$(depthstr)m")
                newstr = replace(newstr,"[sndname]"=>"control")
                newstr = replace(newstr,"[lsfname]"=>"noforcing")
                write(fprm,newstr)
            end

            mv(tprm,nprm,force=true)

        end
    end
end