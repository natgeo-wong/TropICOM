using DrWatson
@quickactivate "TroPrecLS"
using Printf

orun  = projectdir("run","islandsize_modelrun.sh")

slist = [
    5,5*sqrt(2),10,
    10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),
    100,100*sqrt(2),200,200*sqrt(2.5),500,500*sqrt(2),
    1000,1000*sqrt(2),2000,2000*sqrt(2.5),5000.,
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

            nrun = projectdir(
                "run","IslandSize","depth$(depthstr)m","spinup-size$(sizestr)km.sh"
            )
            open(nrun,"w") do fprm
                newstr = replace(oldstr,"[email]"=>"")
                newstr = replace(newstr,"[project]"=>"TroPrecLS")
                newstr = replace(newstr,"[experiment]"=>"size$(sizestr)km")
                newstr = replace(newstr,"[config]"=>"spinup-depth$(depthstr)m")
                newstr = replace(newstr,"[sndname]"=>"control")
                newstr = replace(newstr,"[lsfname]"=>"noforcing")
                write(fprm,newstr)
            end

        end
    end
end