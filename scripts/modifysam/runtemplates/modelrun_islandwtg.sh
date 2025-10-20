#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=0-03:00

#SBATCH --job-name=[sndname]_[expname]
#SBATCH --dependency=singleton

#SBATCH --mail-user=[email]
#SBATCH --mail-type=ALL
#SBATCH --output=./LOGS/samrun.%j.out
#SBATCH --error=./LOGS/samrun.%j.err

module purge
module load openmpi/intel/5.0.8 netcdf-fortran/intel/4.6.2

exproot=[dirname]/exp
prmfile=$exproot/prm/[project]/[expname]/[runname].prm
sndfile=$exproot/snd/[sndname].snd
lsffile=$exproot/lsf/noforcing.lsf

mkdir [sndname]
scp -r CASE_TEMPLATE/* [sndname]/*

prmloc=./[sndname]/prm
sndloc=./[sndname]/snd
lsfloc=./[sndname]/lsf

cp $prmfile $prmloc
cp $sndfile $sndloc
cp $lsffile $lsfloc

scriptdir=$SLURM_SUBMIT_DIR
SAMname=`ls $scriptdir/SAM_*`
echo [sndname] > CaseName

cd ./OUT_3D

for fcom3D in *[expname]*.com3D
do
    rm "$fcom3D"
done

for fcom2D in *[expname]*.com2D
do
    rm "$fcom2D"
done

cd ../OUT_2D

for f2Dcom in *[expname]*.2Dcom
do
    rm "$f2Dcom"
done

cd ../OUT_STAT

for fstat in *[expname]*.stat
do
    rm "$fstat"
done

cd ..

cd $scriptdir
srun $SAMname > ./LOGS/samrun.${SLURM_JOBID}.log

exitstatus=$?
echo SAM stopped with exit status $exitstatus

cd ./OUT_3D

for fcom3D in *[expname]*.com3D
do
    if com3D2nc "$fcom3D" >& /dev/null
    then
        echo "Processing SAM com3D output file $fcom3D ... done"
        rm "$fcom3D"
    else
        echo "Processing SAM com3D output file $fcom3D ... failed"
    fi
done

for fcom2D in *[expname]*.com2D
do
    if com2D2nc "$fcom2D" >& /dev/null
    then
        echo "Processing SAM com2D output file $fcom2D ... done"
        rm "$fcom2D"
    else
        echo "Processing SAM com2D output file $fcom2D ... failed"
    fi
done

cd ../OUT_2D

for f2Dcom in *[expname]*.2Dcom
do
    if 2Dcom2nc "$f2Dcom" >& /dev/null
    then
        echo "Processing SAM 2Dcom output file $f2Dcom ... done"
        rm "$f2Dcom"
    else
        echo "Processing SAM 2Dcom output file $f2Dcom ... failed"
    fi
done

cd ../OUT_STAT

for fstat in *[expname]*.stat
do
    if stat2nc "$fstat" >& /dev/null
    then
        echo "Processing SAM STAT  output file $fstat ... done"
        rm "$fstat"
    else
        echo "Processing SAM STAT  output file $fstat ... failed"
    fi
done

cd ..

exit 0