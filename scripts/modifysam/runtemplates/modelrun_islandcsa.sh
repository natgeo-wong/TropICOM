#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=0-09:00

#SBATCH --account=torch_pr_214_general
#SBATCH --job-name=IslandCSA_[expname]_[runname]

#SBATCH --mail-user=[email]
#SBATCH --mail-type=ALL
#SBATCH --output=./LOGS/samrun.%j.out
#SBATCH --error=./LOGS/samrun.%j.err

module unuse /share/apps/modulefiles
module use /share/apps/intel-2023.2.1/modulefiles/

module purge
module load openmpi/intel/5.0.8 netcdf-fortran/intel/4.6.2

exproot=[dirname]/exp
prmfile=$exproot/prm/[project]/[expname]/[runname]/[config1]-[config2].prm
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

cd $scriptdir
mpirun -np $SLURM_NTASKS $SAMname > ./LOGS/samrun.${SLURM_JOBID}.log

exitstatus=$?
echo SAM stopped with exit status $exitstatus

cd ./OUT_3D

for fcom3D in *[config1]-[config2]*.com3D
do
    if com3D2nc "$fcom3D" >& /dev/null
    then
        echo "Processing SAM com3D output file $fcom3D ... done"
        rm "$fcom3D"
    else
        echo "Processing SAM com3D output file $fcom3D ... failed"
    fi
done

mkdir [config1]-[config2]
for fnc3D in *[config1]-[config2]*.nc
do
    mv $fnc3D [config1]-[config2]/
done

cd ../OUT_2D

for f2Dcom in *[config1]-[config2]*.2Dcom
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

for fstat in *[config1]-[config2]*.stat
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