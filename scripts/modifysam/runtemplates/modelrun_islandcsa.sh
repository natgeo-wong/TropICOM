#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=0-06:00

#SBATCH --job-name=IslandCSA

#SBATCH --mail-user=[email]
#SBATCH --mail-type=ALL
#SBATCH --output=./LOGS/samrun.%j.out
#SBATCH --error=./LOGS/samrun.%j.err

module purge
module load openmpi/intel/5.0.8 netcdf-fortran/intel/4.6.2

exproot=[dirname]/exp
prmfile=$exproot/prm/[project]/[runname].prm
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
srun $SAMname > ./LOGS/samrun.${SLURM_JOBID}.log

exitstatus=$?
echo SAM stopped with exit status $exitstatus

exit 0