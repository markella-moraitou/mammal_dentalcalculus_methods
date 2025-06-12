#!/bin/bash -l

#SBATCH -n 30
#SBATCH --mem 60GB
#SBATCH -t 2:00:00
#SBATCH -J decOM_matrix
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# Create source matrix for decOM

## Change to the directory containing the script
[ -z $SLURM_SUBMIT_DIR ] || cd $SLURM_SUBMIT_DIR

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Project dir: ${PROJDIR}"
echo "Reference dir: ${REFDIR}"

indir=${PROJDIR}/input/ # Directory where the provided input is located
sourcedir=${REFDIR}/decOM_sources/ # Directory to download metagenomes in
outdir=${PROJDIR}/output/ # Output directory

[ -d $REFDIR ] || mkdir $REFDIR # Create directory for all references
[ -d $sourcedir ] || mkdir $sourcedir # Create directory for sources

## Load modules and activate environment
conda activate decOM

cd $sourcedir/sequences

#######################
#### CREATE MATRIX ####
#######################

#echo "Running kmtricks pipeline"
#kmtricks pipeline --file kmtricks.fof --run-dir ../p_sources --mode kmer:pa:bin --restrict-to-list 1

#echo "Running kmtricks aggregate"
#kmtricks aggregate --run-dir ../p_sources --pa-matrix kmer --output ../p_sources/matrices/matrix.pa --format bin

echo "Running kmtricks dump"
kmtricks dump --run-dir ../p_sources --input ../p_sources/matrices/matrix.pa -o ../p_sources/matrices/matrix.pa.txt
