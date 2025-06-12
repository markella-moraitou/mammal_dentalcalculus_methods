#!/bin/bash -l

#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J decOM_sources
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# Download metagenomes to be used as sources for decOM

################
#### SET UP ####
################

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
conda activate oral_mb_evol

# Install ENA Browser tool
#cd $PROJDIR/software
#git clone https://github.com/enasequence/enaBrowserTools.git

alias enaDataGet=$PROJDIR/software/enaBrowserTools/python3/enaDataGet

cd $sourcedir

###############################
##### DOWNLOAD METAGENOMES ####
###############################

[ -d sequences ] || mkdir sequences
cd sequences

# Parse decOM sources and download files
while IFS=',' read -r acc source; do
    # Skip header line
    [[ "$acc" == "SampleID" || "$source" == "Env" ]] && continue
    # Check if any files with the accession number already exist
    if find . -type f -name "*${acc}*" | grep -q .; then
        echo "Files for $acc already exist. Skipping download."
    else
        # Download data
        enaDataGet -f fastq -d . "$acc"
    fi
    # Create kmtricks.fof
    find . -name "${acc}*.gz" | while read file;
    do
        file=$(basename $file)
        echo "${file%.fastq.gz} : $file" >> kmtricks.fof
    done
done <  "${indir}/decOM_sources.csv"

# Move all fastq files to parent dir and remove empty subdirs
find . -type f -name "*.fastq.gz" -exec mv {} . \;
find . -type d -empty -delete
