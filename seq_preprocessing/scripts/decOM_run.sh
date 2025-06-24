#!/bin/bash -l

#SBATCH -n 20
#SBATCH --mem 256GB
#SBATCH -t 20:00:00
#SBATCH -J decOM
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script run decOM using a custom source matrix

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
scriptdir=${PROJDIR}/scripts/ # Directory where the scripts are saved
outdir=${PROJDIR}/output/ # Output directory

sourcedir=${REFDIR}/decOM_sources/p_sources # Sources directory
sinkdir=${outdir}/2_mapped_seqs/ # Sinks directory

## Instructions to install decOM - More in https://github.com/CamilaDuitama/decOM.git
## Installing decOM
#mkdir ${PROJDIR}/software

#cd ${PROJDIR}/software
#git clone https://github.com/CamilaDuitama/decOM.git
#cd decOM
#conda env create -n decOM --file environment.yml
#conda deactivate

# Downloading sources
#wget https://zenodo.org/record/6513520/files/decOM_sources.tar.gz
#tar -xf decOM_sources.tar.gz

export PATH=${PROJDIR}/software/decOM:${PATH}

## Determine samples to be run
sample_list=${indir}/file_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate decOM

## Print some info
echo "Sinks: ${sinkdir}"
echo "Sources: ${sourcedir}"

echo "Sample list:"
echo ${sample_list}
cat $sample_list

echo "Global modules:"
module list
echo "Conda modules:"
conda list

if [ -z "$SLURM_NTASKS" ]; then
  n_proc=1
else
  n_proc=$SLURM_NTASKS
fi
echo "Number of tasks: ${n_proc}"

###################
#### RUN DECOM ####
###################

cd $sinkdir

decom_in=./decOM_input

[ ! -d $sourcedir ] && echo "Sources not found! Download them and provide path." && exit

[ -d $decom_in ] || ( mkdir $decom_in && echo "Creating directory $decom_in" ) # Create subdirectory for the input files

# Get sinks.txt file to use with -p_sinks as well as keys files to use with -p_keys, using the sample list

[ -f $decom_in/sinks.txt ] && rm -f $decom_in/sinks.txt && echo "Removing existing sinks.txt file"

ls *unmapped.fastq.gz | while read input
do
  s=${input%_unmapped.fastq.gz}
  if [[ ! -f $input ]]
  then
    echo "File $input doesn't exist. Skipping"
  else
    echo "Creating input files for $s"
    echo $s >> $decom_in/sinks.txt
    echo "${s} : ${input}" > $decom_in/${s}.fof
  fi
done

# Prep map: decOM_sources.csv doesn't include separate lines for forward and reverse reads,
# so we need to get these from kmtricks.fof

echo "SampleID,Env" > $decom_in/map.csv # Header for the mapping file

cat $indir/decOM_sources.csv | while IFS=, read -r sample env _ _
do
  # find sample in kmtricks.fof and extract id
  grep $sample $sourcedir/kmtricks.fof | while IFS=" : " read -r sample_new file
  do
    echo "${sample_new},${env}" >> $decom_in/map.csv
  done
done

sinks=$decom_in/sinks.txt

[ -d decOM_output ] && rm -r ./decOM_output

decOM-MST -p_sources $sourcedir -m $decom_in/map.csv -p_sinks $sinks -p_keys $decom_in -mem 200GB -t 1 -o ./decOM_output
