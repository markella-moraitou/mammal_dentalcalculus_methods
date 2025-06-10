#!/bin/bash -l

#SBATCH -A snic2022-5-561
#SBATCH -n 12
#SBATCH -t 10:00:00
#SBATCH -J deduplication
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script uses the dedupe module from the BBTools package to remove PCR duplicates

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
seqdir=${outdir}/1_merged_seqs
subdir=${outdir}/2_deduplicated_seqs/ # Output subdirectory for the deduplicated reads reads

[ -d $subdir ] || mkdir $subdir # Create subdirectory if not there

## Determine samples to be run
sample_list=${indir}/file_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info

echo "Saving mapped fastq to: ${subdir}"

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

#################
#### PROCESS ####
#################

cd $subdir 

tail -n+2 ${sample_list} | while IFS=, read -r sample species _ _
do
  # Find input file in the fastq output ($outdir/1_cleaned_seqs)
  input=${seqdir}/${sample}_merged.fastq.gz
  # Check if input files exists if not, print a message and continue
  if [[ ! -f $input ]]
  then
    echo "${input} does not exist"
    continue
  fi
  output=${sample}_deduped.fastq.gz
  
  # Deduplicate
  echo "Deduplicating $sample"
  dedupe.sh -Xmx250g in=$input out=$output ac=false
done

##################################
#### RUN QC AND SAVE METADATA ####
##################################

# Run FastQC
fastqc *_deduped.fastq.gz -f fastq -o . -t $SLURM_NTASKS 

# Generate MultiQC report
multiqc *.zip -f --interactive --title deduped

# Remove FastQC files
rm -f *fastqc*

# Get correct sample names
sed "s/_deduped//g" deduped_multiqc_report_data/multiqc_fastqc.txt > deduped_general_stats.tmp

# Run python script
python $scriptdir/modules/multiqc_to_csv.py deduped_general_stats.tmp ${outdir}/read_count.csv "Total Sequences" "deduped"
python $scriptdir/modules/multiqc_to_csv.py deduped_general_stats.tmp ${outdir}/read_length.csv "avg_sequence_length" "deduped"

rm *tmp
