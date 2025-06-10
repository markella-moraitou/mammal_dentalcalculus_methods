#!/bin/bash -l

#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH -J trimming
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# Cleaning up and merging reads (both this study and previously available)
# Trims poly-G and poly-X tails (resulting from NovaSeq two-colour chemistry), low-quality bases and Ns (sliding windows of 3)
# Removes adapters and barcodes
# Also, filters out sequences less than 30 bp

################
#### SET UP ####
################

## Change to the directory containing the script
[ -z $SLURM_SUBMIT_DIR ] || cd $SLURM_SUBMIT_DIR

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Raw data dir: ${RAWDIR}"
echo "Project dir: ${PROJDIR}"
echo "Reference dir: ${REFDIR}"

indir=${PROJDIR}/input/ # Directory where the provided input is located
scriptdir=${PROJDIR}/scripts/ # Directory where the scripts are saved
outdir=${PROJDIR}/output/ # Output directory
subdir=${outdir}/1_submitted_seqs
demultdir=${outdir}/0_demultiplexed_seqs/

[ -d $subdir ] || mkdir $subdir # Create subdir if not there
[ -d ${scriptdir}/logs ] || mkdir ${scriptdir}/logs # Create directory for script logs if not there

## Determine samples to be run
sample_list=${indir}/file_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info

echo "Saving fastq output to: ${subdir}"
echo "Saving reports to: ${repdir}"

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
cd $subdir/

#################
#### PROCESS ####
#################

#### Samples from this study ####

# Process samples in provided samples list - First filter for the newly generated data
tail -n+2 ${sample_list} | while IFS=, read -r sample species fwd_path rev_path
do
  fwd_file=$demultdir/${sample}_run1_F.fastq.gz
  rev_file=$demultdir/${sample}_run1_R.fastq.gz
  if [[ ! -f $fwd_file || ! -f $rev_file ]]; then
    echo "Input file(s) missing for $sample"
  fi
  
  echo "Processing $sample"
  # Process with fastq
  # --trim_poly_g: removes low quality reads (both sides, mean quality 30, window size 3)and poly-G tails (due to NovaSeq 2-colour chemistry)
  # -3 --cut_tail_window_size 3 --cut_tail_mean_quality 30: trim seqs with qual < 30 in sliding windows from 3'
  # --adapter_fasta: specify adapter sequences (Ns not allowed, so I need to specify one for each index)
  # --trim_front1 7 and --trim_front2 7: remove barcodes by trimming the first 7 nucleotides from both forward and reverse
  # --merge, --overlap_len_require 11 and --overlap_diff_limit 3: merge F and R if there is at least 11 of overlap with at most 3 mismatches (probability of this happening at random <1%)
  # --length_required 30: filters out reads shorter than 30
  # Saving merged reads in one file
  # and unpaired F and R in separate files (regardless if they could not be merged or if the other mate failed the filters)
  fastp \
    --in1 $fwd_file \
    --in2 $rev_file \
    --out1 ${subdir}/${sample}_trimmed_F.fastq.gz \
    --out2 ${subdir}/${sample}_trimmed_R.fastq.gz \
    --trim_poly_g \
    --trim_poly_x \
    -3 --cut_tail_window_size 3 --cut_tail_mean_quality 30 \
    --adapter_fasta ${indir}/adapter_seqs.fasta \
    --trim_front1 7 \
    --trim_front2 7 \
    --json ${sample}.fastp.json  \
    --html ${sample}.fastp.html  \
    --thread $n_proc
done

##################################
#### RUN QC AND SAVE METADATA ####
##################################

## The following bit generates a MultiQC report and then extracts info from it to add to read_count.csv and read_length.csv
## These to files will be updated after every step, to keep track of the pre-processing pipeline

# Since I have to run it twice, first to extract columns "FastQC_mqc-generalstats-fastqc-total_sequences" from a few different MultiQC output tables and add them to read_count.csv and
# then to extract columns "FastQC_mqc-generalstats-fastqc-avg_sequence_length" (again from a few tables) and add them to read_length.tsv,
# I will define the following actions in a function that I can call

## Get QC reports for each type of output

# Run FastQC
fastqc *_trimmed*.fastq.gz -f fastq -o . -t $n_proc

# Generate MultiQC report
multiqc *.zip -f --interactive --title "trimmed"
rm -f *fastqc*

## Get average read length and count per sample (both runs combined)
# Average length and count tables combined for the two runs

echo "sample,submitted_F,submitted_R" > ${outdir}/read_length.csv
echo "sample,submitted_F,submitted_R" > ${outdir}/read_count.csv

ls *_trimmed_F.fastq.gz | while read fwd_file
do
  fwd=$(basename ${fwd_file%.fastq.gz})
  sample=${fwd%_trimmed_F}
  rev=${sample}_trimmed_R
  
  # Get number of reads from multiqc_data
  rc_F=$(grep "^$fwd" trimmed_multiqc_report_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $5}')
  rc_R=$(grep "^$rev" trimmed_multiqc_report_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $5}')
  
  rl_F=$(grep "^$fwd" trimmed_multiqc_report_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $11}')
  rl_R=$(grep "^$rev" trimmed_multiqc_report_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $11}')

  # Add to file
  echo "$sample,$rc_F,$rc_R" >> ${outdir}/read_count.csv
  echo "$sample,$rl_F,$rl_R" >> ${outdir}/read_length.csv
done

