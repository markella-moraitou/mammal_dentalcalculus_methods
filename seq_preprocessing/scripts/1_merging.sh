#!/bin/bash -l

#SBATCH -n 12
#SBATCH -t 3:00:00
#SBATCH -J merging
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# Merging reads
# Also, filters out sequences less than 30 bp

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
seqdir=${outdir}/0_submitted_seqs
subdir=${outdir}/1_merged_seqs

[ -d $subdir ] || mkdir $subdir # Create subdir if not there
[ -d ${scriptdir}/logs ] || mkdir ${scriptdir}/logs # Create directory for script logs if not there

## Determine samples to be run
sample_list=${indir}/file_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info

echo "Saving fastq output to: ${subdir}"

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
  fwd_file=${seqdir}/${sample}_trimmed_F.fastq.gz
  rev_file=${seqdir}/${sample}_trimmed_R.fastq.gz
  if [[ ! -f $fwd_file || ! -f $rev_file ]]; then
    echo "Input file(s) missing for $sample"
  fi
  
  echo "Processing $sample"
  # --merge, --overlap_len_require 11 and --overlap_diff_limit 3: merge F and R if there is at least 11 of overlap with at most 3 mismatches (probability of this happening at random <1%)
  # --length_required 30: filters out reads shorter than 30
  # Saving merged reads in one file
  # and unpaired F and R in separate files (regardless if they could not be merged or if the other mate failed the filters)
  fastp \
    --in1 $fwd_file \
    --in2 $rev_file \
    --merge \
    --merged_out ${subdir}/${sample}_merged.fastq.gz \
    --unpaired1 ${subdir}/${sample}_unpaired_F.fastq.gz \
    --out1 ${subdir}/${sample}_unmerged_F.fastq.gz \
    --unpaired2 ${subdir}/${sample}_unpaired_R.fastq.gz \
    --out2 ${subdir}/${sample}_unmerged_R.fastq.gz \
    --overlap_len_require 11 \
    --overlap_diff_limit 3 \
    --dedup \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --disable_trim_poly_g \
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
fastqc *_merged*.fastq.gz -f fastq -o . -t $n_proc

# Generate MultiQC report
multiqc *.zip -f --interactive --title "merged"
rm -f *fastqc*

# Get correct sample names
sed "s/_merged//g" merged_multiqc_report_data/multiqc_fastqc.txt > merged_general_stats.tmp

# Run python script
python $scriptdir/modules/multiqc_to_csv.py merged_general_stats.tmp ${outdir}/read_count.csv "Total Sequences" "merged"
python $scriptdir/modules/multiqc_to_csv.py merged_general_stats.tmp ${outdir}/read_length.csv "avg_sequence_length" "merged"

rm *tmp
