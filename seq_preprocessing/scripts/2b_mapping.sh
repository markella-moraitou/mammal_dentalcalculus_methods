#!/bin/bash -l

#SBATCH -n 20
#SBATCH -t 5-00:00:00
#SBATCH -J mapping
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script takes adapter and barcode trimmed merged reads and
# maps against the human and host genomes to maintain only unmapped reads

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

conc_refdir=${REFDIR}/concatenated_references/ # Directory to store concatenated references (host-human)
seqdir=${outdir}/1_merged_seqs/ # Directory where the sequences are stored
subdir=${outdir}/2_mapped_seqs/ # Output subdirectory for the mapping reads

[ -d $subdir ] || mkdir $subdir # Create subdirectory if not there

## Determine samples to be run
sample_list=${indir}/file_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info
echo "Taking input from: ${seqdir}"
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
#### MAPPING ####
#################

# This requires the host_genome_per_species.csv file
# Check if it exists, if not exit
if [ ! -f ${outdir}/concatenated_genome_per_species.csv ]
then
  echo "concatenated_genome_per_species.csv not found! Exiting..."
  exit
fi

cd $subdir

echo "Aligning based on ${outdir}/concatenated_genome_per_species.csv"
cat ${outdir}/concatenated_genome_per_species.csv

tail -n+2 ${sample_list} | while IFS=, read -r sample species _ _
do
  species=${species// /"_"} # Replace whitespace with underscore
  
  # Find input file in the fastq output ($outdir/1_cleaned_seqs)
  input=${seqdir}/${sample}_merged.fastq.gz
  # Check if input files exists if not, print a message and continue
  if [[ ! -f $input ]]
  then
    echo "${input} does not exist"
    continue
  fi
  
  # Check if output file exists already, if yes skip
  if [[ ! -f ${sample}_mapped.fastq  ]] && [[ -f ${sample}_mapped.fastq.gz ]]
  then
    echo "Output for ${sample} exists, skipping"
    continue
  fi
  
  # Get reference from ${outdir}/concatenated_genome_per_species.csv 
  ref=$(awk -F "," -v species=$species '$1==species {print $2}' ${outdir}/concatenated_genome_per_species.csv)
  echo "Aligning $sample to $(basename ${ref%.gz} )"
  # Align with bwa aln -- using the parameters suggested in the nfcore-eager pipeline
  # Look for unzipped files, as the reference genomes were indexed before they were compressed
  bwa aln -n 0.04 -k 2 -l 1024 -o 2 -t $n_proc ${ref%.gz} $input > ${sample}_alignment.sai
  # Convert .sai to .sam and then to a sorted .bam
  bwa samse -r "@RG\\tID:${sample}\\tSM:${sample}\\tPL:illumina" ${ref%.gz} ${sample}_alignment.sai ${input} | samtools sort - -@ $n_proc -O bam > ${sample}.bam
  # Generate .csi index
  samtools index -c ${sample}.bam
  # Export unmapped reads as a fastq file
  samtools view -Sb -f 4 -@ 10 ${sample}.bam > ${sample}_unmapped.bam
  bedtools bamtofastq -i ${sample}_unmapped.bam -fq ${sample}_unmapped.fastq
  # Do the same for mapped reads
  samtools view -Sb -F 4 -@ 10 ${sample}.bam > ${sample}_mapped.bam
  bedtools bamtofastq -i ${sample}_mapped.bam -fq ${sample}_mapped.fastq
  # Compressed output
  pigz -p $n_proc ${sample}_unmapped.fastq
  pigz -p $n_proc ${sample}_mapped.fastq
done

##################################
#### RUN QC AND SAVE METADATA ####
##################################

## For unmapped reads ##
# Run FastQC
fastqc *_unmapped.fastq.gz -f fastq -o . -t $n_proc 

# Generate MultiQC report
multiqc *.zip -f --interactive --title unmapped

# Remove FastQC files
rm -f *fastqc*

# Get correct sample names
sed "s/_unmapped//g" unmapped_multiqc_report_data/multiqc_fastqc.txt > unmapped_general_stats.tmp

# Run python script
python $scriptdir/modules/multiqc_to_csv.py unmapped_general_stats.tmp ${outdir}/read_count.csv "Total Sequences" "unmapped"
python $scriptdir/modules/multiqc_to_csv.py unmapped_general_stats.tmp ${outdir}/read_length.csv "avg_sequence_length" "unmapped"

## For mapped reads ##
# Run FastQC
fastqc *_mapped.fastq.gz -f fastq -o . -t $n_proc 

# Generate MultiQC report
multiqc *.zip -f --interactive --title mapped

# Remove FastQC files
rm -f *fastqc*

# Get correct sample names
sed "s/_mapped//g" mapped_multiqc_report_data/multiqc_fastqc.txt > mapped_general_stats.tmp

# Run python script
python $scriptdir/modules/multiqc_to_csv.py mapped_general_stats.tmp ${outdir}/read_count.csv "Total Sequences" "mapped"
python $scriptdir/modules/multiqc_to_csv.py mapped_general_stats.tmp ${outdir}/read_length.csv "avg_sequence_length" "mapped"

rm *tmp

## Get stats on host-mapped, human-mapped and PhiX-mapped

# First get stats per chromosome
ls ${subdir}*_mapped.bam | while read bam
do
  sample=$(basename ${bam%.bam})
  samtools idxstats $bam > ${sample}_idxstats.txt
done

# Summarise mapping stats per genome using in-house script (generates genome_mapping_stats.tsv output in the same folder)
python ${scriptdir}/modules/genome_mapping_stats.py "_idxstats.txt" ${outdir}/contigs_per_reference.csv
