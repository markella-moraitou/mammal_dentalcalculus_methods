#!/bin/bash -l

#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J mapping_refs
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script prepares the reference for mapping
# A large combined reference of all host genomes plus the human and phiX genome
# This allows us not only to remove all host associated reads, but also, with downstream analysis, to identify potentially mislabelled samples

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

host_refdir=${REFDIR}/host_references/ # Directory where the reference genomes are stored
conc_refdir=${REFDIR}/concatenated_references/ # Directory to store concatenated references (host-human)

[ -d $conc_refdir ] || mkdir $conc_refdir # Create directory for concatenated references

human_gen=$( find ${REFDIR}/Homo_sapiens_genome/ncbi_dataset/data/ -name "*.gz" )
phiX_gen=$( find ${REFDIR}/phiX174_genome/ncbi_dataset/data/ -name "*.gz" )

# Check if the human and PhiX genomes exists, if not exit
if [ ! -f $human_gen ]
then
  echo "Human genome not found! Exiting..."
  exit
fi

if [ ! -f $phiX_gen ]
then
  echo "PhiX genome not found! Exiting..."
  exit
fi

## Determine samples to be run
sample_list=${indir}/file_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info
echo "Taking input from: ${host_refdir}"
echo "Saving concatenated references to: ${conc_refdir}"

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

###################################
#### PREPARE MAPPING REFERENCE ####
###################################

# This requires the host_genome_per_species.csv file
# Check if it exists, if not exit
if [ ! -f ${outdir}/host_genome_per_species.csv ]
then
  echo "host_genome_per_species.csv not found! Exiting..."
  exit
fi

## Identify the appropriate host reference genome for each species and create a concatenated version with human and phix, to be used for competitive mapping
echo "species,reference_genome" > ${outdir}/concatenated_genome_per_species.csv

# Iterate over the unique species in the sample list
awk -F "," 'NR>1 && $2!="" {print $2}' ${sample_list} | sort | uniq | while read species
do
  species=${species// /"_"} # Replace whitespace with underscore
  # If the "species" column indicates a blank or control, there is no host reference -- concatenate PhiX and human only
  if [[ $species == "Extraction_blank" ]] || [[ $species == "Library_blank" ]] || [[ $species == "Environmental_control" ]]
  then
    full_ref=${conc_refdir}/concat_human_phix_genomes.fa.gz 
    if [[ ! -f ${full_ref} ]]
    then
      # Concatenate the host reference with the human reference
      echo "Creating concatenated reference for Homo sapiens and PhiX"
      cat $human_gen $phiX_gen> ${full_ref}
    fi
  # For all true samples, the mapping will be done against a concatenated reference of the host and the human genome
  else
    # Find the alternate species being used if there was not ref genome for this species
    alt_species=$( awk -F "," -v species=$species '$1 == species {print $2}' ${outdir}/host_genome_per_species.csv )
    # Find the reference genome using the host_genome_per_species.csv file
    refgen=$(grep ${alt_species}_genome ${outdir}/host_genome_per_species.csv | awk -F "," '{print $3}' )
    # Print a warning if the reference genome was not found in the file
    if [[ -z $refgen ]]
    then
      echo "Genome for ${species} not found in host_genome_per_species.csv"
    fi
    
    full_ref=${conc_refdir}/concat_${alt_species// /"_"}_human_phix_genomes.fa.gz # Name for the concatenated reference
    # Check it it already exists, if not create it
    if [[ ! -f ${full_ref} ]]
    then 
      # Concatenate the host reference with the human reference
      echo "Creating concatenated reference for ${species}, Homo sapiens and PhiX: ${full_ref}"
      cat ${refgen} $human_gen $phiX_gen > ${full_ref}
    fi
  fi
  
  # If index is missing, create it
  if [ ! -f ${full_ref%.gz}.pac ] || [ ! -f ${full_ref%.gz}.amb ] || [ ! -f ${full_ref%.gz}.ann ] || [ ! -f ${full_ref%.gz}.bwt ] || [ ! -f ${full_ref%.gz}.sa ]
  then
    echo "Indexing ${full_ref}"
    # Calculate batch size - we can scale it according to the available memory, the default being 10^6
    batch_size=$( echo "(10^6) * 200 * ${n_proc}" | bc )
    # Decompress reference to index, then recompress
    pigz -d -p ${n_proc} $full_ref
    bwa index -b $batch_size -a bwtsw ${full_ref%.gz}
    pigz -p ${n_proc} ${full_ref%.gz}
  fi
  
  # Print out a file with the full reference per input file
  echo "${species},${full_ref}" >> ${outdir}/concatenated_genome_per_species.csv
done

# First, identify all contigs and which genome they come from
echo "contig,reference" > ${outdir}/contigs_per_reference.csv
cat ${outdir}/host_genome_per_species.csv | while IFS=, read -r orig_species alt_species ref
do
  zcat $ref | awk -v species="$alt_species" -v OFS="," -F " " '/^>/ {print $1, species}' >> ${outdir}/contigs_per_reference.csv
done

# Add human and phix too
zcat $human_gen | awk -v OFS="," -F " " '/^>/ {print $1, "Homo sapiens"}' >> ${outdir}/contigs_per_reference.csv
zcat $phiX_gen | awk -v OFS="," -F " " '/^>/ {print $1, "PhiX"}' >> ${outdir}/contigs_per_reference.csv
