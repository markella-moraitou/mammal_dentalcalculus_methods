#!/bin/bash -l

#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J ncbi_genome_download
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# Get reference genomes per species
# data downloaded in 3 steps
# Download dehydrated zip archive

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
host_refdir=${REFDIR}/host_references/ # Directory where the host reference genomes are stored
outdir=${PROJDIR}/output/ # Output directory

[ -d $REFDIR ] || mkdir $REFDIR # Create directory for all references
[ -d $host_refdir ] || mkdir $host_refdir # Create directory for concatenated references

## Load modules and activate environment
conda activate oral_mb_evol

###########################
##### DOWNLOAD GENOMES ####
###########################

cd $host_refdir

# Download reference genomes using accession number - For species that didn't have a published genome, the genome of an alternate related species is downloaded
{ read
while IFS=',' read -r -a fields
do
    species=${fields[0]// /_}
    taxid=${fields[1]}
    acc=${fields[2]}
    alt_species=${fields[3]// /_}
    for line1 in "${species[@]}" 
    do
    for line2 in "${acc[@]}" # use accession number to find genome
    do
    for line3 in "${alt_species[@]}" # use accession number to find genome
    do
        datasets download genome accession "$line2" --api-key $API_KEY --dehydrated --filename "${line3}_genome.zip" # only reference genomes
    done
    done
    done
done } < ${indir}/ncbi_ref_genomes.csv

# unzip

for f in ./*.zip; do
    # Creating a name for the new directory.
    new_dir="${f##*/}"
    new_dir="${new_dir%.*}"
    # Creating the directory if it doesn't already exist.
    mkdir -p "$new_dir"
    # Unzip contents of "$f" to "$new_dir"
    unzip -d "$new_dir" -- "$f"
done

rm *.zip

# rehydrate to retrieve data

ls -d *_genome/ | while read i
do
  datasets rehydrate --api-key $API_KEY --gzip --directory "$i"
done

# Also download the human and PhiX genomes -- these will be used for competitive mapping

cd $REFDIR

datasets download genome accession "GCA_000001405.29" --api-key $API_KEY --dehydrated --filename "Homo_sapiens_genome.zip"
unzip -d "Homo_sapiens_genome" -- "Homo_sapiens_genome.zip"
datasets rehydrate --api-key $API_KEY --gzip --directory "Homo_sapiens_genome"

datasets download genome accession "GCA_000819615.1" --api-key $API_KEY --dehydrated --filename "phiX174_genome.zip"
unzip -d "phiX174_genome" -- "phiX174_genome.zip"
datasets rehydrate --api-key $API_KEY --gzip --directory "phiX174_genome"

rm *zip

#######################
##### OUTPUT TABLE ####
#######################

# Get a file indicating the reference genome for each species

[ -f ${outdir}/host_genome_per_species.csv ] && rm ${outdir}/host_genome_per_species.csv

{ read
while IFS=, read -r spe _ _ alt
do
  spe=${spe// /_}
  alt=${alt// /_}
  # Find the genome file using the alt column
  genome=$(find "${host_refdir}/${alt}_genome/ncbi_dataset/data/" -type f -name "*_genomic.fna.gz" 2>/dev/null)
  # If no genome found, continue
  if [[ -z $genome ]]
  then
    continue
  # If more than one genome was found, select the first one and print a warning
  elif [[ $(echo "${genome}" | wc -l ) -gt 1 ]]
  then
    genome=$(echo "${genome}" | head -1 )
    echo "More than one reference genome found for ${spe}. Keeping only the first one."
  fi
  echo "${spe},${alt},${genome}" >> ${outdir}/host_genome_per_species.csv
done } < ${indir}/ncbi_ref_genomes.csv 
