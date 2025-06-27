##### CREATE INPUT FOR STARS #####

#### Create input tables from seq preprocessing output

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(stringr)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory paths
indir <- normalizePath(file.path("..", "..", "seq_preprocessing","output"))
outdir <- normalizePath(file.path("..","input")) # output directory

#####  Load input data #####
read_count <- read.csv(file.path(indir, "read_count.csv"))
read_length <- read.csv(file.path(indir, "read_length.csv"))
fastqc <- read.csv(file.path(indir, "0_submitted_seqs", "trimmed_multiqc_report_data", "multiqc_fastqc.txt"), sep = "\t")
host_human <- read.csv(file.path(indir, "2_mapped_seqs", "genome_mapping_stats.tsv"), sep = "\t")
decom <- read.csv(file.path(indir, "2_mapped_seqs", "decOM_output", "decOM_output.csv"))

##### Combine tables #####

# Rename column
colnames(read_count)[2:ncol(read_count)] <- paste0(colnames(read_count)[2:ncol(read_count)], "_count")
colnames(read_count)[1] <- "Ext.ID"

colnames(read_length)[2:ncol(read_length)] <- paste0(colnames(read_length)[2:ncol(read_length)], "_avlength")
colnames(read_length)[1] <- "Ext.ID"

colnames(host_human) <- c("Ext.ID", "host_count", "human_count", "phix_count")

# Filter and rename
decom <- decom %>% select(Sink, starts_with("p_"))
colnames(decom)[which(colnames(decom) == "Sink")] <- "Ext.ID"
decom$Ext.ID <- str_remove(decom$Ext.ID, "_$")

# Get percentage of unique reads per sample
fastqc <- fastqc %>%
        select(Sample, total_deduplicated_percentage) %>% rename(Ext.ID = Sample) %>%
        mutate(Ext.ID = str_remove(Ext.ID, "_trimmed_.*")) %>% group_by(Ext.ID) %>%
        summarise(perc_unique = mean(total_deduplicated_percentage))

##### Combine & Save #####

bioinfo_data <- full_join(read_count, host_human, by = "Ext.ID") %>%
                full_join(read_length, by = "Ext.ID") %>%
                full_join(fastqc, by = "Ext.ID") %>%
                full_join(decom, by = "Ext.ID") %>%
                mutate(submitted_count = submitted_F_count) %>%
                mutate(submitted_F_count = NULL,
                       submitted_R_count = NULL)

write.csv(bioinfo_data, file = file.path(outdir, "bioinfo_metadata.csv"), quote = FALSE, row.names = FALSE)
