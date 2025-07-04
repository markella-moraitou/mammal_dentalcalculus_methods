##### SEQUENCING OUTPUT #####

#### Explore sequencing output

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(scales)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory paths
indir <- normalizePath(file.path("..","input"))
lmdir <- normalizePath(file.path("..","output", "LM"))
outdir <- normalizePath(file.path("..","output", "SEQ")) # output directory

# Create directory for output
dir.create(outdir, recursive  =  TRUE, showWarnings  =  FALSE)

# Theme
source("ggplot_theme.R")
theme_set(theme_bw_alt)

#######################
#####  LOAD INPUT #####
#######################

# Filtered data
filt_dat <- read.table(file.path(lmdir, "filtered_data.tsv"), sep = "\t", quote = "", header = TRUE)

# Load palettes
for (file in dir(file.path(indir, "palettes"))) {
  name <- str_remove(file, ".csv")
  palette <- read.csv(file.path(indir, "palettes", file)) %>% pull(2)
  names(palette) <-  read.csv(file.path(indir, "palettes", file)) %>% pull(1)
  assign(name, palette)
}

# Pool sequencing info
pool_info <- data.frame(name=c("Pool A", "Pool B", "Pool C", "Pool D"),
                        pool_vol=c(264, 327, 484, 387),
                        pool_read_depth=c(1749653386, 2024338024, 1545135146, 1373545594))

#######################
####      PLOT     ####
#######################

# Do samples differ in the number of number of reads they generate?
#indexing_final_vol_ul = 20

seq_data <- filt_dat %>%
  filter(!is.na(submitted_count) & submitted_count!=0) %>%
  select(-c(Collection.notes, Country, Region, Locality, Sampling.comments, Appearance)) %>%
  mutate(failed_indexing = (ind_output_copies < 10^9)) %>%
  # Add info about pool volume and seq output
  full_join(pool_info, by=c("Pool"="name")) %>% 
  # calculate copies added to pool per species copies(t) = copies(0) * vol(t) / vol(0)
  mutate(pooled_copies=ind_output_copies*(Pooled.volume.ul./20)) %>%
  # Then calculate total copies in pool
  group_by(Pool) %>% mutate(pooled_copies_total = sum(pooled_copies)) %>% ungroup %>%
  mutate(pool_percentage=pooled_copies/pooled_copies_total) %>%
  # Calculate the expected number of reads
  mutate(expected_reads=pool_read_depth*pool_percentage) %>%
  # Calculate number of unique reads
  mutate(unique_count = submitted_count * perc_unique/100)

## Is there a correlation between the percentage of unique reads and the amount of input DNA?
# Percentage of unique reads is also affected by sequencing depth,
# for this reason we will limit this question to data with similar read depths

binned_seq_data <- seq_data %>%
  # Remove outliers
  filter(submitted_count < 2e+07 & submitted_count > 1e+06) %>%
  mutate(quartile = ntile(submitted_count, 4)) %>%
  group_by(quartile) %>%
  mutate(label = paste0(
    "Quartile ", quartile, "\n(", 
    format(min(submitted_count), big.mark = ",", scientific = FALSE), " - ", 
    format(max(submitted_count), big.mark = ",", scientific = FALSE), 
    ")"
  ))

input_vs_unique <-
  ggplot(data = binned_seq_data, aes(x = DNA_input_ug, y = perc_unique, colour = submitted_count)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_continuous(trans = "log10", name = "DNA input (μg)")  +
  scale_color_viridis_c(trans = "log10", breaks = c(2e+6, 4e+6, 8e+6, 16e+6), name = "Total sequencing reads",
                        guide = guide_colourbar(position = "top",
                                                theme = theme(
                                                  legend.key.width  = unit(20, "lines"),
                                                  legend.key.height = unit(1, "lines"),
                                                  legend.text = element_text(angle = 90)
                                                ))) +
  ylim(c(50,100)) + ylab("% of unique sequencing reads") +
  ylim(c(50,100)) + ylab("% of unique sequencing reads") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), legend.position = "top")

input_vs_unique

ggsave(input_vs_unique, file = file.path(outdir, "input_vs_unique.png"),
       width = 12, height = 8)

input_vs_unique_bin <-
  ggplot(data = binned_seq_data, aes(x = DNA_input_ug, y = perc_unique, colour = submitted_count)) +
  geom_point(size = 3, alpha = 0.7) +
  facet_grid(rows = vars(label)) +
  scale_x_continuous(trans = "log10", name = "DNA input (μg)") +
  scale_color_viridis_c(trans = "log10", breaks = c(2e+6, 4e+6, 8e+6, 16e+6), name = "Total sequencing reads",
                        guide = guide_colourbar(position = "top",
                                                theme = theme(
                                                  legend.key.width  = unit(20, "lines"),
                                                  legend.key.height = unit(1, "lines"),
                                                  legend.text = element_text(angle = 90)
                                                  ))) +
  ylim(c(50,100)) + ylab("% of unique sequencing reads") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

input_vs_unique_bin

ggsave(input_vs_unique_bin, file = file.path(outdir, "input_vs_unique_binned.png"),
       width = 12, height = 16)

# Above 0.01 ug of DNA how many samples have at least 90% unique reads
binned_seq_data %>% ungroup %>% mutate(high_DNA = (DNA_input_ug > 0.01),
                           complex_library = (perc_unique > 90)) %>%
  dplyr::select(high_DNA, complex_library) %>% table()

# Save table
write.csv(seq_data, file = file.path(outdir, "seq_data.csv"), quote = FALSE, row.names = FALSE)
write.csv(binned_seq_data, file = file.path(outdir, "binned_seq_data.csv"), quote = FALSE, row.names = FALSE)

