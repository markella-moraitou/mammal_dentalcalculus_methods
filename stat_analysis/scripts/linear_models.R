##### LINEAR MODELS #####

#### (Generalised) Linear Mixed-effects models

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggnewscale)
library(vegan)
library(lme4)
library(lmerTest)
library(rstatix)
library(car)
library(blmeco)
library(FDB1)
library(ggfortify)
library(ggpubr)
library(cowplot)
library(multcomp)
library(EcolUtils)
library(scales)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory paths
indir <- normalizePath(file.path("..","input"))
outdir <- normalizePath(file.path("..","output", "LM")) # output directory

# Create directory for output
dir.create(outdir, recursive  =  TRUE, showWarnings  =  FALSE)

# Theme
source("ggplot_theme.R")
theme_set(theme_bw_alt)

#######################
#####  LOAD INPUT #####
#######################

lab_data <- read.csv(file = file.path(indir, "calculus_lab_metadata.csv"))

species_info <- read.csv(file = file.path(indir, "species_info.csv"))
diet <- read.csv(file = file.path(indir, "Lintulaakso_diet_filtered.csv"))

bioinfo_data <- read.csv(file = file.path(indir, "bioinfo_metadata.csv"))

# Load palettes
for (file in dir(file.path(indir, "palettes"))) {
  name <- str_remove(file, ".csv")
  palette <- read.csv(file.path(indir, "palettes", file)) %>% pull(2)
  names(palette) <-  read.csv(file.path(indir, "palettes", file)) %>% pull(1)
  assign(name, palette)
}

#######################
####    PROCESS    ####
#######################

#### Process input table ####

# Merge with species info and diet information from Lintulaakso paper
data <- lab_data %>% left_join(species_info, by = "Species") %>% 
  left_join(diet, by="Species") %>%
  left_join(bioinfo_data, by = "Ext.ID")

colnames(data) <- gsub("..", ".", colnames(data), fixed  =  TRUE)

# Order species will be displayed in
sp_levels <- data %>% dplyr::select(Species, order) %>% unique %>% arrange(order, Species) %>% pull(Species)

data <- as_tibble(data) %>% 
  # Add some columns
  mutate(is.blank  =  grepl("blank", Species),
                                   failed.repeated  =  grepl("F_", Ext.ID)) %>%
  # Combine the rows for samples that come from multiple specimens
  group_by(Ext.ID) %>%
  mutate(combined.sample = (n_distinct(Catalogue.accession)>1)) %>% 
  summarize_all( ~ ifelse(length(unique(.)) > 1, paste(., collapse  =  "/"), .)) %>%
  # Make sure Qubit concentration is read as numbers and that species and diet are factors
  # Qubit measuremnts that were "too low" too be read are replaced with 0.05 which is the sensitivity threshold
  mutate(Qubit.concentration.ng.ul. = as.numeric(Qubit.concentration.ng.ul.) %>% ifelse(is.na(.), 0.049, .),
         # Add a pseudocount on weight
         Sample.weight.g. = case_when(is.na(Sample.weight.g.) ~ 0, TRUE ~ Sample.weight.g.) + 0.001,
         Species = factor(Species, levels = sp_levels),
         calculated_species_main_diet = factor(calculated_species_main_diet, levels = c("Herbivore", "Frugivore", "Animalivore")),
         Sample.weight.mg = Sample.weight.g. * 1000) %>%
  # Make sure batches are read as integers
  mutate(Ext.Date = as.character(Ext.Date),
         LP.date = as.character(LP.date),
         Ind.date = as.character(Ind.date)) %>%
  # Get collection year as numeric (so we can test for the effect of sample age)
  mutate(Year_most_recent = case_when(is.na(as.numeric(Year)) ~ as.numeric(str_remove_all(Year, "<|\\?| \\(received\\)|[0-9][0-9][0-9][0-9]-")), TRUE ~ as.numeric(Year))) %>%
  mutate(Age_approximated = (Year != Year_most_recent)) %>%
  mutate(Age = (2023-Year_most_recent)/100,
        Year_most_recent = NULL)

diet <- diet %>%
  mutate(calculated_species_main_diet = factor(calculated_species_main_diet, levels = c("Herbivore", "Frugivore", "Animalivore")))

# Calculate some variables we need

# First, calculate PCR efficiency the average amount by which the copies are multiplied at every step
# Ind.copies in per ul, but the final vol of the indexing reaction is 20ul
# So the absolute number of indexing copies is Ind.copies*20

elution_vol_ul = 45
barcoding_final_vol_ul = 40
indexing_final_vol_ul = 20

data <- data %>% ungroup %>% 
  mutate(DNA_output_ug = Qubit.concentration.ng.ul.*0.001 * elution_vol_ul, # DNA output in μg: concentration in μg x elution volume in ul
         DNA_input_ug = Qubit.concentration.ng.ul.*0.001 * Extract.vol.used, # DNA input in μg: concentration in μg x vol used for barcoding in ul
         # Same for barcoded copies
         bc_output_copies = LP.copies * barcoding_final_vol_ul,
         bc_input_copies = LP.copies * LP.vol.used,
         ind_output_copies = Ind.copies*indexing_final_vol_ul,
         efficiency = exp(log(ind_output_copies/bc_input_copies)/Ind.cycles) -1)

eff = mean(data$efficiency[!data$is.blank], na.rm=TRUE)

data <- data %>%
  # theoretical indexing copies depending on PCR cycles
  mutate(exp_ind_copies = bc_input_copies * (1 + eff)^Ind.cycles) %>%
  # failed and repeated samples
  filter(failed.repeated == FALSE)

#Save data
write.table(data, file = file.path(outdir, "data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Exclude  as well as blanks
filt_dat <- data %>% filter(is.blank == FALSE & !is.na(Sample.weight.mg))

#######################
#### Summary plots ####
#######################

## Flowchart

flow <- data.frame(x = rep(c(1, 1, 2), 4),
                   y = c(1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8),
                   type = rep(c("Protocol step", "Quantification", "Value"), 4),
                   group = c(rep("Sample", 3), rep("DNA", 3), rep("Barcoding", 3), rep("Indexing", 3)),
                   name = c("Sampling", "Weighing", "Weight (mg)",
                            "DNA extraction", "Qubit dsDNA HS Assay", "DNA amount (μg)",
                            "Barcoding", "Post-barcoding qPCR","Barcoded libraries (molecules)",
                            "Indexing", "Post-indexing qPCR", "Indexed libraries (molecules)"))


flow <- flow %>% 
  mutate(xmin = x - 0.35,
         xmax = x + 0.35,
         ymin = y - 0.25,
         ymax = y + 0.25)

flowchart <-
  ggplot(data = flow, aes(x = x, y = y, fill = type, label = name)) +
  scale_y_reverse() +
  geom_rect(aes(xmin = xmin, ymin = ymin, 
                xmax = xmax, ymax = ymax), alpha = 0.5, colour = "black") +
  geom_text() +
  theme_void_alt +
  scale_fill_manual(name = "", values = c("Protocol step" = "#1F8A02", "Quantification" = "#BBE5AF", "Value" = "white")) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.text.y = element_blank()) +
  labs(title = "DNA extraction and library preparation flowchart",
       subtitle = "indicating measurement steps and variables used within models")


#Can't get arrows to work -- will add them later
ggsave(filename  =  file.path(outdir, "flowchart.png"), flowchart, width  =  10, height = 8)

#### Diet PCA ####
diet_var <- diet[c("cp", "ee", "cf", "ash", "nfe")]
rownames(diet_var) <- filt_dat$Common.name[match(diet$Species, filt_dat$Species)]

ord <- prcomp(diet_var)

# Scale PCs between 0 and 1
scale_01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

ord_df <- data.frame(ord$x) %>% mutate(PC1_scaled = scale(PC1, center=0)[,1],
                                       PC2_scaled = scale(PC2, center=0)[,1])

# Get some info for plotting
var_explained <- round(ord$sdev^2 * 100 / sum(ord$sdev^2), 1) # Variance explained

loadings_matrix <- data.frame(Variables = rownames(ord$rotation[,c(1,2)]), ord$rotation[,c("PC1", "PC2")])
colnames(loadings_matrix)[2:3] <- c("PC1_scaled", "PC2_scaled")

loadings_matrix <- loadings_matrix %>%
        mutate(x = case_when(PC1_scaled > 0 ~ PC1_scaled + 0.2, TRUE ~ PC1_scaled - 0.2), # Scale the arrows
               y = case_when(PC2_scaled > 0 ~ PC2_scaled + 0.1, TRUE ~ PC2_scaled - 0.1)) %>%
               filter(Variables != "ash") # Remove ash as it is not relevant for the PCA

labels <- ord_df %>% as.data.frame %>% mutate(label = case_when(PC1 > 10 ~ rownames(.),
                                                               PC2 == max(PC2) | PC2 == min(PC2) ~ rownames(.),
                                                               rownames(.) == "Olive baboon" ~ rownames(.))) %>%
          mutate(y = case_when(is.na(label) ~ NA,
                              label %in% c("Olive baboon", "South American fur seal", "Orca", "Bornean orangutan") ~ PC2_scaled - 0.1,
                              label %in% c("South American sea lion", "African elephant", "European badger") ~ PC2_scaled + 0.1)) %>%
  dplyr::select(label, y)

# Plot
pca <- ggplot(aes(x = PC1_scaled, y = PC2_scaled, colour=diet$calculated_species_main_diet), data = ord_df) +
  geom_point(size=3, alpha = 0.8) +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  xlab(paste("scaled PC1 -", var_explained[1], "%")) +
  ylab(paste("scaled PC2 -", var_explained[2], "%")) +
  geom_segment(data = loadings_matrix, aes(x = 0, y = 0, xend = PC1_scaled,
                                       yend = PC2_scaled), arrow = arrow(length = unit(0.5, "picas")),
               color = "black") +
  theme_bw_alt + theme(legend.position = "bottom") +
  annotate("text", x = loadings_matrix$x, y = loadings_matrix$y,
           label = loadings_matrix$Variables, size = 8, alpha = 0.8) +
  labs(tag = "A.")

## Add principal components to data
filt_dat$PC1 <- ord_df[match(filt_dat$Common.name, rownames(ord_df)), "PC1"]
filt_dat$PC2 <- ord_df[match(filt_dat$Common.name, rownames(ord_df)), "PC2"]

filt_dat$PC1_scaled <- ord_df[match(filt_dat$Common.name, rownames(ord_df)), "PC1_scaled"]
filt_dat$PC2_scaled <- ord_df[match(filt_dat$Common.name, rownames(ord_df)), "PC2_scaled"]

# Based on this, combine omnivores and frugivores together
diet <- diet %>% mutate(diet_category = case_when(Species == "Papio anubis" ~ "Frugivore",
                                                          TRUE ~ calculated_species_main_diet))
diet$diet_category <- factor(diet$diet_category, levels = c("Animalivore", "Frugivore", "Herbivore"))

filt_dat <- filt_dat %>% left_join(diet)

# Save table
write.table(filt_dat, file = file.path(outdir, "filtered_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Fix PCA
pca <- 
  pca + geom_point(size=3, aes(colour = diet$diet_category)) +
  scale_colour_manual(values = diet2_palette[unique(as.character(diet$diet_category))], name = "Dietary category") +
  geom_text(label = labels$label, size = 6, hjust = 0.5, aes(x = PC1_scaled, y = labels$y, colour = diet$diet_category)) +
    xlim(-1, 4)

ggsave(filename  =  file.path(outdir, "dietary_PCA.png"), pca, width  =  8, height = 8)

# Save just dietary data
diet_data <- filt_dat %>% dplyr::select(Species, Common.name, order,
                                 Species_Lintulaakso_et_al,diet_cluster_name, diet_quality,
                                 cp, ee, cf, ash, nfe, calculated_species_main_diet, PC1, PC2) %>%
  unique

diet_data$PC1_scaled <- ord_df[match(diet_data$Common.name, rownames(ord_df)), "PC1_scaled"]
diet_data$PC2_scaled <- ord_df[match(diet_data$Common.name, rownames(ord_df)), "PC2_scaled"]

# Save table
write.table(diet_data, file = file.path(outdir, "diet_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

#### Diet summary ####
diet_long <- diet %>% left_join(species_info) %>% pivot_longer(c(cp, ee, cf, ash, nfe), values_to = "proportion", names_to = "nutrient") %>%
  mutate(nutrient = factor(nutrient, levels = c("ash", "ee", "cp", "nfe", "cf"))) %>%
  arrange(nutrient)

nutrient_palette <- c(ash = "grey",
                      cf = "#3BA01B",
                      nfe = "#3A459C",
                      cp = "#B41F34",
                      ee = "#BD9120")

diet_barplot <- ggplot(diet_long, aes(x = Common.name, y = proportion, fill = nutrient, group = diet_category)) +
  geom_bar(stat="identity") +
  theme_bw_alt + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top") +
  facet_grid(cols = vars(diet_category), scales = "free", space = "free") +
  scale_fill_manual(values = nutrient_palette,
                    labels = c("ash" = "ash (inorganic)",
                               "cf" = "CF (crude fibre)",
                               "nfe" = "NFE (carbohydrates)",
                               "cp" = "CP (crude protein)",
                               "ee" = "EE (fat)"),
                    name = "") +
  xlab("Species") 

diet_barplot

ggsave(filename  =  file.path(outdir, "diet_barplot.png"), diet_barplot, width  =  12, height = 8)

#### Dataset summary ####

## Species summary
species_summ <-
  filt_dat %>% group_by(Common.name, Species, order, diet_category) %>% 
  summarise(n_extracted = n_distinct(Ext.ID[!is.na(Ext.Date)]),
            n_barcoded = n_distinct(Ext.ID[!is.na(LP.date)]),
            n_indexed = n_distinct(Ext.ID[Ind.date != ""]))

write.csv(species_summ, file = file.path(outdir, "species_summary.csv"), quote  =  FALSE, row.names  =  FALSE)

## Fisher exact test for order and diet
fisher.test(table(species_summ$diet_category, species_summ$order))

## Order and diet
summ_dat <- filt_dat %>% group_by(order, diet_category) %>% 
  summarise(n_samples = n_distinct(Ext.ID), n_species = n_distinct(Species))

write.csv(summ_dat, file = file.path(outdir, "data_summary.csv"), quote  =  FALSE, row.names  =  FALSE)

sums <- summ_dat %>% ungroup %>% summarise(total_species = sum(n_species), total_samples = sum(n_samples))

# Add label colors
summ_dat <- summ_dat %>% mutate(lab_col = case_when(n_samples < 25 ~ "black", TRUE ~ "white"))

p_summary <-
  ggplot(aes(x = diet_category, y = order, size = n_species, fill = n_samples), data = summ_dat) +
  geom_point(shape = 21, colour = "black") +
  theme_void_alt + theme(legend.title.position = "top", axis.text.y = element_text(hjust = 0)) +
  xlab("Dietary category") + ylab("Taxonomic order") +
  labs(tag = "B.") +
  scale_fill_continuous(low = "lightblue", high = "darkblue", name = "Sample number", trans = "log", breaks = c(5, 25, 125), 
                        guide = guide_colourbar(direction = "horizontal")) +
  scale_size(range  =  c(15, 35), breaks = c(min(summ_dat$n_species), mean(unique(summ_dat$n_species)), max(summ_dat$n_species)),
             name = "Species number") +
  new_scale(new_aes = "size") +
  geom_text(aes(label = paste0(n_species, " (", n_samples, ")"), x = diet_category, y = order, colour = lab_col, size = n_species), fontface = "bold") +
  scale_colour_manual(values = summ_dat$lab_col, guide = "none") +
  scale_size_continuous(range = c(4, 8), guide = "none")

ggsave(filename  =  file.path(outdir, "summary_plot.png"), p_summary, width  =  9, height = 8)

## Coloured samples

extract_palette <- c("TRUE"="#660000", "FALSE"="#FFCC66")

# Calculate proportions of pigmented extracts
props <- filt_dat %>% group_by(diet_category, Pigmented_extract) %>% summarise(n_samples = n_distinct(Ext.ID)) %>%
  ungroup(Pigmented_extract) %>% 
  mutate(prop_coloured = round(n_samples * 100 / sum(n_samples), digits = 1)) %>%
  mutate(label=paste(prop_coloured, "%")) %>%
  mutate(yvalues=case_when(Pigmented_extract == TRUE ~ n_samples + 5)) # Y values for text labels

p_col <-
  ggplot(aes(x = diet_category, fill=Pigmented_extract), data=filt_dat) +
  geom_bar() +
  geom_text(data=props, aes(label=label, y=yvalues, colour=Pigmented_extract)) +
  scale_fill_manual(values = extract_palette, name="Extract pigmented") +
  scale_colour_manual(values = extract_palette, name="Extract pigmented") +
  theme_bw_alt + theme(legend.position = "top") +
  ylab("Number of samples") + xlab("Taxonomic order") +
  labs(tag = "B.")

p_col + theme(axis.text.x = element_text(angle = 90))

ggsave(filename  =  file.path(outdir, "coloured_samples.png"), p_col, width  =  7, height = 8)

PC1_col <-
  ggplot(aes(y = PC1, x = Pigmented_extract, fill = Pigmented_extract), data = filt_dat) +
  geom_boxplot()

PC1_col

wilcox.test(filt_dat$PC1[filt_dat$Pigmented_extract==TRUE], filt_dat$PC1[filt_dat$Pigmented_extract==FALSE])

PC2_col <-
  ggplot(aes(y = PC2, x = Pigmented_extract, fill = Pigmented_extract), data = filt_dat) +
  geom_boxplot()

PC2_col

wilcox.test(filt_dat$PC2[filt_dat$Pigmented_extract==TRUE], filt_dat$PC2[filt_dat$Pigmented_extract==FALSE])

#### Assess indexing using blanks ####

cutoff <- 10^9 # Approximate blanks value, cutoff for bad and good samples
bad_samples <- data %>% filter(ind_output_copies > 0 & ind_output_copies < cutoff & !is.blank) %>% nrow
good_samples <- data %>% filter(ind_output_copies > cutoff & !is.blank) %>% nrow
blanks <- data %>% filter(ind_output_copies > 0, is.blank) %>% nrow

hist <-
  data %>% ggplot(aes(x = ind_output_copies, fill = is.blank, colour = Pigmented_extract)) +
  scale_x_continuous(trans = "log10") +
  geom_histogram(bins = 50, linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "darkgrey", "FALSE" = "#FFCC66"), name = "", labels = c("True sample", "Negative control")) +
  scale_colour_manual(values = c("TRUE" = "#660000", "FALSE" = "transparent"), name = "", labels = c("TRUE" = "Pigmented", "FALSE" = "Not pigmented")) +
  xlab("Indexed library molecules") + ylab("count") +
  geom_vline(xintercept = cutoff, linetype = "dashed") +
  annotate("text", x = 10^7, y = 25, label = paste0("Failed indexing:\n n = ", bad_samples), colour = "#C17D21") +
  annotate("text", x = 10^10, y = 30, label = paste0("Successful indexing:\n n = ", good_samples), colour = "#C17D21") +
  annotate("text", x = 10^9, y = 16, label = paste0("Indexed blanks:\n n = ", blanks)) +
  theme(legend.position = "top") +
  labs(tag = "A.")

hist

ggsave(filename  =  file.path(outdir, "indexing_hist.png"), hist, width  =  8, height = 8)

#######################
#### MIXED MODELS  ####
#######################

# Create a set of diagnostic plots for GLMM
diagnose_glmm <- function(model) {
  diag_plots <- list() # Store all diagnostic plots here
  
  # Assumption: residuals shouldn't be correlated to fitted values
  tbl <- data.frame(fitted.values = fitted.values(model), residuals = residuals(model))
  ranf_fitted <- ggplot(data = tbl, aes(x = fitted.values, y = residuals)) + geom_point() + theme_bw() +
    labs(title = "Model residuals to fitted values")
  
  diag_plots[["ranef_fitted"]] <- ranf_fitted
  
  # Assumption: Random effects come from a normal distribution
  ranef <- ggCaterpillar(ranef(model), QQ=FALSE)
  diag_plots <- c(diag_plots, ranef)
  #for (i in 1:length(ranef(model))) {
  #  ranef <- ranef(model)[[i]] # Extract random effects
  #  random_var <- names(ranef(model))[i] # Variable used as random effect
  #  response <- colnames(ranef)[1] # Response variable
  #  # Plot qqplot
  #  ranef_qqplot <- ggplot(data = ranef,
  #                         aes(sample = get(response))) +
  #    geom_qq() + geom_qq_line() + theme_bw() + xlab("Theoretical quantiles") + ylab("Observed quantiles") +
  #    labs(title = paste("Q-Q Plot for random effects of", random_var, "on", response))
  #  diag_plots[[random_var]] <- ranef_qqplot
  #}
  return(diag_plots)
}

# Create a set of diagnostic plots for LMM
diagnose_lmm <- function(model) {
  diag_plots <- list() # Store all diagnostic plots here
 
  # Assumption: residuals shouldn't be correlated to fitted values
  tbl <- data.frame(fitted.values = fitted.values(model), residuals = residuals(model))
  ranf_fitted <- ggplot(data = tbl, aes(x = fitted.values, y = residuals)) + geom_point() + theme_bw() +
    labs(title = "Model residuals to fitted values")
  
  diag_plots[["ranef_fitted"]] <- ranf_fitted
  
  # Assumption: residuals are normal
  res_qqplot <- ggplot(data = data.frame(residuals = residuals(model)),
                         aes(sample = residuals)) +
    geom_qq() + geom_qq_line() + theme_bw() + xlab("Theoretical quantiles") + ylab("Observed quantiles") +
    labs(title = paste("Q-Q Plot for model residuals"))
  
  diag_plots[["res_qqplot"]] <- res_qqplot
  # Assumption: Random effects come from a normal distribution
  ranef <- ggCaterpillar(ranef(model), QQ=FALSE)
  diag_plots <- c(diag_plots, ranef)
  #for (i in 1:length(ranef(model))) {
  #  ranef <- ranef(model)[[i]] # Extract random effects
  #  random_var <- names(ranef(model))[i] # Variable used as random effect
  #  response <- colnames(ranef)[1] # Response variable
  #  # Plot qqplot
  #  ranef_qqplot <- ggplot(data = ranef,
  #                         aes(sample = get(response))) +
  #    geom_qq() + geom_qq_line() + theme_bw() + xlab("Theoretical quantiles") + ylab("Observed quantiles") +
  #    labs(title = paste("Q-Q Plot for random effects of", random_var, "on", response))
  #  diag_plots[[random_var]] <- ranef_qqplot
  #}
  return(diag_plots)
}

#### 1. Weight - Concentration ####
# Does diet affect how much DNA we extract from a given weight?
filt_dat1 <- filt_dat %>% filter(DNA_output_ug < 4)

# Plot
p_m1 <-
  ggplot(aes(x = Sample.weight.mg, y = DNA_output_ug, colour = diet_category, shape=Pigmented_extract), data = filt_dat1) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = diet2_palette, name = "Dietary category") +
  #scale_colour_manual(values = diet_palette, name = "Dietary category") +
  scale_shape_manual(values = c(16, 1), name = "Pigmented extract") +
  xlab("Sample weight (mg)") + ylab("DNA output (μg)") +
  theme_bw_alt +
  labs(tag = "A.")

# Simple linear model
lm1 = lm(DNA_output_ug ~
             Sample.weight.g. + Sample.weight.g.:PC1_scaled + Sample.weight.g.:PC2_scaled + Sample.weight.g.:Pigmented_extract,
           data = filt_dat1)

# Extract weights
weights <- 1 / lm(abs(lm1$residuals) ~ lm1$fitted.values)$fitted.values^2

mm1 = lmer(DNA_output_ug ~
             Sample.weight.g. + Sample.weight.g.:PC1_scaled + Sample.weight.g.:PC2_scaled + Sample.weight.g.:Pigmented_extract +
              (0 + Sample.weight.g. | Species) + (0 + Sample.weight.g.| Ext.Date),
            data = filt_dat1, weights = weights)

# Plot diagnostics
diagn_1 <- diagnose_lmm(mm1)

ggsave(filename  =  file.path(outdir, "diagnostics_mmodel1.png"), plot_grid(plotlist = diagn_1),
       width  =  14, height = 14)

vif(mm1)
shapiro.test(residuals(mm1))
hist(residuals(mm1))
s1 <- summary(mm1)
s1

# I will rescale this for ug of Sample instead of g of sample for plotting
# Output
weight_1 <- s1$coefficients["Sample.weight.g.", "Estimate"]/1000 # Qubit increases so much for every 1ng of sample added
pc1_1 <- s1$coefficients["Sample.weight.g.:PC1_scaled", "Estimate"]/1000 # Reflects the interaction between PC1 and sample weight

# Get PC1 values for plotting 
pc1_val <- data.frame(value=c(-0.5, 2.5)) #PC1 values to plot
pc1_val$lab <- paste0("scaled PC1 = ", pc1_val$value)

# Get formula from significant factors for two pc1 values
# Not carnivorous, low PC1
fun1_1 = function(x) {y = (weight_1 + pc1_1*pc1_val$value[1])*x}
# Carnivorous, high PC1
fun1_2 = function(x) {y = (weight_1 + pc1_1*pc1_val$value[2])*x}

p_m1 <-
  p_m1 + 
  geom_text(aes(x = 20, y = 3), alpha=0.5, colour="black", size = 8,
            label = paste0("y = (", round(weight_1, 2), " + ", round(pc1_1, 3), " * PC1", ") * x")) +
  geom_function(fun = fun1_1, size = 1, colour = diet2_palette["Herbivore"]) +
  geom_function(fun = fun1_2, size = 1, colour = diet2_palette["Animalivore"]) +
  geom_label(aes(x = 30, y = 1.5, label = paste(pc1_val$lab[2],"\n(Animalivorous)")), colour = "black", size=5, alpha=0.5) +
  geom_label(aes(x = 30, y = 0.5, label = paste(pc1_val$lab[1],"\n(Herbi-/Frugivorous)")), colour = "black", size=5, alpha=0.5) 

ggsave(filename  =  file.path(outdir, "plot_mmodel1.png"), p_m1, width  =  15, height = 10)

## Box plot with average yields

filt_dat1 <- filt_dat1 %>% mutate(DNA_yield = DNA_output_ug / Sample.weight.mg)

# Get statistical tests (overall and pairwise)

# Run kruskal wallis test
pval <- kruskal_test(DNA_yield ~ diet_category, data = filt_dat1)$p
effsize <- kruskal_effsize(DNA_yield ~ diet_category, data = filt_dat1)$effsize
annot <- paste0("Kruskal-Wallis:\neffect size = ", round(effsize, 2), ",\np-value = ", pval)

model1_box <-
  ggviolin(data = filt_dat1, x = "diet_category", y = "DNA_yield", fill = "diet_category",
           add = "mean_se") +
  scale_fill_manual(values = diet2_palette, name = "") +
  stat_compare_means(comparisons = list(c("Herbivore", "Frugivore"), c("Frugivore", "Animalivore"), c("Herbivore", "Animalivore")),
                     label = "p.signif", method = "wilcox", size = 8) +
  annotate("text", label = annot, y = 0.1, x = 1.2, size = 6, hjust = 0) +
  ylab("DNA yield (μg per mg of dental calculus)") + xlab("Dietary category") +
  theme_bw_alt + theme(legend.position = "none") +
  labs(tag = "B.")

ggsave(filename  =  file.path(outdir, "model1_boxplot.png"), model1_box, width  =  15, height = 8)

## Also, plot relationship between PC1 and PC2 and DNA concentration
filt_dat1 <- filt_dat1 %>%
  # Bin the PC1 and PC2 values
  mutate(PC1_binned = cut(PC1, breaks = seq(-20, 80, by = 10), ordered_result = TRUE),
         PC2_binned = cut(PC2, breaks = seq(-14, 14, by = 2.5), ordered_result = TRUE))

pc1_box <-
  ggplot(data = filt_dat1, aes(y = DNA_yield, x = PC1_binned)) +
  geom_boxplot()  +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) + theme(axis.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
  xlab("Dietary PC1") + ylab("") + labs(tag = "C.")

pc2_box <-
  ggplot(data = filt_dat1, aes(y = DNA_yield, x = PC2_binned)) +
  geom_boxplot()  +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) + theme(axis.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
  xlab("Dietary PC2") + ylab("") + labs(tag = "D.")

pc_to_conc <- plot_grid(pc1_box, pc2_box)

ggsave(filename  =  file.path(outdir, "PC_to_DNA_yield.png"), pc_to_conc, width  =  8, height = 8)

diet_to_conc <- plot_grid(model1_box, pc1_box, pc2_box,
                          rel_widths = c(2, 1, 1), nrow = 1, 
                          align = "h", axis = "tb")

ggsave(filename  =  file.path(outdir, "diet_to_DNA_yield.png"), diet_to_conc, width  =  16, height = 8)

## Consider effect of age (for subset of samples with collection year)
filt_dat1_age <- filt_dat1 %>% filter(!is.na(Age))

lm1_age = lm(DNA_output_ug ~
             Sample.weight.g. + Sample.weight.g.:PC1_scaled + Sample.weight.g.:PC2_scaled + Sample.weight.g.:Pigmented_extract + Sample.weight.g.:Age,
           data = filt_dat1_age)

# Extract weights
weights <- 1 / lm(abs(lm1_age$residuals) ~ lm1_age$fitted.values)$fitted.values^2

mm1_age = lmer(DNA_output_ug ~
             Sample.weight.g. + Sample.weight.g.:PC1_scaled + Sample.weight.g.:PC2_scaled + Sample.weight.g.:Pigmented_extract + Sample.weight.g.:Age +
              (0 + Sample.weight.g. | Species) + (0 + Sample.weight.g. | Ext.Date),
            data = filt_dat1_age, weights = weights)

# Plot diagnostics
diagn_1_age <- diagnose_lmm(mm1_age)

ggsave(filename  =  file.path(outdir, "diagnostics_mmodel1_age.png"), plot_grid(plotlist = diagn_1),
       width  =  14, height = 14)

vif(mm1_age)
shapiro.test(residuals(mm1_age))
hist(residuals(mm1_age))

summary(mm1_age)

#### 2. Concentration - LP copy numbers ####

filt_dat2 <- filt_dat %>%
  # Get output in 10^10
  mutate(bc_output_e10 = bc_output_copies * 10^(-10)) %>%
  # Remove outliers
  filter(!is.na(LP.copies)) %>%
  # Remove outliers
  filter(bc_output_e10 < 8 & DNA_input_ug < 3)

# Plot
p_m2 <-
  ggplot(aes(x = DNA_input_ug, y = bc_output_e10, colour = diet_category, shape=Pigmented_extract), data = filt_dat2) +
  geom_point(alpha = 0.7, size = 3) +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  #scale_colour_manual(values = diet_palette, name = "Dietary category") +
  scale_shape_manual(values = c(16, 1), name = "Pigmented extract") +
  scale_x_continuous(name = "DNA input amount (μg)") + 
  scale_y_continuous(name = expression(paste("Barcoded library molecules (10"^10, ")"))) +
  theme_bw_alt

# Simple linear model
lm2 = lm(bc_output_e10 ~
           DNA_input_ug + DNA_input_ug:PC1_scaled + DNA_input_ug:PC2_scaled + DNA_input_ug:Pigmented_extract,
         data = filt_dat2)

# Extract weights
weights <- 1 / lm(abs(lm2$residuals) ~ lm2$fitted.values)$fitted.values^2

mm2 <-  lmer(bc_output_e10 ~
               DNA_input_ug + DNA_input_ug:PC1_scaled + DNA_input_ug:PC2_scaled + DNA_input_ug:Pigmented_extract +
               (0 + DNA_input_ug | Species) + (0 + DNA_input_ug | LP.date),
              data = filt_dat2, weights = weights)

# Plot diagnostics
diagn_2 <- diagnose_lmm(mm2)
ggsave(filename  =  file.path(outdir, "diagnostics_mmodel2.png"), plot_grid(plotlist = diagn_2),
       width  =  14, height = 14)

vif(mm2)
shapiro.test(residuals(mm2))
hist(residuals(mm2))

s2 <- summary(mm2)
s2

# Output
int_2 <- s2$coefficients["(Intercept)", "Estimate"]
dna_2 <- s2$coefficients["DNA_input_ug", "Estimate"]
colour_2 <- s2$coefficients["DNA_input_ug:Pigmented_extractTRUE", "Estimate"]

fun2_1 <- function(x) {int_2 + (dna_2 + colour_2)*x}
fun2_2 <- function(x) {int_2 + (dna_2)*x}

p_m2 <-
  p_m2 + 
  geom_text(aes(x = 0.5, y = 7), alpha=0.5, colour="black", size = 8,
            label = paste0("y = ", round(int_2, 1), " + (", round(dna_2, 1), " ", round(colour_2, 1), "*pigmentation)*x")) +
    geom_function(fun = fun2_1, linetype = "dashed", size = 1, xlim = c(0, 1), colour = "black") +
    geom_function(fun = fun2_2, size = 1, xlim = c(0, 1), colour = "black") +
    geom_label(aes(x = 1, y = 3, label = "Not pigmented"), size=8, alpha=0.5, colour = "black") +
    geom_label(aes(x = 1, y = 1, label = "Pigmented"), size=8, alpha=0.5, colour = "black") +
   labs(tag = "A.")

p_m2

ggsave(filename  =  file.path(outdir, "plot_mmodel2.png"), p_m2, width  =  15, height = 8)

## Box plot with average yields

filt_dat2 <- filt_dat2 %>% mutate(Library_yield = bc_output_e10 / DNA_input_ug) %>%
  # For labelling
  mutate(Extract = case_when(Pigmented_extract == TRUE ~ "Pigmented",
                             TRUE ~ "Clear"))

# Run kruskal wallis test
# For clear extracts
pval_c <- kruskal_test(Library_yield ~ diet_category, data = filter(filt_dat2, Extract == "Clear"))$p
effsize_c <- kruskal_effsize(Library_yield ~ diet_category, data = filter(filt_dat2, Extract == "Clear"))$effsize
annot_c <- paste0("Kruskal-Wallis:\neffect size = ", round(effsize_c, 2), ",\np-value = ", pval_c)

# For pigmented extracts
pval_p <- kruskal_test(Library_yield ~ diet_category, data = filter(filt_dat2, Extract == "Pigmented" & diet_category != "Animalivore"))$p
effsize_p <- kruskal_effsize(Library_yield ~ diet_category, data = filter(filt_dat2, Extract == "Pigmented" & diet_category != "Animalivore"))$effsize
annot_p <- paste0("Kruskal-Wallis:\neffect size = ", round(effsize_p, 2), ",\np-value = ", round(pval_p, 3))

ann_text <- data.frame(lab = c(annot_p, annot_c),
                       diet_category = c(0.2, 0.2),
                       Library_yield = c(10^(-4), 10^(-4)),
                       Extract = c("Pigmented", "Clear"))

# Plot
model2_box <-
  ggviolin(data = filter(filt_dat2, Extract != "Pigmented" | diet_category != "Animalivore"),
           x = "diet_category", y = "Library_yield", fill = "diet_category",
           add = "mean_se") +
  scale_fill_manual(values = diet2_palette, name = "") +
  stat_compare_means(comparisons = list(c("Animalivore", "Frugivore"),
                                        c("Animalivore", "Herbivore"),
                                        c("Frugivore", "Herbivore")),
                     label = "p.signif", size = 8) +
  geom_text(data = ann_text, aes(label = lab), size = 6, hjust = 0) +
  facet_grid(cols = vars(Extract)) +
  scale_y_continuous(trans = "log10") +
  ylab(expression(paste("Barcoding yield (", 10^10, "molecules per μg of DNA)"))) + xlab("") +
  theme_bw_alt + theme(legend.position = "none") +
  labs(tag = "A.")

ggsave(filename  =  file.path(outdir, "model2_boxplot.png"), model2_box, width  =  10, height = 10)

# Kruskal-Wallis tests

kruskal.test(data = filt_dat2, Library_yield ~ Pigmented_extract)
kruskal_effsize(data = filt_dat2, Library_yield ~ Pigmented_extract)

kruskal.test(data = filt_dat2, Library_yield ~ diet_category)
kruskal_effsize(data = filt_dat2, Library_yield ~ diet_category)

## Consider effect of age (for subset of samples with collection year)
filt_dat2_age <- filt_dat2 %>% filter(!is.na(Age))

# Simple linear model
lm2_age = lm(bc_output_e10 ~
           DNA_input_ug + DNA_input_ug:PC1_scaled + DNA_input_ug:PC2_scaled + DNA_input_ug:Pigmented_extract + DNA_input_ug:Age,
         data = filt_dat2_age)

# Extract weights
weights <- 1 / lm(abs(lm2_age$residuals) ~ lm2_age$fitted.values)$fitted.values^2

mm2_age <-  lmer(bc_output_e10 ~
               DNA_input_ug + DNA_input_ug:PC1_scaled + DNA_input_ug:PC2_scaled + DNA_input_ug:Pigmented_extract + DNA_input_ug:Age +
               (0 + DNA_input_ug | Species) + (0 + DNA_input_ug | LP.date),
              data = filt_dat2_age, weights = weights)

# Plot diagnostics
diagn_2_age <- diagnose_lmm(mm2_age)

ggsave(filename  =  file.path(outdir, "diagnostics_mmodel2_age.png"), plot_grid(plotlist = diagn_1),
       width  =  14, height = 14)

vif(mm2_age)
shapiro.test(residuals(mm2_age))
hist(residuals(mm2_age))

summary(mm2_age)

#### 3. Theoretic - Real ind numbers ####

filt_dat3 <- filt_dat2 %>% filter(!(is.na(Ind.copies))) %>%
  # Get expected and observed copies in 10^10
  mutate(observed_e10 = ind_output_copies/10^10,
         expected_e10 = exp_ind_copies/10^10) %>%
  # Remove
  filter(observed_e10 < 400)

# Plot
p_m3 <-
  ggplot(aes(x = expected_e10, y = observed_e10, colour = diet_category, shape = Pigmented_extract), data = filt_dat3) +
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  #scale_colour_manual(values = diet_palette, name = "Dietary category") +
  scale_shape_manual(values = c(16, 1), name = "Pigmented extract") +
  scale_x_continuous(name = expression(paste("Expected indexed molecules (", 10^10, ")"))) + 
  scale_y_continuous(name = expression(paste("Observed indexed molecules (", 10^10, ")")))+
  theme_bw_alt

p_m3

# Simple linear model
lm3 = lm(observed_e10 ~
           expected_e10 + expected_e10:PC1_scaled + expected_e10:PC2_scaled + expected_e10:Pigmented_extract,
         data = filt_dat3)

# Extract weights
weights <- 1 / lm(abs(lm3$residuals) ~ lm3$fitted.values)$fitted.values^2

mm3 <-  lmer(observed_e10 ~
               expected_e10 + expected_e10:PC1_scaled + expected_e10:PC2_scaled + expected_e10:Pigmented_extract +
               (0 + expected_e10 | Species) + (0 + expected_e10 | Ind.date),
             data = filt_dat3, weights = weights)

# Plot diagnostics
diagn_3 <- diagnose_lmm(mm3)
ggsave(filename  =  file.path(outdir, "diagnostics_mmodel3.png"), plot_grid(plotlist = diagn_3),
       width  =  14, height = 14)

vif(mm3)
shapiro.test(residuals(mm3))
hist(residuals(mm3))

s3 <- summary(mm3)
s3

# Output
int_3 <- s3$coefficients["(Intercept)", "Estimate"]
exp_3 <- s3$coefficients["expected_e10", "Estimate"]
colour_3 <- s3$coefficients["expected_e10:Pigmented_extractTRUE", "Estimate"]

# Write functions for herbivorous, frugivorous and carnivorous, with or without pigmentation

fun3_1 <- function(x) {int_3 + (exp_3 + colour_3)*x}
fun3_2 <- function(x) {int_3 + exp_3*x}

p_m3 <-
  p_m3 + 
  geom_text(aes(x = 100, y = 300), alpha=0.5, colour="black", size = 8,
            label = paste0("y = ", round(int_3, 1), " + (", round(exp_3, 2), " ", round(colour_3, 2), "*pigmentation) * x")) +
  geom_function(fun = fun3_1, size = 1, linetype = "dashed", colour = "black") +
  geom_function(fun = fun3_2, size = 1, colour = "black") +
  geom_label(aes(x = 150, y = 250), size=8, colour = "black", alpha=0.5, label = "not pigmented") +
  geom_label(aes(x = 150, y = 150), size=8, colour = "black", alpha=0.5, label = "pigmented extract") +
  labs(tag = "B.")

p_m3

ggsave(filename  =  file.path(outdir, "plot_mmodel3.png"), p_m3, width  =  15, height = 8)

## Box plot with PCR efficiency

# Run kruskal wallis test
# For clear extracts
pval_c <- kruskal_test(efficiency ~ diet_category, data = filter(filt_dat3, Extract == "Clear"))$p
effsize_c <- kruskal_effsize(efficiency ~ diet_category, data = filter(filt_dat3, Extract == "Clear"))$effsize
annot_c <- paste0("Kruskal-Wallis:\neffect size = ", round(effsize_c, 2), ",\np-value = ", round(pval_c, 3))

# For pigmented extracts
pval_p <- kruskal_test(efficiency ~ diet_category, data = filter(filt_dat3, Extract == "Pigmented" & diet_category != "Animalivore"))$p
effsize_p <- kruskal_effsize(efficiency ~ diet_category, data = filter(filt_dat3, Extract == "Pigmented" & diet_category != "Animalivore"))$effsize
annot_p <- paste0("Kruskal-Wallis:\neffect size = ", round(effsize_p, 2), ",\np-value = ", round(pval_p, 2))

ann_text <- data.frame(lab = c(annot_p, annot_c),
                       diet_category = c(0.2, 0.2),
                       efficiency = c(-0.5, -0.5),
                       Extract = c("Pigmented", "Clear"))

model3_box <-
  ggviolin(data = filter(filt_dat3, (Extract == "Clear" | diet_category != "Animalivore")),
           x = "diet_category", y = "efficiency", fill = "diet_category",
           add = "mean_se") +
  scale_fill_manual(values = diet2_palette, name = "") +
  stat_compare_means(comparisons = list(c("Animalivore", "Frugivore"),
                                        c("Animalivore", "Herbivore"),
                                        c("Frugivore", "Herbivore")),
                     label = "p.signif", size = 8) +
  geom_text(data = ann_text, aes(label = lab), size = 6, hjust = 0) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = c(-0.5, 0, 0.5, 1, 1.5, 2, eff)) +
  ylab("Indexing PCR efficiency") + xlab("") +
  facet_grid(cols = vars(Extract)) +
  geom_hline(yintercept = eff, linetype = "dashed") +
  theme_bw_alt + theme(legend.position = "none") +
  labs(tag = "B.")

ggsave(filename  =  file.path(outdir, "model3_boxplot.png"), model3_box, width  =  10, height = 10)

ggsave(filename  =  file.path(outdir, "model2_3_boxplots.png"),
       plot_grid(model2_box, model3_box + xlab("Diet category"), nrow = 2, ncol = 1, align = "v", axis = "lr"),
       width  =  10, height = 16)

# Kruskal-Wallis tests

kruskal.test(data = filt_dat3, efficiency ~ Pigmented_extract)
kruskal_effsize(data = filt_dat3, efficiency ~ Pigmented_extract)

kruskal.test(data = filt_dat3, efficiency ~ diet_category)
kruskal_effsize(data = filt_dat3, efficiency ~ diet_category)

## Consider effect of age (for subset of samples with collection year)
filt_dat3_age <- filt_dat3 %>% filter(!is.na(Age))

# Simple linear model
lm3_age = lm(observed_e10 ~
           expected_e10 + expected_e10:PC1_scaled + expected_e10:PC2_scaled + expected_e10:Pigmented_extract + expected_e10:Age,
         data = filt_dat3_age)

# Extract weights
weights <- 1 / lm(abs(lm3_age$residuals) ~ lm3_age$fitted.values)$fitted.values^2

mm3_age <-  lmer(observed_e10 ~
               expected_e10 + expected_e10:PC1_scaled + expected_e10:PC2_scaled + expected_e10:Pigmented_extract + expected_e10:Age +
               (0 + expected_e10 | Species) + (0 + expected_e10 | Ind.date),
             data = filt_dat3_age, weights = weights)

# Plot diagnostics
diagn_3_age <- diagnose_lmm(mm3_age)

ggsave(filename  =  file.path(outdir, "diagnostics_mmodel3_age.png"), plot_grid(plotlist = diagn_1),
       width  =  14, height = 14)

vif(mm3_age)
shapiro.test(residuals(mm3_age))
hist(residuals(mm3_age))

summary(mm3_age)

#### 4. Weight - Ind numbers (not sure abou this bit) ####

# Remove some outliers
filt_dat4 <- filt_dat3 %>%
  mutate(ind_e10_yield = observed_e10/Sample.weight.mg)

# Plot
p_m4 <-
  ggplot(aes(x = Sample.weight.mg, y = observed_e10, colour = diet_category, shape = Pigmented_extract), data = filt_dat4) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(aes(x = Sample.weight.mg, y = observed_e10), inherit.aes = FALSE) +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  scale_shape_manual(values = c(16, 1), name = "Pigmented extract") +
  scale_x_continuous(name = "Sample weight (mg)") + 
  scale_y_continuous(name = "Indexed molecules (10^10)") +
  theme_minimal() + theme(plot.background  =  element_rect(fill = "white", colour = "white")) 

p_m4

# Simple linear model
lm4 = lm(observed_e10 ~
           Sample.weight.mg + Sample.weight.mg:PC1_scaled + Sample.weight.mg:PC2_scaled + Sample.weight.mg:Pigmented_extract,
         data = filt_dat4)

# Extract weights
weights <- 1 / lm(abs(lm4$residuals) ~ lm4$fitted.values)$fitted.values^2

mm4 <-  lmer(observed_e10 ~
               Sample.weight.mg + Sample.weight.mg:PC1_scaled + Sample.weight.mg:PC2_scaled + Sample.weight.mg:Pigmented_extract +
               (0 + Sample.weight.mg | Species),
             data = filt_dat4, weights = weights)

# Plot diagnostics
diagn_4 <- diagnose_lmm(mm4)

ggsave(filename  =  file.path(outdir, "diagnostics_mmodel4.png"), plot_grid(plotlist = diagn_4),
       width  =  14, height = 14)

vif(mm4)
shapiro.test(residuals(mm4))
hist(residuals(mm4))

s4 <- summary(mm4)
s4

ggsave(filename  =  file.path(outdir, "plot_mmodel4.png"), p_m4, width  =  12, height = 8)

# Box plot with indexed copies per sample
model4_box <-
  ggviolin(data = filter(filt_dat4, (Extract == "Clear" | diet_category != "Animalivore")),
           x = "diet_category", y = "ind_e10_yield", fill = "diet_category",
           add = "mean_se") +
  scale_fill_manual(values = diet2_palette, name = "") +
  stat_compare_means(comparisons = list(c("Animalivore", "Frugivore"),
                                        c("Animalivore", "Herbivore"),
                                        c("Frugivore", "Herbivore")),
                     label = "p.signif") +
  ylab("Indexed library molecules (10^10) per mg of sample") + xlab("") +
  facet_grid(cols = vars(Extract)) +
  geom_hline(yintercept = eff, linetype = "dashed") +
  theme_bw_alt + theme(legend.position = "none")

model4_box

## Consider effect of age (for subset of samples with collection year)
filt_dat4_age <- filt_dat4 %>% filter(!is.na(Age))

# Simple linear model
lm4_age = lm(observed_e10 ~
           Sample.weight.mg + Sample.weight.mg:PC1_scaled + Sample.weight.mg:PC2_scaled + Sample.weight.mg:Pigmented_extract + Sample.weight.mg:Age,
         data = filt_dat4_age)

# Extract weights
weights <- 1 / lm(abs(lm4_age$residuals) ~ lm4_age$fitted.values)$fitted.values^2

mm4_age <-  lmer(observed_e10 ~
               Sample.weight.mg + Sample.weight.mg:PC1_scaled + Sample.weight.mg:PC2_scaled + Sample.weight.mg:Pigmented_extract + Sample.weight.mg:Age +
               (0 + Sample.weight.mg | Species),
             data = filt_dat4_age, weights = weights)

# Plot diagnostics
diagn_4_age <- diagnose_lmm(mm4_age)

ggsave(filename  =  file.path(outdir, "diagnostics_mmodel4_age.png"), plot_grid(plotlist = diagn_1),
       width  =  14, height = 14)

vif(mm4_age)
shapiro.test(residuals(mm4_age))
hist(residuals(mm4_age))

summary(mm4_age)

#### Host vs Microbiome reads ####

filt_dat_h <- filt_dat %>% filter(!is.na(host_count)) %>%
  dplyr::select(Ext.ID, PC1, PC2, Pigmented_extract, Age, Species, Common.name, order, diet_category, host_count, unmapped_count) %>%
  # Calculate host proportion ratio
  mutate(filtered_count = host_count + unmapped_count) %>% # This is the number of reads after preprecessing but before mapping
  mutate(host_perc = host_count / filtered_count) %>%
  # Order
  arrange(order, Species) %>%
  mutate(order = factor(order, level=unique(order))) %>%
  mutate(Species = factor(Species, level=unique(Species))) %>%
  # Shorten order names for labels
  mutate(order_short = case_when(order == "Diprotodontia" ~ "Dipr.",
                                 order == "Proboscidea" ~ "Prob.",
                                 order == "Sirenia" ~ "Sir.",
                                 order == "Carnivora" ~ "Carniv.",
                                 TRUE ~ order))

mean_host <- mean(filt_dat_h$host_perc)

# Plot boxplot
p_h <-
  ggplot(aes(y = host_perc + 0.01, x = Common.name, colour = diet_category), data = filt_dat_h) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(alpha = 0.2, width = 0.2, height = 0.2) +
  geom_hline(yintercept = mean_host, linetype = "dashed") +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  theme_bw_alt + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top",
                       panel.border = element_rect(color = "black", fill = NA)) +
  facet_grid(cols = vars(order_short), scales = "free", space = "free") +
  ylab("Host read content") +
  scale_y_continuous(trans = "log10", labels = percent_format(accuracy = 2), breaks = c(0.01, 0.1, 1, mean_host))

p_h

ggsave(filename  =  file.path(outdir, "host_vs_microbiome.png"), p_h, width  =  16, height = 12)

ggplot(aes(y = host_perc + 0.01, x = order, colour = order), data = filt_dat_h) +
  geom_point(position = position_jitter(0.1), alpha = 0.5, size = 2) +
  geom_boxplot(alpha = 0.5)+
  scale_y_continuous(trans = "log10")

## GLMM Model on order level
gmm_h1 <- glm(cbind(host_count, filtered_count) ~ diet_category + order,
              family = quasibinomial, data = filt_dat_h)

s_h1 <- summary(gmm_h1)
s_h1

car::vif(gmm_h1) # Check for collinearity

# Analysis of deviance
car::Anova(gmm_h1, type = "II") # analysis of deviance

# Rerun without nonsignificant var
gmm_h2 <- glm(cbind(host_count, filtered_count) ~ order,
              family = quasibinomial, data = filt_dat_h)

# GLHT and plot
glht_h2 <- summary(glht(gmm_h2, mcp(order="Tukey")))

df_glht2 <- data.frame(odds_ratio = exp(glht_h2$test$coeff), # Exponentiate to get odds ratio
                       p.value = glht_h2$test$pvalues,
                       comparison = names(glht_h2$test$coeff)) %>%
  # Split comparison column
  separate(comparison, into = c("Order1", "Order2"), sep = " - ") %>%
  arrange(Order2, Order1) %>%
  mutate(label = case_when(p.value < 0.001 ~ "***",
                           p.value < 0.01 ~ "**",
                           p.value < 0.05 ~ "*",
                           p.value < 0.1 ~ "."))

# Save table
write.csv(df_glht2, file = file.path(outdir, "host_order_comparison.csv"), quote = FALSE, row.names = FALSE)

# Plot GLHT results as an output
glht_tiles <- ggplot(data = df_glht2, aes(x = Order1, y = Order2, fill = odds_ratio, label = label)) +
  geom_tile(colour = "black") +
  scale_fill_gradient2(high = "red", low = "blue", midpoint = 1,
                       transform = "log", breaks = c(0.01, 0.1, 1, 5), name = "Ratio") + # midpoint = 1 because odds ratio = 1 means no change
  xlab("1st term of comparison") + ylab("2nd term of comparison") +
  theme_bw_alt + theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1),
                       panel.border = element_rect(color = "black", fill = NA)) +
  scale_y_discrete(position = "right") +
  geom_text(size = 6)

glht_tiles

ggsave(filename  =  file.path(outdir, "tiles_host_vs_microbiome.png"), glht_tiles, width  =  9, height = 7)

## Consider sample age

## GLMM Model on order level
gmm_h_age <- glm(cbind(host_count, filtered_count) ~ order + Age,
              family = quasibinomial, data = filter(filt_dat_h, !is.na(Age)))

summary(gmm_h_age)

car::vif(gmm_h_age) # Check for collinearity

# Analysis of deviance
car::Anova(gmm_h_age, type = "II") # analysis of deviance

#### Source composition ####

## Plot composition by sample

decom_comp <- filt_dat %>% dplyr::select(Ext.ID, Common.name, order, diet_category, Age, starts_with("p_")) %>%
  # Keep only samples with decOM results
  filter(!is.na(p_Unknown)) %>%
  # Remove species with less than 5 samples
  group_by(Common.name) %>% filter(n_distinct(Ext.ID) >= 4) %>% 
  #Calculate sum of ancient and modern oral microbiome 
  mutate(oral_prop = (p_OralH + p_OralTM + p_OralMM),
         contam_prop = (p_Skin + p_Sediment.Soil)) %>%
  ungroup() %>%
  # Group orders with small number of samples
  group_by(order) %>%
  mutate(order_grouped=case_when(n_distinct(Ext.ID) < 25 ~ "Other",
                                 TRUE ~ order)) %>% ungroup %>%
  # Pivot longer
  pivot_longer(starts_with("p_"), names_to = "Source", values_to = "Proportion") %>%
  arrange(desc(Common.name), oral_prop) %>%
  mutate(Common.name = factor(Common.name, levels = unique(Common.name)))

# Turn orders, sources and extraction IDs to factors
ord_levels <- unique(sort(as.character(decom_comp$order_grouped)))
ord_levels <- append(ord_levels[-which(ord_levels == "Other")], "Other")

# Order extraction IDs by species, and get one label per species
id_levels <- unique(dplyr::select(decom_comp, c(Ext.ID, Common.name, oral_prop))) %>%
  mutate(rows = rownames(.)) %>%
  # Get species labels
  group_by(Common.name) %>%
  # Get a species label in the middle of the samples from that species
  mutate(label = case_when(row_number() == round(n_distinct(Ext.ID)/2, digits = 0) ~ Common.name, TRUE ~ "")) %>%
  # Get positions for horizontal lines to separate species
  mutate(line = case_when(row_number() == n_distinct(Ext.ID) ~ Ext.ID, TRUE ~ ""))

new_labels <- id_levels$label
names(new_labels) <- id_levels$Ext.ID

decom_comp <- 
  decom_comp %>% mutate(order_grouped = factor(order_grouped, levels = ord_levels),
                        Source = factor(Source, levels = c("p_OralH", "p_OralTM", "p_OralMM", "p_Rumen", "p_Sediment.Soil", "p_Skin", "p_Unknown")),
                        Ext.ID = factor(Ext.ID, level = id_levels$Ext.ID)) %>%
  arrange(order_grouped, Ext.ID, Source)

# Plot stacked barplot
source_palette <- list("p_OralH"="#EB3838", "p_OralTM"="#A20404", "p_OralMM"="#FF5F5F", "p_Rumen"="#C4D307", "p_Sediment.Soil"="#6E7C15", "p_Skin"="#FFCC73", "p_Unknown"="#AAAAAA")

decom_bar <-
  ggplot(data = decom_comp, aes(y = Ext.ID, x = Proportion, fill = Source, group = Common.name)) +
  geom_bar(stat = "identity", colour = NA) +
  scale_fill_manual(values = source_palette, name = "",
                    labels = c("oral (human)", "oral (terrestrial mammal)", "oral (marine mammal)", "rumen", "sediment/soil", "skin", "unknown")) +
  facet_grid(rows = vars(order_grouped), scales = "free", space = "free", switch = "y") +
  theme_bw_alt + theme(axis.ticks.y = element_blank(),
                       axis.text.y = element_blank(),
                       panel.background = element_rect(fill = NA, colour = NA),
                       panel.border = element_rect(color = "black", fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.position = "top", legend.direction = "vertical") +
  ylab("") + xlab("") +
  # Replace y labels
  scale_y_discrete(labels = new_labels)

# Get diet metadata as a tile plot
diet_meta_plot <-
  ggplot(data = decom_comp, aes(y = Ext.ID, x=1, fill = diet_category, group = Common.name)) +
  geom_tile() + scale_fill_manual(values = diet2_palette, name = "") +
  facet_grid(rows = vars(order_grouped), scales = "free", space = "free", switch = "y") +
    theme_bw_alt + theme(axis.ticks = element_blank(),
                         axis.text.x = element_blank(),
                         panel.background = element_rect(fill = NA, colour = NA),
                         panel.border = element_rect(color = NA, fill = NA),
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.position = "top", legend.direction = "vertical",
                         strip.background = element_blank(),
                         strip.text = element_blank()) + scale_y_discrete(labels = new_labels, position = "right") +
    xlab("") + ylab("")

ggsave(filename  =  file.path(outdir, "decom_barplot.png"),
       plot_grid(decom_bar, diet_meta_plot, ncol = 2, align = "h", rel_widths = c(6, 4)),
       width  =  10, height = 16)

## Plot composition by species
decom_comp_s <- decom_comp %>% group_by(Common.name, order_grouped, Source, diet_category) %>% summarise(Proportion = mean(Proportion)) %>%
  mutate(order_grouped = case_when(order_grouped == "Perissodactyla" ~ "Periss.",
                                   order_grouped == "Carnivora" ~ "Carn.",
                                   TRUE ~ order_grouped)) %>%
  arrange(order_grouped, desc(Common.name), Source)

# Fix order factor
ord_levels <- unique(decom_comp_s$order_grouped)
ord_levels <- append(ord_levels[-which(ord_levels == "Other")], "Other")

decom_comp_s$order_grouped <- factor(decom_comp_s$order_grouped, levels = ord_levels)

decom_bar_s <-
  ggplot(data = decom_comp_s, aes(y = Common.name, x = Proportion, fill = Source, group = Common.name)) +
  geom_bar(stat = "identity", colour = NA) +
  scale_fill_manual(values = source_palette, name = "",
                    labels = c("oral (human)", "oral (terrestrial mammal)", "oral (marine mammal)", "rumen", "sediment/soil", "skin", "unknown")) +
  facet_grid(rows = vars(order_grouped), scales = "free", space = "free", switch = "y") +
  theme_bw_alt + theme(axis.ticks.y = element_blank(),
                       axis.text.y = element_blank(),
                       panel.background = element_rect(fill = NA, colour = NA),
                       panel.border = element_rect(color = "black", fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.position = "bottom", legend.direction = "vertical",
                       legend.location = "plot" ) +
  ylab("")

## Plot oral to contaminant ratio
oral_contam_ratio <- decom_comp %>%
  mutate(ratio = (oral_prop + 0.1)/(contam_prop  + 0.1))

ratio_box <- 
  ggplot(data = oral_contam_ratio, aes(y = Common.name, x = ratio)) +
  geom_boxplot() +
  geom_point() +
  scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10)) +
  facet_grid(rows = vars(order_grouped), scales = "free", space = "free") +
  theme_bw_alt + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       panel.border = element_rect(colour = "black"),
                       strip.background = element_blank(),
                       strip.text = element_blank()) +
  # Add vertical line on 1
  geom_vline(aes(xintercept = 1)) +
  xlab("oral-contaminant ratio") + ylab("")

## Diet metadata
# Get diet metadata as a tile plot
diet_meta_plot <-
  ggplot(data = decom_comp_s, aes(y = Common.name, x=1, fill = diet_category)) +
  geom_tile() + scale_fill_manual(values = diet2_palette, name = "") +
  facet_grid(rows = vars(order_grouped), scales = "free", space = "free", switch = "y") +
  theme_bw_alt + theme(axis.ticks = element_blank(),
                       axis.text.x = element_text(colour = "white"),
                       panel.background = element_rect(fill = NA, colour = NA),
                       panel.border = element_rect(color = NA, fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.position = "bottom", legend.direction = "vertical",
                       strip.background = element_blank(),
                       strip.text = element_blank()) +
  scale_y_discrete(position = "right") +
  xlab("") + ylab("")

## Grid
decom_grid <- plot_grid(decom_bar_s, ratio_box, diet_meta_plot, align = "h", axis = "tb", ncol = 3, rel_widths = c(2, 1, 1))

ggsave(filename  =  file.path(outdir, "decom_species_barplot.png"), decom_grid, width  =  16, height = 12)

## PERMANOVA

# Remove Unknown from data
decom_perm <- decom_comp %>% filter(Source != "p_Unknown") %>%
  # Recalculate proportions
  group_by(Ext.ID) %>% mutate(Proportion = Proportion * 100 / (sum(Proportion))) %>%
  # Remove NAs
  filter(!is.na(Proportion)) %>%
  pivot_wider(names_from = "Source", values_from = "Proportion") %>% ungroup %>%
  # Sum ancient and modern oral
  mutate(p_Oral = p_OralH + p_OralTM + p_OralMM) %>% dplyr::select(-p_OralH, -p_OralTM, -p_OralMM)

# Get factor to be used in the test
composition <- data.frame(dplyr::select(decom_perm, starts_with("p_")))
rownames(composition) <- decom_perm$Ext.ID
order <- factor(decom_perm$order, levels = sort(unique(decom_perm$order)))
diet <- factor(decom_perm$diet_category, levels = sort(unique(decom_perm$diet_category)))

set.seed(123)
perm <-adonis2(composition ~ order + diet,
               method = "euclidean", by = "margin")
perm

# Pairwise comparisons
d <- vegdist(composition, method = "euclidean")

pp_order <- adonis.pair(d, Factor = order, corr.method = "holm")
pp_diet <- adonis.pair(d, Factor = diet, corr.method = "holm")

pp_df <-
  pp_order %>%
  # Split combination column
  separate(combination, into = c("Order1", "Order2"), sep = " <-> ") %>%
  arrange(Order1, Order2) %>%
  mutate(Order1 = factor(Order1, levels = unique(Order1)),
         Order2 = factor(Order2, levels = unique(Order2))) %>%
  # Add an asterisk where p-value is below 0.05
  mutate(label = case_when(P.value.corrected < 0.001 ~ "***",
                           P.value.corrected < 0.01 ~ "**",
                           P.value.corrected < 0.05 ~ "*",
                           P.value.corrected < 0.1 ~ ".",
                           TRUE ~ NA))

# Plot heatmap with pairwise PERMANOVA results
ggplot(data = pp_df, aes(x = Order1, y = Order2, fill = MeanSqs, label = label)) +
  geom_tile(colour = "black") +
  scale_fill_gradient2(high = "red", mid = "orange", low = "white", midpoint = quantile(pp_df$MeanSqs, 0.75),
                       name = "Means Squares") +
  xlab("1st term of comparison") + ylab("2nd term of comparison") +
  theme_bw_alt + theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1),
                       panel.border = element_rect(color = "black", fill = NA)) +
  scale_y_discrete(position = "right") +
  geom_text(size = 6)

## Test differences in oral proportion
# GLM

## GLMM Model on order level
decom_comp$order <- factor(decom_comp$order, ordered = FALSE) # Turn order to factor

gmm_c <- glm(oral_prop/100 ~ order + diet_category,
             family = quasibinomial, data = decom_comp)

s_c <- summary(gmm_c)
s_c

car::vif(gmm_c) # Check for collinearity

# Analysis of deviance
car::Anova(gmm_c, type = "II") # analysis of deviance

# GLHT for orders and plot
glht_c <- summary(glht(gmm_c, mcp(order="Tukey")), test = adjusted(type = "holm"))

df_glht_c <- data.frame(odds_ratio = exp(glht_c$test$coeff), # Exponentiate to get odds ratio
                       p.value = glht_c$test$pvalues,
                       comparison = names(glht_c$test$coeff)) %>%
  # Split comparison column
  separate(comparison, into = c("Order1", "Order2"), sep = " - ") %>%
  arrange(Order2, Order1) %>%
  mutate(label = case_when(p.value < 0.001 ~ "***",
                           p.value < 0.01 ~ "**",
                           p.value < 0.05 ~ "*",
                           p.value < 0.1 ~ "."))

# Save table
write.csv(df_glht_c, file = file.path(outdir, "composition_order_comparison.csv"), quote = FALSE, row.names = FALSE)

# Plot GLHT results as an output
glht_tiles <- ggplot(data = df_glht_c, aes(x = Order1, y = Order2, fill = odds_ratio, label = label)) +
  geom_tile(colour = "black") +
  scale_fill_gradient2(high = "red", low = "blue", midpoint = 1,
                       transform = "log", breaks = c(0.1, 0.3, 1, 3), name = "Ratio") + # midpoint = 1 because odds ratio = 1 means no change
  xlab("1st term of comparison") + ylab("2nd term of comparison") +
  theme_bw_alt + theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1),
                       panel.border = element_rect(color = "black", fill = NA)) +
  scale_y_discrete(position = "right") +
  geom_text(size = 6)

glht_tiles

ggsave(filename  =  file.path(outdir, "tiles_oral_prop.png"), glht_tiles, width  =  9, height = 7)

# GLHT for diets
glht_c2 <- summary(glht(gmm_c, mcp(diet_category="Tukey")))

df_glht_c2 <- data.frame(odds_ratio = exp(glht_c2$test$coeff), # Exponentiate to get odds ratio
                        p.value = glht_c2$test$pvalues,
                        comparison = names(glht_c2$test$coeff))

## Test for age in subset
gmm_c_age <- glm(oral_prop/100 ~ order + diet_category + Age,
             family = quasibinomial, data = filter(decom_comp, !is.na(Age)))

summary(gmm_c_age)

car::vif(gmm_c_age) # Check for collinearity

# Analysis of deviance
car::Anova(gmm_c_age, type = "II") # analysis of deviance

### Repeats ####
# Focus on samples that have been repeated after dilution

# Exclude failed and repeated samples as well as blanks
filt_rep <- data %>% mutate(replicate_pair = str_remove(Ext.ID, "^F_")) %>% group_by(replicate_pair) %>%
  filter(n_distinct(Ext.ID) > 1 & !is.na(Sample.weight.mg)) %>%
  # indicate first and second replicate
  mutate(replicate = case_when(row_number() == 1 ~ "First",
                               row_number() == 2 ~ "Second")) %>% ungroup

filt_rep <- filt_rep %>% group_by(replicate_pair) %>%
  mutate(line_colour = Extract.vol.used[replicate == "Second"]) %>%
  fill(line_colour, .direction = "up") %>%
  # Calculate ratio
  mutate(ratio = LP.copies[replicate == "Second"] / LP.copies[replicate == "First"])

ggplot(data = filt_rep, aes(y = LP.copies, x = replicate, colour = Extract.vol.used, group = replicate_pair)) +
  geom_point(size = 5, alpha = 0.5) +
  scale_colour_gradient(high = "darkblue", low = "turquoise") +
  geom_line(aes(colour = line_colour))

ggplot(data = filt_rep, aes(y = ratio, x = as.character(Extract.vol.used))) +
  geom_point() +
  scale_y_continuous(trans = "log10")
