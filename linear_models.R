##### LINEAR MODELS #####

#### (Generalised) Linear Mixed-effects models

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(here)
library(stringr)
library(ggplot2)
library(ggnewscale)
library(vegan)
library(lme4)
library(lmerTest)
library(blmeco)
library(FDB1)
library(ggfortify)
library(ggpubr)
library(cowplot)
library(multcomp)
library(EcolUtils)

#### VARIABLES AND WORKING DIRECTORY ####
# Make sure working directory is correctly set
wd <- here()
setwd(wd)

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
         Sample.weight.g. = Sample.weight.g. + 0.001,
         Species = factor(Species, levels = sp_levels),
         calculated_species_main_diet = factor(calculated_species_main_diet, levels = c("Herbivore", "Frugivore", "Animalivore", "Undefined")),
         Sample.weight.ug = Sample.weight.g. * 1000)

diet <- diet %>%
  mutate(calculated_species_main_diet = factor(calculated_species_main_diet, levels = c("Herbivore", "Frugivore", "Animalivore", "Undefined")))

# Calculate some variables we need
# E.g. consider extraction volume used for library prep to get input ng of DNA
data <- data %>% mutate(DNA_input_ng = Qubit.concentration.ng.ul. * Extract.vol.used,
                        # and LP vol used for indexing to get input copy number
                        libraries_input_copies = LP.copies * LP.vol.used,
                        # theoretical indexing copies depending on PCR cycles
                        theor_ind_copies = libraries_input_copies * 2^Ind.cycles) %>%
  # failed and repeated samples
  filter(failed.repeated == FALSE)

#Save data
write.table(data, file = file.path(outdir, "data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Exclude  as well as blanks
filt_dat <- data %>% filter(is.blank == FALSE & !is.na(Sample.weight.g.))

#######################
#### Summary plots ####
#######################

## Flowchart

flow <- data.frame(x = rep(c(1, 1, 2), 4),
                   y = c(1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8),
                   type = rep(c("Protocol step", "Quantification", "Value"), 4),
                   group = c(rep("Sample", 3), rep("DNA", 3), rep("Barcoding", 3), rep("Indexing", 3)),
                   name = c("Sampling", "Weighing", "Weight (ug)",
                            "DNA extraction", "Qubit quantification", "DNA concentration (ng/ul)",
                            "Barcoding", "Post-barcoding qPCR","Barcoded libraries (molecules/ul)",
                            "Indexing", "Post-indexing qPCR", "Indexed libraries (molecules/ul)"))

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
  theme(legend.position = "bottom") +
  labs(title = "DNA extraction and library preparation flowchart",
       subtitle = "indicating measurement steps and variables used within models")

flowchart

#Can't get arrows to work -- will add them later
ggsave(filename  =  file.path(outdir, "flowchart.png"), flowchart, width  =  10, height = 8)

#### Diet PCA ####
diet_var <- diet[c("cp", "ee", "cf", "ash", "nfe")]
rownames(diet_var) <- filt_dat$Common.name[match(diet$Species, filt_dat$Species)]

ord <- prcomp(diet_var)

# Get some info for plotting
var_explained <- round(ord$sdev^2 * 100 / sum(ord$sdev^2), 1) # Variance explained

loadings_matrix <- data.frame(Variables = rownames(ord$rotation[,c(1,2)]), ord$rotation[,c("PC1", "PC2")])

labels <- ord$x %>% as.data.frame %>% mutate(label = case_when(PC1 > 10 ~ rownames(.),
                                                               PC2 == max(PC2) | PC2 == min(PC2) ~ rownames(.),
                                                               rownames(.) == "Olive baboon" ~ rownames(.))) %>%
  pull(label)

#labels <- ord$x %>% as.data.frame %>% mutate(label = rownames(.)) %>% pull(label)

# Plot
pca <- ggplot(aes(x = PC1, y = PC2, colour=diet$calculated_species_main_diet), data = data.frame(ord$x)) +
  geom_point(size=3) +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  geom_segment(data = loadings_matrix, aes(x = 0, y = 0, xend = (PC1*18),
                                       yend = (PC2*18)), arrow = arrow(length = unit(0.5, "picas")),
               color = "black") +
  theme_bw_alt + theme(legend.position = "bottom") +
  annotate("text", x = (loadings_matrix$PC1*20), y = (loadings_matrix$PC2*20),
           label = loadings_matrix$Variables, size = 6) +
  xlim(c(-35, 80)) + # Adjust position of text a bit
  labs(tag = "A.")

pca

## Add principal components to data
filt_dat$PC1 <- ord$x[match(filt_dat$Common.name, rownames(ord$x)), 1]
filt_dat$PC2 <- ord$x[match(filt_dat$Common.name, rownames(ord$x)), 2]

# Save table
write.table(filt_dat, file = file.path(outdir, "filtered_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Based on this, combine omnivores and frugivores together
diet <- diet %>% mutate(diet_category = case_when(Species == "Papio anubis" ~ "Frugivore",
                                                          TRUE ~ calculated_species_main_diet))
diet$diet_category <- factor(diet$diet_category, levels = c("Animalivore", "Frugivore", "Herbivore"))

filt_dat <- filt_dat %>% left_join(diet)

# Fix PCA
set.seed(219)
pca <- 
  pca + geom_point(size=3, aes(colour = diet$diet_category)) +
  scale_colour_manual(values = diet2_palette[unique(as.character(diet$diet_category))], name = "Dietary category") +
  geom_text(label = labels, size = 5, hjust = 1,
            position = position_jitter(height = 1), aes(x = PC1 - 2, colour = diet$diet_category))

pca

ggsave(filename  =  file.path(outdir, "dietary_PCA.png"), pca, width  =  8, height = 8)

# Save just dietary data
diet_data <- filt_dat %>% dplyr::select(Species, Common.name, order,
                                 Species_Lintulaakso_et_al,diet_cluster_name, diet_quality,
                                 cp, ee, cf, ash, nfe, calculated_species_main_diet, PC1, PC2) %>%
  unique

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
  theme_void_alt +
  xlab("Dietary category") + ylab("Taxonomic order") +
  labs(tag = "B.") +
  scale_fill_continuous(low = "lightblue", high = "darkblue", name = "Sample number", trans = "log", breaks = c(5, 25, 125)) +
  scale_size(range  =  c(12, 30), breaks = c(min(summ_dat$n_species), mean(unique(summ_dat$n_species)), max(summ_dat$n_species)),
             name = "Species number") +
  new_scale(new_aes = "size") +
  geom_text(aes(label = paste0(n_species, " (", n_samples, ")"), x = diet_category, y = order, colour = lab_col, size = n_species), fontface = "bold") +
  scale_colour_manual(values = summ_dat$lab_col, guide = "none") +
  scale_size_continuous(range = c(3, 6), guide = "none")

p_summary

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
  labs(tag = "A.")

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

cutoff <- 10^8 # Approximate blanks value, cutoff for bad and good samples
bad_samples <- data %>% filter(Ind.copies > 0 & Ind.copies < cutoff & !is.blank) %>% nrow
good_samples <- data %>% filter(Ind.copies > cutoff & !is.blank) %>% nrow
blanks <- data %>% filter(Ind.copies > 0, is.blank) %>% nrow

hist <-
  data %>% ggplot(aes(x = Ind.copies, fill = is.blank, colour = Pigmented_extract)) +
  scale_x_continuous(trans = "log10") +
  geom_histogram(bins = 50, linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "darkgrey", "FALSE" = "#FFCC66"), name = "", labels = c("true sample", "blank")) +
  scale_colour_manual(values = c("TRUE" = "#660000", "FALSE" = "transparent"), name = "", labels = c("Pigmented", "Not pigmented")) +
  xlab("Indexing copies") + ylab("count") +
  geom_vline(xintercept = cutoff, linetype = "dashed") +
  annotate("text", x = 10^6, y = 25, label = paste0("Failed indexing:\n n = ", bad_samples), colour = "#C17D21") +
  annotate("text", x = 10^9, y = 30, label = paste0("Successful indexing:\n n = ", good_samples), colour = "#C17D21") +
  annotate("text", x = 10^8, y = 16, label = paste0("Indexed blanks:\n n = ", blanks)) +
  theme(legend.position = "top") +
  labs(tag = "B.")

hist

ggsave(filename  =  file.path(outdir, "indexing_hist.png"), hist, width  =  8, height = 8)

# Get relationship between blank LP copies and ind. copies
blanks <- data %>% filter(is.blank == TRUE)

ggplot(aes(x = log(theor_ind_copies), y = log(Ind.copies)), data = blanks) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(intercept = 0, slope = 1)

blanks$Ind.copies %>% summary

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
filt_dat1 <- filt_dat %>% filter(Qubit.concentration.ng.ul. < 75)

# Plot
p_m1 <-
  ggplot(aes(x = Sample.weight.ug, y = Qubit.concentration.ng.ul., colour = diet_category, shape=Pigmented_extract), data = filt_dat1) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = diet2_palette, name = "Dietary category") +
  #scale_colour_manual(values = diet_palette, name = "Dietary category") +
  scale_shape_manual(values = c(16, 1), name = "Pigmented extract") +
  xlab("Sample weight (ug)") + ylab("DNA concentration (ng/ul)") +
  theme_bw_alt +
  labs(tag = "A.")

p_m1

# Model
gmm1 <- glmer(Qubit.concentration.ng.ul. ~ Sample.weight.ug + PC1 + PC2 + Pigmented_extract + (0 + Sample.weight.ug|Species) + (0 + Sample.weight.ug|Ext.Date),
              data = filt_dat1, family = Gamma(link="log"))

s1 <- summary(gmm1)

# Output
s1
int_1 <- s1$coefficients["(Intercept)", "Estimate"]
weight_1 <- s1$coefficients["Sample.weight.ug", "Estimate"] # Qubit increases so much for every 1ng of sample added
pc2_1 <- s1$coefficients["PC2", "Estimate"] # Qubit increases so much for every unit of PC2

pc2_val <- data.frame(value=c(10, -10)) #PC2 values to plot
pc2_val$lab <- paste0("PC2 = ", pc2_val$value)

p_m1 <-
  p_m1 + 
  geom_text(aes(x = 20, y = 65), alpha=0.5, colour="black", size = 8,
            label = paste0("y = exp(", round(int_1, 2), " + ", round(weight_1, 2), "* weight + ", round(pc2_1, 2), " * PC2)")) +
  geom_function(fun = function(x) {y = exp(int_1 + x*weight_1 + pc2_val$value[1]*pc2_1)}, size = 1, colour = "black") +
  geom_label(aes(x = 30, y = 35, label = paste(pc2_val$lab[1],"(more carbohydrates)")), colour = "black", size=5, alpha=0.5) +
  geom_function(fun = function(x) {y = exp(int_1 + x*weight_1 + pc2_val$value[2]*pc2_1)}, size = 1, colour = "black") +
  geom_label(aes(x = 30, y = 9, label = paste(pc2_val$lab[2],"(more fibre)")), colour = "black", size=5, alpha=0.5)
  
p_m1

ggsave(filename  =  file.path(outdir, "plot_mmodel1.png"), p_m1, width  =  12, height = 8)

# Plot diagnostics
diagn_1 <- diagnose_glmm(gmm1)
ggsave(filename  =  file.path(outdir, "diagnostics_mmodel1.png"), plot_grid(plotlist = diagn_1),
       width  =  14, height = 14)

## Box plot with average yields

filt_dat1 <- filt_dat1 %>% mutate(DNA_yield = (Qubit.concentration.ng.ul. * 40) / Sample.weight.ug)

model1_box <-
  ggboxplot(data = filt_dat1, x = "diet_category", y = "DNA_yield", fill = "diet_category") +
  scale_fill_manual(values = diet2_palette, name = "") +
  stat_compare_means(comparisons = list(c("Herbivore", "Frugivore"), c("Frugivore", "Animalivore"), c("Herbivore", "Animalivore")),
                     label = "p.signif", method = "wilcox") +
  stat_compare_means(label.y = 120) +
  ylab("DNA yield (ng per ug of sample)") + xlab("Dietary category") +
  theme_bw_alt + theme(legend.position = "none") +
  labs(tag = "B.")

model1_box
ggsave(filename  =  file.path(outdir, "model1_boxplot.png"), model1_box, width  =  14, height = 8)

#### 2. Concentration - LP copy numbers ####

filt_dat2 <- filt_dat %>% filter(!is.na(LP.copies) & DNA_input_ng < 3000)

# Plot
p_m2 <-
  ggplot(aes(x = rank(DNA_input_ng), y = rank(LP.copies), colour = diet_category, shape=Pigmented_extract), data = filt_dat2) +
  geom_point(alpha = 0.7, size = 3) +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  #scale_colour_manual(values = diet_palette, name = "Dietary category") +
  scale_shape_manual(values = c(16, 1), name = "Pigmented extract") +
  scale_x_continuous(name = "rank: DNA input amount (ng)") + 
  scale_y_continuous(name = "rank: Barcoded library molecules (per ul)") +
  theme_bw_alt

p_m2

# Model
mm2 <-  lmer(data = filt_dat2, rank(LP.copies) ~ rank(DNA_input_ng) + PC1 + PC2 + Pigmented_extract + (1|Species) + (1|LP.date))

s2 <- summary(mm2)

# Output
s2
int_2 <- s2$coefficients["(Intercept)", "Estimate"]
dna_2 <- s2$coefficients["rank(DNA_input_ng)", "Estimate"] # Qubit increases so much for every 1ng of sample added
colour_2 <- s2$coefficients["Pigmented_extractTRUE", "Estimate"]

p_m2 <-
  p_m2 + 
  geom_text(aes(x = 300, y = 550), alpha=0.5, colour="black", size = 8,
            label = paste0("y = ", round(int_2, 1), " + ", round(dna_2, 2), " * DNA ", " + (", round(colour_2, 2), ") * (Coloured==TRUE)")) +
  geom_function(fun = function(x) {y = int_2 + x*dna_2}, size = 2, colour = "black") +
  geom_function(fun = function(x) {y = int_2 + x*dna_2 + colour_2}, size = 2, colour = "black", linetype = "dashed") +
  geom_label(aes(x = 45, y = 190, label = "not pigmented"), colour = "black", size=5, alpha=0.5) +
  geom_label(aes(x = 45, y = 120, label = "pigmented extract"), colour = "black", size=5, alpha=0.5)  +
  labs(tag = "A.")

p_m2

ggsave(filename  =  file.path(outdir, "plot_mmodel2.png"), p_m2, width  =  15, height = 8)

# Plot diagnostics
diagn_2 <- diagnose_lmm(mm2)
ggsave(filename  =  file.path(outdir, "diagnostics_mmodel2.png"), plot_grid(plotlist = diagn_2),
       width  =  14, height = 14)

## Box plot with average yields

filt_dat2 <- filt_dat2 %>% mutate(Library_yield = LP.copies * LP.vol.used / DNA_input_ng) %>%
  # For labelling
  mutate(Extract = case_when(Pigmented_extract == TRUE ~ "Pigmented",
                             TRUE ~ "Clear"))
model2_box <-
  ggboxplot(data = filt_dat2, x = "Extract", y = "Library_yield", fill = "diet_category") +
  scale_fill_manual(values = diet2_palette, name = "") +
  stat_compare_means(comparisons = list(c("Clear", "Pigmented")),
                     label = "p.signif") +
  scale_y_continuous(trans = "log10") +
  ylab("Barcoding yield (molecules per ng of DNA)") + xlab("Dietary category") +
  facet_grid(cols = vars(diet_category)) +
  theme_bw_alt + theme(legend.position = "none") +
  labs(tag = "A.")

model2_box

ggsave(filename  =  file.path(outdir, "model2_boxplot.png"), model2_box, width  =  12, height = 8)

# Kruskal-Wallis tests

kruskal.test(data = filt_dat2, Library_yield ~ Pigmented_extract)
kruskal.test(data = filt_dat2, Library_yield ~ diet_category)

#### 3. Theoretic - Real ind numbers ####

filt_dat3 <- filt_dat2 %>% filter(!(is.na(Ind.copies)))

# Plot
p_m3 <-
  ggplot(aes(x = rank(theor_ind_copies), y = rank(Ind.copies), colour = diet_category, shape = Pigmented_extract), data = filt_dat) +
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  #scale_colour_manual(values = diet_palette, name = "Dietary category") +
  scale_shape_manual(values = c(16, 1), name = "Pigmented extract") +
  scale_x_continuous(name = "rank: Expected indexed molecules") + 
  scale_y_continuous(name = "rank: Real indexed molecules") +
  theme_bw_alt

p_m3

# Model
mm3 <-  lmer(data = filt_dat3, rank(Ind.copies) ~ rank(theor_ind_copies) + PC1 + PC2 + Pigmented_extract + (1|Species) + (1|Ind.date))
s3 <- summary(mm3)

# Output
s3
int_3 <- s3$coefficients["(Intercept)", "Estimate"]
theory_3 <- s3$coefficients["rank(theor_ind_copies)", "Estimate"] # Qubit increases so much for every 1ng of sample added
colour_3 <- s3$coefficients["Pigmented_extractTRUE", "Estimate"]

p_m3 <-
  p_m3 + 
  geom_text(aes(x = 300, y = 550), alpha=0.5, colour="black", size = 8,
            label = paste0("y = ", round(int_3, 1), " + libraries * ", round(theory_3, 2), " + (Coloured==TRUE) *", round(colour_3, 2))) +
  geom_function(fun = function(x) {y = int_3 + x*theory_3}, size = 2, colour = "black") +
  geom_function(fun = function(x) {y = int_3 + x*theory_3 + colour_3}, size = 2, colour = "black", linetype = "dashed") +
  geom_label(aes(x = 45, y = 220), size=5, colour = "black", alpha=0.5, label = "not pigmented") +
  geom_label(aes(x = 45, y = 120), size=5, colour = "black", alpha=0.5, label = "pigmented extract") +
  labs(tag = "B.")

p_m3

ggsave(filename  =  file.path(outdir, "plot_mmodel3.png"), p_m3, width  =  15, height = 8)

# Plot diagnostics
diagn_3 <- diagnose_lmm(mm3)
ggsave(filename  =  file.path(outdir, "diagnostics_mmodel3.png"), plot_grid(plotlist = diagn_3),
       width  =  14, height = 14)

## Box plot with average yields

filt_dat3 <- filt_dat3 %>% mutate(real_to_theor = Ind.copies / (theor_ind_copies + 1)) # Adding pseudocount of 1

model3_box <-
  ggboxplot(data = filt_dat3, x = "Extract", y = "real_to_theor", fill = "diet_category") +
  scale_fill_manual(values = diet2_palette, name = "") +
  stat_compare_means(comparisons = list(c("Clear", "Pigmented")),
                     label = "p.signif") +
  scale_y_continuous(trans = "log10") +
  ylab("Indexing efficacy (Observed to Expected molecules)") + xlab("Dietary category") +
  facet_grid(cols = vars(diet_category)) +
  theme_bw_alt + theme(legend.position = "none") +
  labs(tag = "B.")

model3_box

ggsave(filename  =  file.path(outdir, "model3_boxplot.png"), model3_box, width  =  12, height = 8)

# Kruskal-Wallis tests

kruskal.test(data = filt_dat3, real_to_theor ~ Pigmented_extract)
kruskal.test(data = filt_dat3, real_to_theor ~ diet_category)

#### 4. Weight - Ind numbers ####

# Remove some outliers
filt_dat4 <- filt_dat3 %>% filter(Ind.copies < 3*10^11)

# Plot
p_m4 <-
  ggplot(aes(x = Sample.weight.ug, y = Ind.copies, colour = diet_category, shape = Pigmented_extract), data = filt_dat4) +
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  scale_shape_manual(values = c(16, 1), name = "Pigmented extract") +
  scale_x_continuous(name = "Sample weight (ug)") + 
  scale_y_continuous(name = "Indexing copies", trans = "log10") +
  theme_minimal() + theme(plot.background  =  element_rect(fill = "white", colour = "white"))

p_m4

# Model
gmm4 <- glmer(Ind.copies ~ Sample.weight.ug + PC1 + PC2 + Pigmented_extract + (1|Species),
      data = filt_dat3, family = Gamma(link="log"))

s4 <- summary(gmm4)

# Output
s4
int_4 <- s4$coefficients["(Intercept)", "Estimate"]
colour_4 <- s4$coefficients["Pigmented_extractTRUE", "Estimate"]

p_m4 <-
  p_m4 + 
  geom_text(aes(x = 20, y = 1.8*10^11), alpha=0.5, colour="black",
            label = paste0("exp(y = ", round(int_4, 1), " + (Pigmented==TRUE) *", round(colour_4, 2), ")")) +
  geom_function(fun = function(x) {y = exp(int_4)}, size = 2, colour = "black") +
  geom_function(fun = function(x) {y = exp(int_4 + colour_4)}, size = 2, colour = "black", linetype = "dashed") +
  geom_label(aes(x = 30, y = 2*10^10), size=3, colour = "black", alpha=0.5, label = "not pigmented") +
  geom_label(aes(x = 30, y = 2*10^9), size=3, colour = "black", alpha=0.5, label = "pigmented extract") +
  labs(title = paste("Modeling final output of protocol based on sample weight. N =", nrow(filt_dat3)),
       subtitle = "Model lines predicted using GLMM for pigmented and not pigmented extracts")

p_m4

ggsave(filename  =  file.path(outdir, "plot_mmodel4.png"), p_m4, width  =  12, height = 8)

# Plot diagnostics
diagn_4 <- diagnose_glmm(gmm4)
ggsave(filename  =  file.path(outdir, "diagnostics_mmodel4.png"), plot_grid(plotlist = diagn_4),
       width  =  14, height = 14)

# Indexing output per weight bracket


#### Source composition ####

## Plot composition by sample
source_palette <- list("p_aOral"="#E54457", "p_mOral"="#AB0A1D", "p_Sediment.Soil"="#B2AD0B", "p_Skin"="#FFCC73", "p_Unknown"="#AAAAAA")

decom_comp_w <- filt_dat %>% dplyr::select(Ext.ID, Species, order, diet_category, starts_with("p_")) %>%
  # Keep only samples with decOM results
  filter(!is.na(p_Unknown)) %>%
  # Remove species with less than 5 samples
  group_by(Species) %>% filter(n_distinct(Ext.ID) >= 5) %>% ungroup() %>%
  # Group orders with small number of samples
  group_by(order) %>%
  mutate(order_grouped=case_when(n_distinct(Ext.ID) < 20 ~ "Other",
                                 TRUE ~ order)) %>% ungroup

decom_comp <-
  decom_comp_w%>%
  # Pivot longer
  pivot_longer(starts_with("p_"), names_to = "Source", values_to = "Proportion")

# Turn orders, sources and extraction IDs to factors
ord_levels <- unique(decom_comp$order_grouped)
ord_levels <- append(ord_levels[-which(ord_levels == "Other")], "Other")

# Order extraction IDs by species, and get one label per species
id_levels <- unique(dplyr::select(decom_comp, c(Ext.ID, Species))) %>% arrange(Species) %>%
  mutate(rows = rownames(.)) %>%
  # Get species labels
  group_by(Species) %>%
  # Get a species label in the middle of the samples from that species
  mutate(label = case_when(row_number() == round(n_distinct(Ext.ID)/2, digits = 0) ~ Species, TRUE ~ "")) %>%
  # Get positions for horizontal lines to separate species
  mutate(line = case_when(row_number() == n_distinct(Ext.ID) ~ Ext.ID, TRUE ~ ""))

new_labels <- id_levels$label
names(new_labels) <- id_levels$Ext.ID

decom_comp <- 
  decom_comp %>% mutate(order_grouped = factor(order_grouped, levels = ord_levels)) %>%
  mutate(Source = factor(Source, levels = names(source_palette))) %>%
  # Order samples by species
  mutate(Ext.ID = factor(Ext.ID, levels = id_levels$Ext.ID))

# Plot stacked barplot
decom_bar <-
  ggplot(data = decom_comp, aes(y = Ext.ID, x = Proportion, fill = Source, group = Species)) +
  geom_bar(stat = "identity", colour = NA) +
  scale_fill_manual(values = source_palette, name = "",
                    labels = c("ancient oral", "modern oral", "sediment/soil", "skin", "unknown")) +
  facet_grid(rows = vars(order_grouped), scales = "free", space = "free", switch = "y") +
  theme_bw_alt + theme(axis.ticks.y = element_blank(),
                       panel.background = element_rect(fill = NA, colour = NA),
                       panel.border = element_rect(color = "black", fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.position = "top") +
  ylab("Samples") + xlab("") +
  # Add lines to separate species
  geom_hline(data = id_levels, aes(yintercept = line), colour = "white", size = 1) +
  # Replace y labels
  scale_y_discrete(labels = new_labels, position = "right")

decom_bar

ggsave(filename  =  file.path(outdir, "decom_barplot.png"), decom_bar, width  =  8, height = 16)

## Plot composition by species
decom_comp_s <- decom_comp %>% group_by(Species, order_grouped, Source) %>% summarise(Proportion = mean(Proportion)) %>%
  mutate(order_grouped = case_when(order_grouped == "Perissodactyla" ~ "Periss.",
                                   order_grouped == "Carnivora" ~ "Carn.",
                                   TRUE ~ order_grouped))

decom_bar_s <-
  ggplot(data = decom_comp_s, aes(y = Species, x = Proportion, fill = Source)) +
  geom_bar(stat = "identity", colour = NA) +
  scale_fill_manual(values = source_palette, name = "",
                    labels = c("ancient oral", "modern oral", "sediment/soil", "skin", "unknown")) +
  facet_grid(rows = vars(order_grouped), scales = "free", space = "free", switch = "y") +
  theme_bw_alt + theme(axis.ticks.y = element_blank(),
                       panel.background = element_rect(fill = NA, colour = NA),
                       panel.border = element_rect(color = "black", fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.position = "bottom") +
  ylab("") +
  labs(tag = "A.")

decom_bar_s

## Plot oral to contaminant ratio
oral_contam_ratio <- decom_comp_w %>%
  mutate(p_Oral = (p_aOral + p_mOral),
         p_Contam = (p_Skin + p_Sediment.Soil)) %>%
  mutate(ratio = (p_Oral + 1)/(p_Contam + 1))

ratio_box <- 
  ggplot(data = oral_contam_ratio, aes(y = Species, x = ratio)) +
  geom_boxplot() +
  geom_point() +
  scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10)) +
  facet_grid(rows = vars(order_grouped), scales = "free", space = "free") +
  labs(tag = "B.") +
  theme_bw_alt + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       panel.border = element_rect(colour = "black"),
                       strip.background = element_blank(),
                       strip.text = element_blank()) +
  # Add vertical line on 1
  geom_vline(aes(xintercept = 1)) +
  xlab("oral-contaminant\nratio") + ylab("")

ratio_box

## Grid
decom_grid <- plot_grid(decom_bar_s, ratio_box, align = "h", axis = "tb", ncol = 2, rel_widths = c(15, 5))

ggsave(filename  =  file.path(outdir, "decom_species_barplot.png"), decom_grid, width  =  12, height = 12)

## PERMANOVA

# Remove Unknown from data
decom_perm <- decom_comp %>% filter(Source != "p_Unknown") %>%
  # Recalculate proportions
  group_by(Ext.ID) %>% mutate(Proportion = Proportion * 100 / (sum(Proportion))) %>%
  # Remove NAs
  filter(!is.na(Proportion)) %>%
  pivot_wider(names_from = "Source", values_from = "Proportion") %>% ungroup %>%
  # Sum ancient and modern oral
  mutate(p_Oral = p_mOral + p_aOral) %>% dplyr::select(-p_mOral, -p_aOral)

# Get factor to be used in the test
matrix <- dplyr::select(decom_perm, starts_with("p_"))
order <- factor(decom_perm$order, levels = sort(unique(decom_perm$order)))
diet <- factor(decom_perm$diet_category, levels = sort(unique(decom_perm$diet_category)))

perm <-adonis2(matrix ~ order + diet,
               method = "euclidean")
perm

# Pairwise comparisons
d <- vegdist(matrix, method = "euclidean")

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

## Test differences in oral/contam ratio
# GLM

## GLMM Model on order level
oral_contam_ratio$order <- factor(oral_contam_ratio$order) # Turn order to factor

gmm_c <- glm(cbind(p_Oral + 1, p_Contam + 1) ~ order + diet_category,
             family = quasibinomial, data = oral_contam_ratio)

s_c <- summary(gmm_c)
s_c

vif(gmm_c) # Check for collinearity

# Analysis of deviance
car::Anova(gmm_c, type = "II") # analysis of deviance

# GLHT for orders and plot
glht_c <- summary(glht(gmm_c, mcp(order="Tukey")))

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
                       transform = "log", breaks = c(0.3, 1, 3), name = "Ratio") + # midpoint = 1 because odds ratio = 1 means no change
  xlab("1st term of comparison") + ylab("2nd term of comparison") +
  theme_bw_alt + theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1),
                       panel.border = element_rect(color = "black", fill = NA)) +
  scale_y_discrete(position = "right") +
  geom_text(size = 6)

glht_tiles

ggsave(filename  =  file.path(outdir, "tiles_oral_vs_contaminant.png"), glht_tiles, width  =  9, height = 7)

# GLHT for diets
glht_c2 <- summary(glht(gmm_c, mcp(diet_category="Tukey")))

df_glht_c2 <- data.frame(odds_ratio = exp(glht_c2$test$coeff), # Exponentiate to get odds ratio
                        p.value = glht_c2$test$pvalues,
                        comparison = names(glht_c2$test$coeff))

### Repeats ####
# Focus on samples that have been repeated after dilution

# Exclude failed and repeated samples as well as blanks
filt_rep <- data %>% mutate(replicate_pair = str_remove(Ext.ID, "^F_")) %>% group_by(replicate_pair) %>%
  filter(n_distinct(Ext.ID) > 1 & !is.na(Sample.weight.ug)) %>%
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

#### Host vs Microbiome reads ####

filt_dat_h <- filt_dat %>% filter(!is.na(host_count)) %>%
  dplyr::select(Ext.ID, PC1, PC2, Pigmented_extract, Species, order, diet_category, host_count, unmapped_count) %>%
  # Calculate host proportion ratio
  mutate(filtered_count = host_count + unmapped_count) %>% # This is the number of reads after preprecessing but before mapping
  mutate(host_perc = host_count * 100 / filtered_count) %>%
  mutate(host_perc_log = log10(host_perc + 0.01)) %>%
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

# Plot boxplot
p_h <-
  ggplot(aes(y = host_perc + 0.01, x = Species, colour = diet_category), data = filt_dat_h) +
  geom_point(position = position_jitter(0.1), alpha = 0.5, size = 2) +
  geom_boxplot(alpha = 0.5) +
  scale_colour_manual(values = diet2_palette, name = "Dietary category") +
  theme_bw_alt + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top",
                       panel.border = element_rect(color = "black", fill = NA)) +
  facet_grid(cols = vars(order_short), scales = "free", space = "free") +
  ylab("% host reads") +
  scale_y_continuous(trans = "log10")

p_h

ggsave(filename  =  file.path(outdir, "host_vs_microbiome.png"), p_h, width  =  14, height = 12)

ggplot(aes(y = host_perc + 0.01, x = order, colour = order), data = filt_dat_h) +
  geom_point(position = position_jitter(0.1), alpha = 0.5, size = 2) +
  geom_boxplot(alpha = 0.5)+
  scale_y_continuous(trans = "log10")

## GLMM Model on order level
gmm_h1 <- glm(cbind(host_count, filtered_count) ~ diet_category + order,
                family = quasibinomial, data = filt_dat_h)

s_h1 <- summary(gmm_h1)
s_h1

vif(gmm_h1) # Check for collinearity

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
