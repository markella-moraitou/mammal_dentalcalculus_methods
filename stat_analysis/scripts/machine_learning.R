##### MACHINE LEARNING #####

#### ML approach to identify most important factors determining dental calculus processing success

################
#### SET UP ####
################

#### LOAD PACKAGES ####
#remotes::install_github("mlr-org/mlr3extralearners@*release")
#remotes::install_url('https://github.com/catboost/catboost/releases/download/v1.2.5/catboost-R-windows-x86_64-1.2.5.tgz', INSTALL_opts = c("--no-multiarch", "--no-test-load"))
library(tidyr)
library(dplyr)
library(renv)
library(here)
library(ggplot2)
library(mlr3)
library(mlr3tuning)
library(mlr3mbo)
library(mlr3learners)
library(mlr3extralearners)
library(mlr3viz)
library(mlr3tuningspaces)
library(catboost)
library(patchwork)
library(paradox)
library(stringr)
library(pdp)
library(rlist)

#### VARIABLES AND WORKING DIRECTORY ####
# Make sure working directory is correctly set
wd <- here()
setwd(wd)

# Directory paths
indir <- normalizePath(file.path("..","output", "LM"))
outdir <- normalizePath(file.path("..","output", "ML")) # output directory

# Create directory for output
dir.create(outdir, recursive  =  TRUE, showWarnings  =  FALSE)

# Set seed
set.seed(12)

# Theme
source("ggplot_theme.R")
theme_set(theme_bw_alt)

# Load palettes
for (file in dir(file.path("..", "input", "palettes"))) {
  name <- str_remove(file, ".csv")
  palette <- read.csv(file.path("..", "input", "palettes", file)) %>% pull(2)
  names(palette) <-  read.csv(file.path("..", "input", "palettes", file)) %>% pull(1)
  assign(name, palette)
}

#### Function for tuning ####
# Write function for tuning and saving results: this will be used to tune several candidate learners before benchmarking

tune_and_save <- function(learner, label, task) {
  # Define resampling strategy
  set.seed(12)
  rcv = rsmp("repeated_cv", repeats = 10, folds = 10)
  
  # Define terminator
  term = trm("evals", n_evals = 100)
  
  # Define tuner: Bayesian optimisation
  tuner = tnr("mbo")
  
  # Tune learner on task
  instance = TuningInstanceBatchSingleCrit$new(task, learner, rcv, measure, term)
  tuner$optimize(instance)
  instance$result
  
  # Plot tuning results
  p <- wrap_plots(autoplot(instance))
  ggsave(p, filename  =  file.path(outdir, paste0(label, "_tuning.png")), width  =  13, height = 8)
  
  # Add tuning to best_params list and save list
  best_params[[label]] <<- instance$result_learner_param_vals
  saveRDS(best_params, file.path(outdir, "best_params.RDS"))
  list.save(best_params, file.path(outdir, "best_params.yml"))
  # Return best parameters
  return(instance$result_learner_param_vals)
}

#####  Load input data #####

# Filtered data
filt_dat <- read.table(file.path(indir, "filtered_data.tsv"), sep = "\t", quote = "", header = TRUE)

# Filter table
ml_data <- filt_dat %>% dplyr::select(Ext.ID, Species, order, PC1_scaled, PC2_scaled, diet_category, Pigmented_extract, DNA_input_ug, LP.vol.used, ind_output_copies) %>%
  # Define factors
  mutate(order = factor(order),
         Pigmented_extract = factor(Pigmented_extract)) %>%
  filter(!is.na(ind_output_copies)) %>%
  # Get log transformed indexing copies
  mutate(Ind.copies.e11 = ind_output_copies*10^(-11)) %>%
  # Define successfull indexing as above 10^9 copies
  mutate(failed_indexing = (ind_output_copies < 10^9)) %>%
  mutate(ind_output_copies = NULL)

# Make sure order is treated as unordered factor
ml_data$order <- factor(ml_data$order, ordered = FALSE)

# Formula
predictors = " ~ order + Species + PC1_scaled + PC2_scaled + Pigmented_extract + DNA_input_ug + LP.vol.used"

#best_params <- readRDS(file.path(outdir, "best_params.RDS"))

#saved_learners <- dir(path = outdir, pattern = "regr.*.RDS")
#for (filename in saved_learners) {
#  objname <- filename %>% str_remove("regr.") %>% str_remove(".RDS")
#  obj <- readRDS(file = file.path(outdir, filename))
#  assign(objname, obj)
#}

#saved_learners <- dir(path = outdir, pattern = "classif.*.RDS")
#for (filename in saved_learners) {
#  objname <- filename %>% str_remove("classif.") %>% str_remove(".RDS")
#  obj <- readRDS(file = file.path(outdir, filename))
#  assign(objname, obj)
#}

# Create table for renaming variables in plots
relabel <- data.frame("DNA_input_ug" = "DNA input (μg)",
                      "PC1_scaled" = "PC1 (carnivory ), scaled",
                      "PC2_scaled" = "PC2 (herbi- vs frugivory), scaled",
                      "Pigmented_extract" = "Pigmented extract",
                      "order" = "Taxonomic order",
                      "diet_category" = "Diet category",
                      "LP.vol.used" = "Library volume used (μl)")

#########################
#### REGRESSION TASK ####
#########################

# Create task with indexing copies log as the target
task_regr <- as_task_regr(as.formula(paste("Ind.copies.e11", predictors)), data = ml_data)

#### Split training - testing ####
# Set strata (so both train and test datasets get a similar distribution of these factors)
task_regr$set_col_roles("Pigmented_extract", c("feature", "stratum"))
task_regr$set_col_roles("Species", c("stratum"))

set.seed(12)
splits = partition(task_regr, ratio = 0.67)

# Plot datasets to verify stratification
task_regr_split <- task_regr$data()
task_regr_split$Sample <- ml_data$Ext.ID
task_regr_split$dataset <- ifelse(task_regr$row_ids %in% splits$train, "train",
                                  ifelse(task_regr$row_ids %in% splits$test, "test", NA))


ggplot(aes(x = dataset, fill = order), data = task_regr_split) + geom_bar()

# Save
write.csv(task_regr_split, file = file.path(outdir, "task_split.csv"), quote = FALSE, row.names = FALSE)

# Use only the training data for tuning
task_train = task_regr$clone()$filter(splits$train)
task_test = task_regr$clone()$filter(splits$test)

# Plot task features
task_regr_plot <- autoplot(task_regr, type = "pairs")
ggsave(task_regr_plot, filename  =  file.path(outdir, "regression_task_data.png"), width  =  12, height = 12)

#### Prepare for tuning and benchmarking ####

# Use Mean Absolute Error as measure
measure = msr("regr.mae")

# Best parameters
best_params = list()

#### Featureless and lm ####
lrn_featureless = lrn("regr.featureless")
lrn_lm = lrn("regr.lm")

#### Tune catboost ####
#lrn_catboost = lrn("regr.catboost",
#                   iterations = to_tune(100,500),
#                   learning_rate = to_tune(0.01, 0.3),
#                   depth = to_tune(4, 10),
#                   l2_leaf_reg = to_tune(1, 10),
#                   bagging_temperature = to_tune(0, 1),
#                   #subsample = to_tune(0.5, 1),
#                   random_strength = to_tune(0, 1)
#                   )

#lrn_catboost$param_set$values <- tune_and_save(lrn_catboost, "regr.catboost", task_train)

#### Tune kknn ####
lrn_kknn = lrn("regr.kknn",
               k = to_tune(1, 30),
               distance = to_tune(0, 5),
               kernel = to_tune(c("rectangular", "triangular", "epanechnikov")))

lrn_kknn$param_set$values <- tune_and_save(lrn_kknn, "regr.kknn", task_train)

#### Tune nnet ####
lrn_nnet = lrn("regr.nnet",
               size = to_tune(1, 20),
               decay = to_tune(-10, 10),
               maxit = to_tune(100, 1000))

lrn_nnet$param_set$values <- tune_and_save(lrn_nnet, "regr.nnet", task_train)

#### Tune randomForest ####
lrn_randomForest = lrn("regr.randomForest", 
                       ntree = to_tune(100, 1000), 
                       mtry = to_tune(1, task_regr$n_features),
                       nodesize = to_tune(5, 30), 
                       maxnodes = to_tune(10, 30), 
                       replace = FALSE
                       )

lrn_randomForest$param_set$values <- tune_and_save(lrn_randomForest, "regr.randomForest", task_train)

#### Tune ranger ####
lrn_ranger = lrn("regr.ranger", 
                 num.trees = to_tune(100, 1000), 
                 mtry = to_tune(1, task_regr$n_features), 
                 min.node.size = to_tune(5, 30), 
                 sample.fraction = to_tune(0.5, 1), 
                 max.depth = to_tune(4, 20), 
                 splitrule = to_tune(c("variance", "extratrees", "maxstat")),
                 importance = "impurity")

# Tune on training dataset
lrn_ranger$param_set$values <- tune_and_save(lrn_ranger, "regr.ranger", task_train)

#### Benchmark ####
rcv = rsmp("repeated_cv", repeats = 10, folds = 10)

# Get benchmark design: Benchmark using the training set
learners <- ls(pattern = "^lrn_.*")

design = benchmark_grid(
  tasks = task_train,
  learners = lapply(learners, get),
  resamplings = rcv
  )

# Save plot
set.seed(12)
bmr = benchmark(design)
p <- autoplot(bmr, measure = measure) +
  theme(theme_bw_alt) +
  scale_fill_viridis_d() +
  theme(plot.background = element_rect(fill = "white")) +
  theme(legend.location = "none")

ggsave(p, filename  =  file.path(outdir, "benchmarking_regr.png"), width  =  13, height = 8)

# Save aggregate scores
bmr_scores <- bmr$aggregate(measure)
write.csv(dplyr::select(bmr_scores, c(nr, learner_id, measure$id)), file = file.path(outdir, "bmr_regr_aggregate_scores.csv"), quote = FALSE, row.names = FALSE)

#### Evaluate using test dataset ####
evaluations <- list()

for (name in learners) {
  set.seed(12)
  obj <- get(name)
  obj$train(task_train)
  evaluations[paste(name, "test", sep = "_")] <- obj$predict(task_test)$score(measure)
  evaluations[paste(name, "train", sep = "_")] <- obj$predict(task_train)$score(measure)
  p <- obj$predict(task_test) %>% autoplot
  p <- p + theme(plot.background = element_rect(fill = "white"))
  ggsave(p, filename  =  file.path(outdir, paste0("regr.", name, "_predictions.png")), width  =  13, height = 8)
}

# Plot
evaluations_df <- unlist(evaluations) %>% as.data.frame()
colnames(evaluations_df) <- measure$label
evaluations_df$learner <- rownames(evaluations_df) %>% str_remove(., "lrn_") %>% paste0("regr.", .) %>% str_remove("_.*")
evaluations_df$dataset <- rownames(evaluations_df) %>% str_remove(".*_")

write.csv(evaluations_df, file = file.path(outdir, "model_evaluation_regr.csv"), quote = FALSE, row.names = FALSE)

p <- ggplot(aes(y = !!sym(measure$label), x = dataset, fill = learner), data = evaluations_df) +
  facet_grid(cols = vars(learner)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() + theme(legend.position = "none") 

ggsave(p, filename  =  file.path(outdir, "model_evaluation_regr.png"), width  =  13, height = 8)

## Save all learners
for (l in learners) {
  name = paste0("regr.", l, ".RDS")
  obj = get(l)
  saveRDS(object = obj, file = file.path(outdir, name))
}

#### Extract best model ####
final_regr <- lrn_ranger
model <- final_regr$model

## Variable importance
importance <- ranger::importance(model)
importance_df <- data.frame(
  Variable = names(importance),
  Importance = as.numeric(importance)
) %>% arrange(-Importance)

write.csv(importance_df, file = file.path(outdir, "regr.final_learner_importance.csv"), quote = FALSE, row.names = FALSE)

# Plot
importance_df$Variable <- factor(importance_df$Variable, levels = rev(importance_df$Variable))
pImp <- ggplot(data = importance_df, aes(x = Importance, y = Variable)) + geom_bar(stat = "identity", fill = "#AA3C39") +
  scale_y_discrete(labels = relabel)

## Partial dependence plots

# Get partial dependence for each variable and overlay with real data

for (i in 1:length(task_regr$feature_names)) {
  varname <- task_regr$feature_names[i]
  vartype <- task_regr$feature_types[i]$type
  # Get partial dependence
  pd <- partial(model, pred.var = varname, train = task_regr$data())

  # Plot
  p <- ggplot(data = pd, aes(y = yhat, x = .data[[varname]])) +
    ylab("Predicted indexed molecules (log10)") +
    xlab(relabel[varname])
  if (vartype == "numeric" | vartype == "integer") { p <- p + geom_line() }
  if (vartype == "factor") { p <- p + geom_bar(stat = "identity") }
  assign(paste0("p", i), p)
}

p_grid <- cowplot::plot_grid(pImp + labs(tag = "A."), 
                   p1 + ylab("Prediction") + labs(tag = "B."),
                   p4 + ylab("Prediction") + labs(tag = "C."),
                   p2 + ylab("Prediction") + labs(tag = "D."))

ggsave(p_grid, filename  =  file.path(outdir, "regr_final_learner_pd.png"), width  =  18, height = 8)

#############################
#### CLASSIFICATION TASK ####
#############################

# Create task with indexing copies log as the target
task_classif <- as_task_classif(as.formula(paste("failed_indexing", predictors)), data = ml_data)
task_classif$set_col_roles("Species", c("stratum"))

# Plot task features
task_classif_plot <- autoplot(task_classif, type = "pairs")

ggsave(task_classif_plot, filename  =  file.path(outdir, "classification_task_data.png"), width  =  12, height = 12)

#### Split training - testing ####

# Use only the training data for tuning
task_train = task_classif$clone()$filter(splits$train)
task_test = task_classif$clone()$filter(splits$test)

#### Prepare for tuning and benchmarking ####

# Remove previous learners
rm(list = ls(pattern = "lrn_.*"))

# Classification measure
measure = msr("classif.logloss")

#### Featureless ####
lrn_featureless = lrn("classif.featureless", predict_type = "prob")

#### Tune catboost ####
#lrn_catboost = lrn("classif.catboost", predict_type = "prob",
#                   iterations = to_tune(100,500),
#                   learning_rate = to_tune(0.01, 0.3),
#                   depth = to_tune(4, 10),
#                   l2_leaf_reg = to_tune(1, 10),
#                   bagging_temperature = to_tune(0, 1),
#                   #subsample = to_tune(0.5, 1),
#                   random_strength = to_tune(0, 1))

#lrn_catboost$param_set$values <- tune_and_save(lrn_catboost, "classif.catboost", task_train)

#### Tune kknn ####
lrn_kknn = lrn("classif.kknn", predict_type = "prob",
               k = to_tune(1, 30),
               distance = to_tune(0, 5),
               kernel = to_tune(c("rectangular", "triangular", "epanechnikov")))

lrn_kknn$param_set$values <- tune_and_save(lrn_kknn, "classif.kknn", task_train)

#### Tune nnet ####
lrn_nnet = lrn("classif.nnet", predict_type = "prob",
               size = to_tune(1, 20),
               decay = to_tune(-10, 10),
               maxit = to_tune(100, 1000))

lrn_nnet$param_set$values <- tune_and_save(lrn_nnet, "classif.nnet", task_train)

#### Tune randomForest ####
#lrn_randomForest = lrn("classif.randomForest", predict_type = "prob", 
#                       ntree = to_tune(100, 1000), 
#                       mtry = to_tune(1, task_classif$n_features),
#                       nodesize = to_tune(5, 30), 
#                       maxnodes = to_tune(1, 30), 
#                       replace = FALSE
#                       )

#lrn_randomForest$param_set$values <- tune_and_save(lrn_randomForest, "classif.randomForest", task_train)

#### Tune ranger ####
lrn_ranger = lrn("classif.ranger", predict_type = "prob", 
                 num.trees = to_tune(100, 1000), 
                 mtry = to_tune(1, task_classif$n_features), 
                 min.node.size = to_tune(5, 30), 
                 sample.fraction = to_tune(0.5, 1), 
                 max.depth = to_tune(4, 25), 
                 splitrule = to_tune(c("gini", "extratrees", "hellinger")),
                 importance = "impurity")

# Tune on training dataset
lrn_ranger$param_set$values <- tune_and_save(lrn_ranger, "classif.ranger", task_train)

#### Tune rpart ####
lrn_rpart = lrn("classif.rpart", predict_type = "prob",
                cp = to_tune(0.001, 0.1),  # Complexity parameter
                minbucket = to_tune(5, 30),  # Minimum number of observations in any terminal node
                maxdepth = to_tune(1, 30)  # Maximum depth of any node of the final tree
)

lrn_rpart$param_set$values <- tune_and_save(lrn_rpart, "classif.rpart", task_train)

#### Benchmark ####
rcv = rsmp("repeated_cv", repeats = 10, folds = 10)

# Get benchmark design: Benchmark using the training set
learners <- ls(pattern = "^lrn_.*")

design = benchmark_grid(
  tasks = task_train,
  learners = lapply(learners, get),
  resamplings = rcv
)

# Save plot
set.seed(12)
bmr = benchmark(design)
p <- autoplot(bmr, measure = measure) +
  theme(theme_bw_alt) +
  scale_fill_viridis_d() +
  theme(plot.background = element_rect(fill = "white")) +
  theme(legend.location = "none")

ggsave(p, filename  =  file.path(outdir, "benchmarking_classif.png"), width  =  13, height = 8)

# Save aggregate scores
bmr_scores <- bmr$aggregate(measure)
write.csv(dplyr::select(bmr_scores, c(nr, learner_id, measure$id)), file = file.path(outdir, "bmr_classif_aggregate_scores.csv"), quote = FALSE, row.names = FALSE)

#### Evaluate using test dataset ####
evaluations <- list()

for (name in learners) {
  set.seed(12)
  obj <- get(name)
  obj$train(task_train)
  evaluations[paste(name, "test", sep = "_")] <- obj$predict(task_test)$score(measure)
  evaluations[paste(name, "train", sep = "_")] <- obj$predict(task_train)$score(measure)
  # Convert the confusion matrix to a data frame for ggplot
  confusion_df <- as.data.frame(obj$predict(task_test)$confusion)
  
  # Plotting the confusion matrix
  p <- 
    ggplot(confusion_df, aes(x = response, y = truth, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "white") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(title = "Confusion Matrix", x = "Predicted", y = "Actual")

  ggsave(p, filename  =  file.path(outdir, paste0("classif.", name, "_predictions.png")), width  =  13, height = 8)
}

# Plot
evaluations_df <- unlist(evaluations) %>% as.data.frame()
colnames(evaluations_df) <- measure$label
evaluations_df$learner <- rownames(evaluations_df) %>% str_remove(., "lrn_") %>% paste0("classif.", .) %>% str_remove("_.*")
evaluations_df$dataset <- rownames(evaluations_df) %>% str_remove(".*_")

write.csv(evaluations_df, file = file.path(outdir, "model_evaluation_classif.csv"), quote = FALSE, row.names = FALSE)

p <- ggplot(aes(y = !!sym(measure$label), x = dataset, fill = learner), data = evaluations_df) +
  facet_grid(cols = vars(learner)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() + theme(legend.position = "none") 

ggsave(p, filename  =  file.path(outdir, "model_evaluation_classif.png"), width  =  13, height = 8)

## Save all learners

for (l in learners) {
  name = paste0("classif.", l, ".RDS")
  obj = get(l)
  saveRDS(object = obj, file = file.path(outdir, name))
}

#### Extract best model ####
final_learner <- lrn_ranger
model <- final_learner$model

# How good is the model?
predictions = final_learner$predict(task_test)
predictions$confusion
predictions$score(msr("classif.fpr"))
predictions$score(msr("classif.fnr"))

## Examine false postive rate
# Get predictions and label them as false positive and negative
predictions_df <- data.frame(row_ids = predictions$row_ids,
                             truth = predictions$truth,
                             response = predictions$response,
                             prob = predictions$prob[,"TRUE"]) %>%
  mutate(Prediction = case_when(truth == response & response == TRUE ~ "TP",
                              truth == response & response == FALSE ~ "TN",
                              truth != response & response == TRUE ~ "FP",
                              truth != response & response == FALSE ~ "FN"))

# Add indexing output (copies)
predictions_df <- cbind(predictions_df, ml_data[predictions_df$row_ids,])
 
ggplot(data = predictions_df, aes(x = Ind.copies.e11, fill = Prediction)) +
  geom_histogram() +
  scale_x_continuous(trans = "log10")

ggplot(data = predictions_df, aes(x = Ind.copies.e11, y = prob, colour = response, group = response)) +
         geom_point(position = position_jitter()) +
  ylab("Predicted probability of failing indexing") +
  scale_x_continuous(trans="log10") +
  geom_vline(xintercept = 10^-2) + geom_hline(yintercept = 0.5)

## Variable importance
importance <- model$variable.importance
importance_df <- data.frame(
  Variable = names(importance),
  Importance = as.numeric(importance)
) %>% arrange(-Importance)

write.csv(importance_df, file = file.path(outdir, "classif_final_learner_importance.csv"), quote = FALSE, row.names = FALSE)

# Plot
importance_df$Variable <- factor(importance_df$Variable, levels = rev(importance_df$Variable))
pImp <- ggplot(data = importance_df, aes(x = Importance, y = Variable)) + geom_bar(stat = "identity", fill = "#AA3C39") +
  scale_y_discrete(labels = relabel)

## Partial dependence plots

for (i in 1:length(task_classif$feature_names)) {
  varname <- task_classif$feature_names[i]
  vartype <- task_classif$feature_types[i]$type
  # Get partial dependence
  pd <- partial(model, pred.var = varname, train = task_classif$data())
  
  # Plot
  p <- ggplot(data = pd, aes(y = yhat, x = !!sym(varname))) +
    ylab("Probability of failed indexing") +
    xlab(relabel[varname])
  
  if (vartype == "numeric" | vartype == "integer") { p <- p + geom_line() }
  if (vartype == "factor") { p <- p + geom_bar(stat = "identity") }
  assign(paste0("p", i), p)
}

p_grid <- cowplot::plot_grid(pImp + labs(tag = "A."), 
                             p1 + labs(tag = "B."),
                             p3 + labs(tag = "C."),
                             p4 + labs(tag = "D."))

ggsave(p_grid, filename  =  file.path(outdir, "classif_final_learner_pd.png"), width  =  18, height = 8)

## Get predictions

# Create hypothetical data for one artiodactyl, one frugivorous primate and one carnivoran
hypothetical_species <- data.frame(species = c("hypothetical species 1", "hypothetical species 2", "hypothetical species 3"),
                                   order = c("Artiodactyla", "Primates", "Carnivora"),
                                   PC1_scaled = c(-0.5, -0.5, 2.5),
                                   PC2_scaled = c(-1, 1, 0),
                                   diet_category = c("Herbivore", "Frugivore", "Animalivore"))

bracket_size = 0.05 # will also use this later for the real data
max = 0.5
concs <- seq(0, max, bracket_size)
sp <- hypothetical_species$species
pigm <- factor(c(TRUE, FALSE), levels = c("FALSE", "TRUE"))

# Get all combinations of data
new_data <- expand.grid(DNA_input_ug = concs,
                        species = sp,
                        LP.vol.used = 18,
                        Pigmented_extract = pigm) %>%
  left_join(hypothetical_species) %>%
  # Remove pigmented carnivores: very uncommon
  filter(!(order == "Carnivora" & Pigmented_extract == TRUE))

# Get prediction
new_data$pred.success <- final_learner$predict_newdata(new_data)$data$prob[, "FALSE"]*100

# Get real data for comparison
# Keep only samples that roughlt correspond to these hypothetical taxa
real_data <- ml_data %>% filter((order == "Carnivora" & diet_category == "Animalivore") |
                                  (order == "Primates" & diet_category == "Frugivore") |
                                  (order == "Artiodactyla" & diet_category == "Herbivore")) %>%
  filter(DNA_input_ug < max+bracket_size/2) %>%
  # create brackets for the DNA input variable
  mutate(binned_input = round(DNA_input_ug/bracket_size, 0)*bracket_size) %>% group_by(binned_input, order, Pigmented_extract, diet_category) %>%
  summarise(total = n_distinct(Ext.ID), success = n_distinct(Ext.ID[failed_indexing==FALSE])) %>%
  mutate(real.success = round(success * 100 / total, digits = 1)) %>%
  # Add diet labels
  left_join(unique(select(new_data, c(order, diet_category))))

combined_data <- left_join(new_data, real_data, by = c("DNA_input_ug" = "binned_input",
                                                       "Pigmented_extract" = "Pigmented_extract",
                                                       "diet_category" = "diet_category",
                                                       "order" = "order"))

write.csv(combined_data, file = file.path(outdir, "classif_predictions_hypothetical_data.csv"), quote = FALSE, row.names = FALSE)

# Plot
labels_vec <- c("Animalivore" = "animalivorous carnivoran",
  "Frugivore" = "frugivorous primate",
  "Herbivore" = "herbivorous artiodactyl",
  "Omnivore" = "NA")

p <-
  ggplot(data = combined_data, aes(y = pred.success, x = DNA_input_ug)) +
  # Add predicted values for hypothetical data
  geom_line(aes(colour = diet_category, linetype = Pigmented_extract), linewidth = 1.5) +
  scale_linetype(name = "Hypothetical samples", labels = c(`FALSE` = "Clear extract", `TRUE` = "Pigmented extract")) +
  scale_colour_manual(values = diet2_palette, name = "Host type", labels = labels_vec) +
  # Add real data
  ggnewscale::new_scale_colour() +
  geom_point(aes(y = real.success, shape = Pigmented_extract, colour = diet_category),
             size = 5, alpha = 0.7, position = position_jitter(width = 0.01, height = 0.02)) +
  scale_shape_manual(values = c(16, 8), name = "Real samples", labels = c(`FALSE` = "Clear extract", `TRUE` = "Pigmented extract")) +
  scale_colour_manual(values = diet2_palette, name = "Host type", labels = labels_vec) +
  # Other
  ylab("Probability of indexing success (%)") + xlab("DNA input (μg)") +
  geom_label(x = 100, y = 45, label = "Pigmented extracts", size=5, alpha=0.5, colour = "black") +
  geom_label(x = 100, y = 80, label = "Clear extracts", size=5, alpha=0.5, colour = "black") +
  labs(tag = "C.") +
  theme(
    legend.position = "top",
    legend.direction = "vertical",
    legend.text = element_text(size = 20)
  ) +
  guides(fill = guide_legend(ncol = 1))
p

ggsave(p, filename  =  file.path(outdir, "classif_predictions_hypothetical.png"), width  =  11, height = 16)
