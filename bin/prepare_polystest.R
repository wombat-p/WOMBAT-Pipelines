#!/usr/bin/Rscript
## Script to prepare data for statistical testing with PolySTest

# Read command line parameters for normalization and threads
args <- commandArgs(trailingOnly = TRUE)
normalization_method <- strsplit(grep("--normalization", args, value = TRUE), split = "=")[[1]][[2]]
threads <- strsplit(grep("--threads", args, value = TRUE), split = "=")[[1]][[2]]


source("stat_funcs.R")

############### RUNNING THE SCRIPT ####################
# Read experimental design
exp_design <- read_expdesign("/data/tmp/exp_design.txt")

# Read proteins
proteins <- read_proteins("/data/tmp/stand_prot_quant.csv")

# Read peptides
peptides <- read_peptides("/data/tmp/stand_pep_quant.csv")

# Check if all names are present
all_names <- check_names(exp_design, proteins)
all_names <- check_names(exp_design, peptides)

# Average over technical replicates
proteins <- average_technical_replicates(proteins, exp_design)
peptides <- average_technical_replicates(peptides, exp_design)

# arrange data for PolySTest
proteins <- arrange_data_for_polystest(proteins)
peptides <- arrange_data_for_polystest(peptides)

write.csv(proteins, "stand_prot_polystest_in.csv", row.names = TRUE)
write.csv(peptides, "stand_pep_polystest_in.csv", row.names = TRUE)

# Write PolySTest yaml file
write_polystest_yaml(proteins, exp_design, TRUE, normalization_method, threads)
write_polystest_yaml(peptides, exp_design, FALSE, normalization_method, threads)
