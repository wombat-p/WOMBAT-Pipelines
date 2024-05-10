## Script to run statistical testing with ROTS on peptide and protein data from WOMBAT-P
## FURTHER FEATURES:  parallelize the comparisons
##                   add normalization methods
source("stat_funcs.R")

############### RUNNING THE SCRIPT ####################
# Read experimental design
exp_design <- read_expdesign("exp_design.txt")

# Read proteins
proteins <- read_proteins("stand_prot_quant.csv")

# Read peptides
peptides <- read_peptides("stand_pep_quant.csv")

# Check if all names are present
all_names <- check_names(exp_design, proteins)
all_names <- check_names(exp_design, peptides)

# Average over technical replicates
proteins <- average_technical_replicates(proteins, exp_design)
peptides <- average_technical_replicates(peptides, exp_design)

# filter for proteins with at least 2 valid values per group
groups <- factor(sub("abundance_", "", sub("_[0-9]*$", "", colnames(proteins))))
proteins <- filter_data(proteins, groups)
peptides <- filter_data(peptides, groups)

# Create pairwise compparisons of all vs all groups
comparisons <- combn(levels(groups), 2)

# run ROTS for each comparison
proteins <- run_rots(proteins, groups, comparisons)
peptides <- run_rots(peptides, groups, comparisons)

write.csv(proteins, "stand_prot_quant_merged.csv", row.names = TRUE)
write.csv(peptides, "stand_pep_quant_merged.csv", row.names = TRUE)

