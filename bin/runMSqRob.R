library(limma)
library(QFeatures)
library(msqrob2)

source("stat_funcs.R")

## reading cmd arguments
out_args <- read_args(args)
normalization_method <- out_args$normalization_method
min_peptides <- out_args$min_peptides


## Reading files
# Experimental design
exp_design <- read_expdesign("exp_design.txt")

peptidesTable <- read_peptides("stand_pep_quant.csv")
proteinsTable <- read_proteins("stand_prot_quant.csv")    

# Check if all names are present
all_names <- check_names(exp_design, peptidesTable)
all_names <- check_names(exp_design, proteinsTable)

# Make QFeatures
peptides <- make_qfeatures(peptidesTable, all_names, exp_design, is_peptides = TRUE)
proteins <- make_qfeatures(proteinsTable, all_names, exp_design, is_peptides = FALSE)

# Normalization
peptides <- normalize_qfeat(peptides, normalization_method, is_peptides=TRUE)
proteins <- normalize_qfeat(proteins, normalization_method, is_peptides=FALSE)

# filter proteins for min_peptides
proteins <- filter_min_peps(proteins, min_peptides)

# Filter peptides for proteins to have at least as many non-NA values as experimental conditions
rowData(peptides)$peptideNorm$numvalues <- rowSums(!is.na(assay(peptides, "peptideNorm")))
rowData(proteins)$proteinNorm$numvalues <- rowSums(!is.na(assay(proteins, "proteinNorm")))
proteins <- filterFeatures(proteins, i="proteinNorm", ~ numvalues >= length(unique(exp_design[,"exp_condition"])))
peptides <- filterFeatures(peptides, i="peptideNorm", ~ numvalues >= length(unique(exp_design[,"exp_condition"])))

# Running MSqRob with error handling
run_msqrob <- function(object, i, formula) {
  tryCatch({
    msqrob(object = object, i = i, formula = formula)
  }, error = function(e) {
    stop("Failed to run msqrob on", i, "with error:", e$message, "\n")
  })
}

if (max(exp_design$techrep) == 2) {
   proteins <- run_msqrob(object = proteins, i = "proteinNorm", formula = ~ (1|exp_condition) + (1|biorep))
   peptides <- run_msqrob(object = peptides, i = "peptideNorm", formula = ~ (1|exp_condition) + (1|biorep))
} else if (max(exp_design$techrep) > 2) {
   proteins <- run_msqrob(object = proteins, i = "proteinNorm", formula = ~ exp_condition + (1|biorep), ridge=T)
   peptides <- run_msqrob(object = peptides, i = "peptideNorm", formula = ~ exp_condition + (1|biorep), ridge=T)
} else {
   proteins <- run_msqrob(object = proteins, i = "proteinNorm", formula = ~ exp_condition)
   peptides <- run_msqrob(object = peptides, i = "peptideNorm", formula = ~ exp_condition)
}


# Now make the contrast matrix
contrasts <- create_contrasts(proteins)

# Test the hypotheses
proteins <- hypothesisTest(object=proteins, i="proteinNorm", contrast=contrasts)
peptides <- hypothesisTest(object=peptides, i="peptideNorm", contrast=contrasts)

## adding the new columns to the data frame
stand_prot_out <- final_assembly(proteins, contrasts, FALSE)
stand_pep_out <- final_assembly(peptides, contrasts, TRUE)

write.csv(stand_prot_out, "stand_prot_quant_merged.csv", row.names = F)
write.csv(stand_pep_out, "stand_pep_quant_merged.csv", row.names = F)
