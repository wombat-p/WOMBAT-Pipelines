library(limma)
library(QFeatures)
library(msqrob2)

## reading cmd arguments
args <- commandArgs(trailingOnly = TRUE)
normalization_method <- strsplit(grep("--normalization", args, value = TRUE), split = "=")[[1]][[2]]
min_peptides <- strsplit(grep("--min_peptides", args, value = TRUE), split = "=")[[1]][[2]]
if (!any(normalization_method == c("sum", "median", "mean", "quantiles", "none"))) {
  stop("Invalid normalization method, should be one of: sum, median, mean, quantiles, none")
}
if (any(normalization_method == c("mean", "median"))) {
  normalization_method <- paste0("center.", normalization_method)
}
if (any(normalization_method == c("quantiles"))) {
  normalization_method <- "quantiles.robust"
}
if (min_peptides < 1) {
  stop("Minimum number of peptides should be at least 1")
}

## Reading files
# Experimental design
exp_annotation <- read.csv("exp_design.txt", sep = "\t")
# Check for existing column names
for (col in c("raw_file", "exp_condition", "biorep", "techrep", "fraction")) {
  if (!any(col %in% colnames(exp_annotation))) {
    stop(paste("Missing column in experimental design file:", col))
  }
}

# Remove all fractions larger than one as they have been summed into one sample
exp_annotation <- exp_annotation[exp_annotation$fraction == 1, ]
exp_annotation$raw_file <- tools::file_path_sans_ext(exp_annotation$raw_file)
exp_annotation$exp_condition <- make.names(exp_annotation$exp_condition)
# needed to ensure factors
exp_annotation$biorep <- as.character(exp_annotation$biorep)
exp_annotation$run <- paste0("abundance_", exp_annotation$exp_condition, "_", exp_annotation$biorep)


# Peptides
peptidesTable <- read.csv("stand_pep_quant.csv")

ecols <- grep(
  "^abundance_",
  names(peptidesTable)
)

peptides <- readQFeatures(
  table = peptidesTable,
  fnames = "modified_peptide",
  ecol = ecols,
  name = "peptideRaw", sep = ","
)

colData(peptides)$exp_condition[exp_annotation$run] <- exp_annotation$exp_condition
colData(peptides)$biorep[exp_annotation$run] <- exp_annotation$biorep
colData(peptides)$exp_condition <- factor(colData(peptides)$exp_condition)
colData(peptides)$biorep <- factor(colData(peptides)$biorep
peptides <- zeroIsNA(peptides, "peptideRaw")
peptides <- logTransform(peptides, base = 2, i = "peptideRaw", name = "peptideLog")

# Proteins
proteinTable <- read.csv("stand_prot_quant.csv")

ecols <- grep(
  "^abundance_",
  names(proteinTable)
)
proteins <- readQFeatures(
  table = proteinTable,
  fnames = "protein_group",
  ecol = ecols,
  name = "proteinRaw", sep = ","
)

colData(proteins)$exp_condition[exp_annotation$run] <- exp_annotation$exp_condition
colData(proteins)$biorep[exp_annotation$run] <- exp_annotation$biorep
colData(proteins)$exp_condition <- factor(colData(proteins)$exp_condition)
colData(proteins)$biorep <- factor(colData(proteins)$biorep)

## Normalization
# take exp and log for sum only
if (normalization_method == "none") {
  peptides <- addAssay(peptides, peptides[["peptideLog"]],
    name = "peptideNorm"
  )
  addAssayLinkOneToOne(peptides, "peptideLog", "peptideNorm")
  proteins <- addAssay(proteins, proteins[["proteinRaw"]],
    name = "proteinNorm"
  )
  addAssayLinkOneToOne(proteins, "proteinRaw", "proteinNorm")
} else {
  if (normalization_method == "sum") {
    ttt <- peptides[["peptideLog"]]
    assay(ttt) <- 2^assay(ttt)
    peptides[["peptideLog"]] <- ttt
    ttt <- proteins[["proteinRaw"]]
    assay(ttt) <- 2^assay(ttt)
    proteins[["proteinRaw"]] <- ttt
  }
  peptides <- normalize(peptides,
    i = "peptideLog",
    name = "peptideNorm",
    method = normalization_method
  )
  proteins <- normalize(proteins,
    i = "proteinRaw",
    name = "proteinNorm",
    method = normalization_method
  )
  if (normalization_method == "sum") {
    ttt <- peptides[["peptideNorm"]]
    assay(ttt) <- log2(assay(ttt))
    peptides[["peptideNorm"]] <- ttt
    ttt <- proteins[["proteinNorm"]]
    assay(ttt) <- log2(assay(ttt))
    proteins[["proteinNorm"]] <- ttt
  }
}


# filter proteins for min_peptides
prots <- proteins[["proteinNorm"]]
nums <- rowData(proteins)[["proteinRaw"]]
nums <- as.matrix(nums[,grep("^number_of_peptides_", names(nums))])
assay(prots)[nums < min_peptides ] <- NA

# Filter peptides for proteins to have at least as many non-NA values as experimental conditions
rowData(peptides)$peptideNorm$numvalues <- rowSums(!is.na(assay(peptides, "peptideNorm")))
rowData(proteins)$proteinNorm$numvalues <- rowSums(!is.na(assay(proteins, "proteinNorm")))
proteins <- filterFeatures(proteins, i="proteinNorm", ~ numvalues >= length(unique(exp_annotation[,"exp_condition"])))
peptides <- filterFeatures(peptides, i="peptideNorm", ~ numvalues >= length(unique(exp_annotation[,"exp_condition"])))

# Running MSqRob with error handling
run_msqrob <- function(object, i, formula) {
  tryCatch({
    msqrob(object = object, i = i, formula = formula)
  }, error = function(e) {
    stop("Failed to run msqrob on", i, "with error:", e$message, "\n")
  })
}

if (max(exp_annotation$techrep) == 2) {
   proteins <- run_msqrob(object = proteins, i = "proteinNorm", formula = ~ (1|exp_condition) + (1|biorep))
   peptides <- run_msqrob(object = peptides, i = "peptideNorm", formula = ~ (1|exp_condition) + (1|biorep))
} else if (max(exp_annotation$techrep) > 2) {
   proteins <- run_msqrob(object = proteins, i = "proteinNorm", formula = ~ exp_condition + (1|biorep), ridge=T)
   peptides <- run_msqrob(object = peptides, i = "peptideNorm", formula = ~ exp_condition + (1|biorep), ridge=T)
} else {
   proteins <- run_msqrob(object = proteins, i = "proteinNorm", formula = ~ exp_condition)
   peptides <- run_msqrob(object = peptides, i = "peptideNorm", formula = ~ exp_condition)
}


# Now make the contrast matrix
conditions <- levels(colData(proteins)$exp_condition)

# Create a design matrix
design <- model.matrix(~ 0 + colData(proteins)$exp_condition)
colnames(design) <- conditions

# Choose the type of contrast: all-vs-all (TODO, use optional parameter for all-vs-first)
   contrast_formulas <- combn(conditions, 2, FUN = function(x) paste(x[1], "-", x[2]), simplify = FALSE)
    contrast_names <- combn(conditions, 2, FUN = function(x) paste(x[1], "vs", x[2], sep="_"), simplify = FALSE)
    contrasts <- setNames(contrast_formulas, contrast_names)
    contrast_matrix <- makeContrasts(levels = design, contrasts = contrasts)
    rownames(contrast_matrix) <- paste0("(Intercept)exp_condition", rownames(contrast_matrix))

# Test the hypotheses
proteins <- hypothesisTest(object=proteins, i="proteinNorm", contrast=contrast_matrix)
peptides <- hypothesisTest(object=peptides, i="peptideNorm", contrast=contrast_matrix)

## adding the new columns to the data frame
stand_prot_out <- assay(proteins[["proteinNorm"]])
add_cols <- rowData(proteins[["proteinNorm"]])
stand_prot_out <- cbind(add_cols[, grep("^protein_group", colnames(add_cols))], stand_prot_out, add_cols[, grep("^number_of_peptides_", colnames(add_cols))])
for (i in unlist(contrasts)) {
    ttt <- add_cols[, i]
    colnames(ttt) <- paste0(c("log_ratios_", "standard_error_", "degrees_freedom_",
                              "t_", "differential_abundance_pvalue_", "differential_abundance_qvalue_"), names(i))
    stand_prot_out <- cbind(stand_prot_out, ttt[,c(1,5,6)])
}
stand_pep_out <- assay(peptides[["peptideNorm"]])
add_cols <- rowData(peptides[["peptideNorm"]])
stand_pep_out <- cbind(add_cols[, grep("^modified_peptide", colnames(add_cols))], stand_pep_out, add_cols[, grep("^number_of_psms_", colnames(add_cols))])
for (i in unlist(contrasts)) {
    ttt <- add_cols[, i]
    colnames(ttt) <- paste0(c("log_ratios_", "standard_error_", "degrees_freedom_",
                              "t_", "differential_abundance_pvalue_", "differential_abundance_qvalue_"), names(i))
    stand_pep_out <- cbind(stand_pep_out, ttt[,c(1,5,6)])
}
write.csv(stand_prot_out, "stand_prot_quant_merged.csv", row.names = F)
write.csv(stand_pep_out, "stand_pep_quant_merged.csv", row.names = F)
