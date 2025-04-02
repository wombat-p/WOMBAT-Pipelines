library(tidyr)
library(dplyr)

# reading cmd arguments
args <- commandArgs(trailingOnly = TRUE)
normalyzerMethod <- strsplit(grep("--method", args, value = TRUE), split = "=")[[1]][[2]]
comps <- strsplit(grep("--comps", args, value = TRUE), split = "=")[[1]]
if (length(comps) > 1) {
  comps <- comps[[2]]
} else {
  comps <- ""
}
exp_file <- strsplit(grep("--exp_design", args, value = TRUE), split = "=")[[1]][[2]]
compfile <- strsplit(grep("--comp_file", args, value = TRUE), split = "=")[[1]][[2]]

# merging experimental design file from sdrf parser with actual design
final_exp <- read.csv(exp_file, sep = "\t")

# Reduce to unique rows
final_exp <- final_exp[!duplicated(final_exp[, !grepl("Assay|Run", colnames(final_exp))]), , drop = F]

# Create column for (biological) replicate number if not existing already
if (is.null(final_exp$biorep)) {
  final_exp$biorep <- 1
  for (i in unique(final_exp$group)) {
    ttt <- final_exp[final_exp$group == i, "biorep"]
    final_exp[final_exp$group == i, "biorep"] <- 1:length(ttt)
  }
}

# Write to file again
write.table(final_exp, exp_file, sep = "\t", quote = F, row.names = F)

# comparison set to everything versus first
if (comps == "") {
  # compfile <- "Normalyzer_comparisons.txt"
  if (file.exists(compfile)) {
    print("Reading comparisons from file")
    comps <- readChar(compfile, file.info(compfile)$size)
    comps <- gsub("[\r\n]", "", comps)
    comps <- unlist(strsplit(comps, ","))
  } else {
    print("Comparing against first condition")
    comps <- unique(final_exp[, "group"])
    comps <- paste0(comps[2:length(comps)], "-", comps[1])
  }
} else {
  comps <- unlist(strsplit(comps, ","))
}

# Reduce protein accessions from long format (e.g. "sp|P12345|A1BG_HUMAN;sp|P12346|A1BG_HUMAN") to a string of only the accession numbers
reduce_prot_accs <- function(accessions) {
  tout <- sapply(accessions, function(y) {
    tgroup <- unlist(strsplit(y, ";"))
    tgroup <- lapply(tgroup, function(x) {
      if (is.na(x)) {
        return(NA)
      }
      parts <- strsplit(x, "\\|")[[1]]

      # If it has exactly 3 parts, the second is the accession
      if (length(parts) == 3) {
        return(parts[2])
      } else {
        # Otherwise, return the original entry
        return(x)
      }
    })
    return(paste(tgroup, collapse = ","))
  })

  return(tout)
}


## Read evidence file to create ion standard file and make proforma files
evidence_data <- read.delim("evidence.txt", stringsAsFactors = FALSE)
# Remove reverse hits
evidence_data <- evidence_data[evidence_data$Reverse != "+", ]
# change "Modified sequence" to proforma format
# for that remove the "site" starts with a space and ends with a bracket
evidence_data$Modified.sequence <- gsub("\\s\\([^()]*\\)", "", evidence_data$Modified.sequence)
# Remove underscores
evidence_data$Modified.sequence <- gsub("_", "", evidence_data$Modified.sequence)
# change parenthesis to brackets
evidence_data$Modified.sequence <- gsub("\\(", "[", evidence_data$Modified.sequence)
evidence_data$Modified.sequence <- gsub("\\)", "]", evidence_data$Modified.sequence)
# When the first character is a bracket, add a "-" after the first "]"
evidence_data$Modified.sequence <- sapply(evidence_data$Modified.sequence, function(x) {
  if (grepl("^\\[", x)) {
    x <- sub("\\]", "]-", x)
  }
  return(x)
})
evidence_data$Proteins <- reduce_prot_accs(evidence_data$Proteins)


# Create ion standard file
std_ion_output <- evidence_data[, c("Modified.sequence", "Proteins", "Experiment", "Charge", "Intensity")]
colnames(std_ion_output) <- c("modified_peptide", "protein_group", "exp_conditions", "charge", "abundance")
# Change to wide format
std_ion_wide <- as.data.frame(pivot_wider(
  std_ion_output,
  names_from = exp_conditions,
  values_from = abundance,
  names_prefix = "abundance_",
  values_fn = ~ sum(.x, na.rm = TRUE)
))
# Filter out the rows with all abundance values being NA or zero
std_ion_wide <- std_ion_wide[rowSums(std_ion_wide[, grep("^abundance_", colnames(std_ion_wide))], na.rm = T) > 0, ]


write.csv(std_ion_wide, "stand_ion_quant_merged.csv", row.names = F)

# Create peptidoform level file from std_ion_output
pep <- as.data.frame(std_ion_wide %>%
  group_by(modified_peptide, protein_group) %>%
  summarise(
    across(
      starts_with("abundance_"),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ))
pep_vals <- pep[, grep("^abundance_", colnames(pep))]
pep_vals[pep_vals == 0] <- NA
pep[, grep("^abundance_", colnames(pep))] <- pep_vals

pep_mq <- pep
colnames(pep_mq)[1:2] <- c("Sequence", "Proteins")
colnames(pep_mq)[2:ncol(pep_mq)] <- sub("^abundance_", "Intensity ", colnames(pep_mq)[2:ncol(pep_mq)])
pep_mq$Mass <- pep_mq$Leading.razor.protein <- pep_mq$Charges <- pep_mq$PEP <- NA
write.table(pep_mq, "peptide_file.txt", row.names = F, sep = "\t", quote = F)

## run Normalyzer
if (min(table(final_exp[, "group"])) > 1 & length(unique(final_exp[, "group"])) > 1) {
  NormalyzerDE::normalyzer(jobName = "NormalyzerProteins", designPath = "Normalyzer_design.tsv", dataPath = "protein_file.txt", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir = "./", requireReplicates = F)
  NormalyzerDE::normalyzer(jobName = "NormalyzerPeptides", designPath = "Normalyzer_design.tsv", dataPath = "peptide_file.txt", zeroToNA = TRUE, inputFormat = "maxquantpep", outputDir = "./", requireReplicates = F)
  print("Now running differential expression analysis")
  print(paste0("./NormalyzerProteins/", normalyzerMethod, "-normalized.txt"))
  NormalyzerDE::normalyzerDE(jobName = "NormalyzerProteins", comparisons = comps, designPath = "Normalyzer_design.tsv", dataPath = paste0("./NormalyzerProteins/", normalyzerMethod, "-normalized.txt"), outputDir = "./", leastRepCount = "0")
  NormalyzerDE::normalyzerDE(jobName = "NormalyzerPeptides", comparisons = comps, designPath = "Normalyzer_design.tsv", dataPath = paste0("./NormalyzerPeptides/", normalyzerMethod, "-normalized.txt"), outputDir = "./", leastRepCount = "0")
} else {
  NormalyzerDE::normalyzer(jobName = "NormalyzerProteins", designPath = "Normalyzer_design.tsv", dataPath = "protein_file.txt", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir = "./", requireReplicates = F, skipAnalysis = T)
  NormalyzerDE::normalyzer(jobName = "NormalyzerPeptides", designPath = "Normalyzer_design.tsv", dataPath = "peptide_file.txt", zeroToNA = TRUE, inputFormat = "maxquantpep", outputDir = "./", requireReplicates = F, skipAnalysis = T)
  print("No statistical testing as at least one sample group with only 1 replicate or only one sample group")
  write.csv(NA, "NormalyzerProteins/Normalyzer_stats.tsv")
  write.csv(NA, "NormalyzerPeptides/Normalyzer_stats.tsv")
}

## Preparing for standardized format
# Reading files
peptides <- read.csv("peptide_file.txt", sep = "\t", row.names = 1)
proteins <- read.csv("protein_file.txt", sep = "\t", row.names = 1)
norm_peptides <- read.csv(paste0("NormalyzerPeptides/", normalyzerMethod, "-normalized.txt"), sep = "\t", row.names = 1)
norm_proteins <- read.csv(paste0("NormalyzerProteins/", normalyzerMethod, "-normalized.txt"), sep = "\t", row.names = 1)
stats_peptides <- stats_proteins <- NULL
if (file.exists("NormalyzerPeptides/NormalyzerPeptides_stats.tsv")) {
  stats_peptides <- read.csv("NormalyzerPeptides/NormalyzerPeptides_stats.tsv", sep = "\t", row.names = 1)
  stats_proteins <- read.csv("NormalyzerProteins/NormalyzerProteins_stats.tsv", sep = "\t", row.names = 1)
} else {
  stats_peptides <- read.csv(paste0("NormalyzerPeptides/", normalyzerMethod, "-normalized.txt"), sep = "\t", row.names = 1)
  stats_proteins <- read.csv(paste0("NormalyzerProteins/", normalyzerMethod, "-normalized.txt"), sep = "\t", row.names = 1)
}


# changing column names
peptides$missed_cleavages <- peptides$Missed.cleavages
peptides$charge <- peptides$Charges
peptides$protein_group <- reduce_prot_accs(peptides$Proteins)
if (any(grepl("PValue$", colnames(stats_peptides)))) {
  pval_cols <- colnames(stats_peptides)[grep("AdjPVal$", colnames(stats_peptides))]
  colnames(stats_peptides)[grep("AdjPVal$", colnames(stats_peptides))] <-
    paste0("differential_abundance_qvalue_", sub("_AdjPVal$", "", pval_cols))
  pval_cols <- colnames(stats_proteins)[grep("AdjPVal$", colnames(stats_proteins))]
  colnames(stats_proteins)[grep("AdjPVal$", colnames(stats_proteins))] <-
    paste0("differential_abundance_qvalue_", sub("_AdjPVal$", "", pval_cols))
  pval_cols <- colnames(stats_peptides)[grep("PValue$", colnames(stats_peptides))]
  colnames(stats_peptides)[grep("PValue$", colnames(stats_peptides))] <-
    paste0("differential_abundance_pvalue_", sub("_PValue$", "", pval_cols))
  pval_cols <- colnames(stats_proteins)[grep("PValue$", colnames(stats_proteins))]
  colnames(stats_proteins)[grep("PValue$", colnames(stats_proteins))] <-
    paste0("differential_abundance_pvalue_", sub("_PValue$", "", pval_cols))
  pval_cols <- colnames(stats_peptides)[grep("_log2FoldChange$", colnames(stats_peptides))]
  colnames(stats_peptides)[grep("_log2FoldChange$", colnames(stats_peptides))] <-
    paste0("log_fold_change_", sub("_log2FoldChange$", "", pval_cols))
  pval_cols <- colnames(stats_proteins)[grep("_log2FoldChange$", colnames(stats_proteins))]
  colnames(stats_proteins)[grep("_log2FoldChange$", colnames(stats_proteins))] <-
    paste0("log_fold_change_", sub("_log2FoldChange$", "", pval_cols))
}

colnames(peptides) <- unlist(sub("^Experiment\\.", "number_of_psms_", colnames(peptides)))

# colnames(peptides) <- unlist(sub("^LFQ\\.intensity\\.","abundance_", colnames(peptides)))
# for (i in 1:nrow(final_exp)) {
#   colnames(peptides) <- unlist(sub(final_exp[i,1], final_exp[i,2], colnames(peptides)))
#   colnames(proteins) <- unlist(sub(final_exp[i,1], final_exp[i,2], colnames(proteins)))
# }

# changing generic names to experimental design
colnames(proteins) <- unlist(sub("^Razor\\.\\.\\.unique\\.peptides\\.", "number_of_peptides_", colnames(proteins)))
proteins$protein_group <- rownames(proteins)
for (s in 1:nrow(final_exp)) {
  substitute_from <- paste0(make.names(final_exp$source_name[s]), "_Tr_", final_exp$technical_replicate[s])
  substitute_with <- paste0(make.names(final_exp$group[s]), "_", final_exp$biorep[s])
  # avoid setting X in front
  if (grepl("^X", substitute_with)) substitute_with <- sub("^X", "", substitute_with)
  colnames(proteins) <- sub(substitute_from, substitute_with, colnames(proteins))
  colnames(peptides) <- sub(substitute_from, substitute_with, colnames(peptides))
  colnames(norm_proteins) <- sub(paste0("^", substitute_from), paste0("abundance_", substitute_with), colnames(norm_proteins))
  colnames(norm_peptides) <- sub(paste0("^", substitute_from), paste0("abundance_", substitute_with), colnames(norm_peptides))
}
colnames(proteins) <- make.unique(colnames(proteins))
colnames(peptides) <- make.unique(colnames(peptides))
colnames(norm_proteins) <- make.unique(colnames(norm_proteins))
colnames(norm_peptides) <- make.unique(colnames(norm_peptides))


# getting relevant columns
norm_peptides[, grep("^abundance_", colnames(norm_peptides), value = T)] <- 2^(norm_peptides[, grep("^abundance_", colnames(norm_peptides), value = T)])
if (!any(grepl("^differential_abundance", colnames(stats_peptides)))) {
  proteins <- cbind(
    proteins[rownames(norm_proteins), c(
      "protein_group",
      grep("^number_of_peptides_", colnames(proteins), value = T)
    )],
    norm_proteins[, grep("^abundance_", colnames(norm_proteins), value = T)]
  )
  peptides <- cbind(
    modified_peptide = rownames(norm_peptides),
    peptides[rownames(norm_peptides), c(
      "protein_group",
      grep("^number_of_psms_", colnames(peptides), value = T)
    ), drop = F],
    norm_peptides[, grep("^abundance_", colnames(norm_peptides), value = T)]
  )
} else {
  proteins <- cbind(
    proteins[rownames(norm_proteins), c(
      "protein_group",
      grep("^number_of_peptides_", colnames(proteins), value = T)
    )],
    stats_proteins[rownames(norm_proteins), grep("^differential_abundance", colnames(stats_proteins), value = T)],
    norm_proteins[, grep("^abundance_", colnames(norm_proteins), value = T)]
  )
  peptides <- cbind(
    modified_peptide = rownames(norm_peptides),
    peptides[rownames(norm_peptides), c(
      "protein_group",
      grep("^number_of_psms_", colnames(peptides), value = T)
    ), drop = F],
    norm_peptides[, grep("^abundance_", colnames(norm_peptides), value = T)],
    stats_peptides[
      rownames(norm_peptides),
      grep("^differential_abundance", colnames(stats_peptides), value = T)
    ]
  )
}

proteins$protein_group <- reduce_prot_accs(proteins$protein_group)

write.csv(proteins, "stand_prot_quant_merged.csv", row.names = F)
write.csv(peptides, "stand_pep_quant_merged.csv", row.names = F)
exp_design_out <- final_exp[, c("Run", "group")]
colnames(exp_design_out) <- c("raw_file", "exp_conditions")
write.table(exp_design_out, "exp_design_calcb.tsv", quote = F, sep = "\t", row.names = F)

cat("Done\n")
