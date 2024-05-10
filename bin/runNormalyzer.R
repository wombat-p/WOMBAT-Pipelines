## Script for running Normalyzer on peptide and protein data from WOMBAT-P

# Getting general methods for statistical analysis
source("stat_funcs.R")

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

# Read experimental design
exp_design <- read_expdesign(exp_file)

# Read proteins
proteins <- read_proteins("/data/tmp/stand_prot_quant.csv")

# Read peptides
peptides <- read_peptides("/data/tmp/stand_pep_quant.csv")

# Check if all names are present
all_names <- check_names(exp_design, proteins)
all_names <- check_names(exp_design, peptides)

# Write experimental design for Normalyzer
# take advantage of the fact that all_names is in the same order as the rows in the protein and peptide tables
groups <- match(all_names, colnames(proteins))
final_exp <- data.frame(sample = all_names, group = exp_design$exp_condition, quote = F)
write.table(final_exp, "Normalyzer_design.tsv", sep = "\t", row.names = F, quote = F)

# Write protein and peptide files for Normalyzer
write.table(cbind(protein_group = rownames(proteins), proteins),
  "protein_file.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(cbind(modified_peptide = rownames(peptides), peptides),
  "peptide_file.txt",
  sep = "\t", quote = F, row.names = F
)

# comparison set to everything versus first
if (comps == "") {
  # compfile <- "Normalyzer_comparisons.txt"
  if (file.exists(compfile)) {
    print("Reading comparisons from file")
    comps <- readChar(compfile, file.info(compfile)$size)
    comps <- gsub("[\r\n]", "", comps)
    comps <- unlist(strsplit(comps, ","))
    names(comps) <- sub("-", "_vs_", comps)
  } else {
    print("Comparing against first condition")
    comps <- unique(final_exp[, "group"])
    cmp <- comps[seq_len(length(comps) - 1) + 1]
    ref <- comps[1]
    comps <- paste0(cmp, "-", ref)
    names(comps) <- paste0(cmp, "_vs_", ref)
  }
} else {
  comps <- unlist(strsplit(comps, ","))
}


## run Normalyzer
if (min(table(final_exp$group)) > 1 & length(unique(final_exp$group)) > 1) {
  NormalyzerDE::normalyzer(
    jobName = "NormalyzerProteins",
    designPath = "Normalyzer_design.tsv",
    dataPath = "protein_file.txt", zeroToNA = TRUE,
    outputDir = "./", requireReplicates = FALSE
  )
  NormalyzerDE::normalyzer(
    jobName = "NormalyzerPeptides",
    designPath = "Normalyzer_design.tsv",
    dataPath = "peptide_file.txt", zeroToNA = TRUE,
    outputDir = "./", requireReplicates = FALSE
  )

  print("Now running differential expression analysis")
  print(paste0("./NormalyzerProteins/", normalyzerMethod, "-normalized.txt"))
  NormalyzerDE::normalyzerDE(
    jobName = "NormalyzerProteins",
    comparisons = comps, designPath = "Normalyzer_design.tsv",
    dataPath = paste0(
      "./NormalyzerProteins/",
      normalyzerMethod, "-normalized.txt"
    ),
    outputDir = "./", leastRepCount = "0"
  )
  NormalyzerDE::normalyzerDE(
    jobName = "NormalyzerPeptides",
    comparisons = comps, designPath = "Normalyzer_design.tsv",
    dataPath = paste0(
      "./NormalyzerPeptides/",
      normalyzerMethod, "-normalized.txt"
    ),
    outputDir = "./", leastRepCount = "0"
  )
} else {
  NormalyzerDE::normalyzer(
    jobName = "NormalyzerProteins",
    designPath = "Normalyzer_design.tsv",
    dataPath = "protein_file.txt", zeroToNA = TRUE,
    inputFormat = "maxquantprot",
    outputDir = "./", requireReplicates = F, skipAnalysis = T
  )
  NormalyzerDE::normalyzer(
    jobName = "NormalyzerPeptides",
    designPath = "Normalyzer_design.tsv",
    dataPath = "peptide_file.txt", zeroToNA = TRUE,
    inputFormat = "maxquantpep",
    outputDir = "./", requireReplicates = F, skipAnalysis = T
  )
  print("No statistical testing as at least one sample group
    with only 1 replicate or only one sample group")
  write.csv(NA, "NormalyzerProteins/Normalyzer_stats.tsv")
  write.csv(NA, "NormalyzerPeptides/Normalyzer_stats.tsv")
}

## Preparing for standardized format
# Reading files
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
peptides$protein_group <- peptides$Proteins
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

# adding _vs_ to comparisons
proteins$protein_group <- rownames(proteins)
comps[1:length(comps)] <- make.names(comps)
for (s in length(comps)) {
  col_pos <- grep(comps[s], colnames(proteins))
  colnames(proteins)[col_pos] <- sub(comps[s], names(comps)[s], colnames(proteins)[col_pos])
  col_pos <- grep(comps[s], colnames(peptides))
  colnames(peptides)[col_pos] <- sub(comps[s], names(comps)[s], colnames(peptides)[col_pos])
}



write.csv(proteins, "stand_prot_quant_merged.csv", row.names = F)
write.csv(peptides, "stand_pep_quant_merged.csv", row.names = F)
exp_design_out <- final_exp[, c("sample", "group")]
colnames(exp_design_out) <- c("raw_file", "exp_conditions")
write.table(exp_design_out, "exp_design_calcb.tsv", quote = F, sep = "\t", row.names = F)

cat("Done\n")
