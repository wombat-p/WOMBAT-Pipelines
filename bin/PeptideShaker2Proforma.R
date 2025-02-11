library(matrixStats)
library(stringi)

# Reading files
ions <- read.csv("peptideshaker_filtered_out.txt", sep = "\t")
peptides <- read.csv("peptideshaker_peptides_out.txt", sep = "\t")
proteins <- read.csv("peptideshaker_proteins_out.txt", sep = "\t")

print(names(peptides))

# Creating modified sequences
modify_sequence <- function(fmods, vmods, sequence) {
  modified_peptides <- lapply(paste0(fmods, ";", vmods), function(x) strsplit(x, ";")[[1]])
  modified_peptides <- lapply(modified_peptides, function(y) {
    mm <- lapply(y, function(x) {
      if (any(!is.na(x))) {
        # Extract string in parentheses
        modpos <- gregexpr("\\((.*?)\\)", x)
        modpos <- regmatches(x, modpos)
        modpos <- gsub("[()]", "", modpos) # Remove parentheses

        modpos <- unlist(strsplit(modpos, ","))[[1]]
        # Remove from x
        x <- gsub("\\(.*?\\)", "", x)
        # Split by " of "
        x <- strsplit(x, " of ")[[1]][1]
        c(x, modpos)
      } else {
        NA
      }
    })
    mods <- NULL
    for (i in 1:length(mm)) {
      if (!is.na(mm[[i]][1])) {
        mods <- rbind(mods, mm[[i]])
      }
    }
    return(mm)
  })


  modified_sequence <- sequence
  for (i in 1:length(modified_peptides)) {
    x <- unlist(modified_peptides[[i]])
    if (!is.null(x)) {
      modified_sequence[i] <- stri_sub_replace_all(modified_sequence[i],
        replacement = paste0("[", x[1], "]"),
        from = as.numeric(x[2]) + 1,
        to = as.numeric(x[2])
      )
    }
  }

  return(modified_sequence)
}

peptides$modified_sequence <- modify_sequence(peptides$Variable.Modifications, peptides$Fixed.Modifications, peptides$Sequence)
head(peptides, 20)
ions$modified_sequence <- modify_sequence(ions$Variable.Modifications, ions$Fixed.Modifications, ions$Sequence)

# Reduce protein accessions from long format (e.g. "sp|P12345|A1BG_HUMAN;sp|P12346|A1BG_HUMAN") to a string of only the accession numbers
reduce_prot_accs <- function(accessions) {
  tout <- sapply(accessions, function(y) {
    tgroup <- unlist(strsplit(y, "; "))
    tgroup <- lapply(tgroup, function(x) {
      if (is.na(x)) {
        return(NA)
      }
      if (stri_count_fixed(x, "|") != 2) {
        return(x)
      }
      return(unlist(strsplit(x, "\\|"))[2])
    })
    return(paste(tgroup, collapse = ","))
  })

  return(tout)
}

proteins$samesets_accessions <- reduce_prot_accs(proteins$samesets_accessions)
peptides$samesets_accessions <- reduce_prot_accs(peptides$samesets_accessions)
ions$samesets_accessions <- reduce_prot_accs(ions$samesets_accessions)

stand_peps <- data.frame(
  "modified_peptide" = peptides$modified_sequence, protein_group = peptides$samesets_accessions,
  peptides[, grep("^number_of_psms", colnames(peptides)), drop = F],
  2^peptides[, grep("^abundance", colnames(peptides)), drop = F],
  peptides[, grep("^log_ratios", colnames(peptides)), drop = F],
  peptides[, grep("^differential_abundance_qvalue", colnames(peptides)), drop = F]
)

stand_ions <- data.frame(
  modified_peptide = ions$modified_sequence,
  protein_group = ions$samesets_accessions,
  charge = ions$master_quant_peptide_ion_charge,
  ions[, grep("^number_of_psms", colnames(ions)), drop = F],
  ions[, grep("^abundance", colnames(ions)), drop = F]
)

# deleting charge states with lower intensities to maintain max. 1 (modified) peptide sequence
stand_peps <- stand_peps[order(rowMeans(peptides[, grep("^abundance", colnames(peptides))], na.rm = T), decreasing = T), ]
stand_peps <- stand_peps[!duplicated(stand_peps$modified_peptide), ]
stand_peps <- stand_peps[order(stand_peps$protein_group), ]
write.csv(stand_peps, "stand_pep_quant_merged.csv", row.names = F)

# Still needs more adjustments of colnames, ...
write.csv(stand_ions, "stand_ions_quant_merged.csv", row.names = F)


# Converting column names
colnames(proteins) <- sub("^log\\.ratios\\.", "log_ratios_", colnames(proteins))
for (s in unique(exp_design$exp_condition)) colnames(proteins) <- sub(paste0("^", s, "\\."), paste0("abundance_", s, "_"), colnames(proteins))
colnames(proteins) <- sub("^FDR\\.PolySTest\\.[X]?", "differential_abundance_qvalue_", colnames(proteins))
stand_prots <- data.frame(
  protein_group = proteins$samesets_accessions,
  proteins[, grep("^number_of_peptides", colnames(proteins)), drop = F],
  proteins[, grep("^abundance_", colnames(proteins)), drop = F],
  proteins[, grep("^log_ratios", colnames(proteins)), drop = F],
  proteins[, grep("^differential_abundance_qvalue_", colnames(proteins)), drop = F]
)
write.csv(stand_prots, "stand_prot_quant_merged.csv", row.names = F)
