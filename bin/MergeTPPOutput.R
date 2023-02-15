# read in the experimental design file
xp_design <- read.csv("exp_design.txt", sep = "\t")
# sort the experimental design file by the first column
exp_design <- xp_design[order(xp_design[, 1]), , drop = FALSE]
# set the row names to the first column
rownames(exp_design) <- exp_design[, 1]
# add a third column with indices for replicate numbers
exp_design <- data.frame(exp_design, replicates = 1)
counts <- vector("numeric", length(unique(exp_design[, 2])))
names(counts) <- unique(exp_design[, 2])
for (exp in 1:nrow(exp_design)) {
  counts[exp_design[exp, 2]] <- counts[exp_design[exp, 2]] + 1
  exp_design$replicates[exp] <- counts[exp_design[exp, 2]]
}


## Protein tablee
# read and merge all output files from StPeter output converted to csv
quant_out <- quant_pep_out <- list()
keep_columns_once <- c("protein_name", "n_indistinguishable_proteins")


keep_columns_all <- c(
  "probability", "total_number_peptides",
  "total_number_distinct_peptides",
  "Spectral.index",
  "Normalized.spectral.index", "number.of.quantified.peptides",
  "Normalized.spectral.abundance.factor"
)
keep_pep_columns_once <- c(
  "modified_peptide", "peptide_sequence", "protein_name", "charge",
  "nsp_adjusted_probability", "calc_neutral_pep_mass",
  "StPeterQuant_peptide.charge"
)
keep_pep_columns_all <- c(
  "fileName", "pct_spectrum_ids", "confidence", "initial_probability",
  "nsp_adjusted_probability", "fpkm_adjusted_probability", "weight",
  "is_nondegenerate_evidence", "n_instances", "exp_tot_instances",
  "is_contributing_evidence",
  "StPeterQuant_peptide.sequence", "StPeterQuant_peptide.SI",
  "StPeterQuant_peptide.SC"
)

prot_convert_names <- c(
  protein_group = "protein_name",
  number_of_peptides = "total_number_distinct_peptides",
  abundance = "Normalized.spectral.index"
)

pep_convert_names <- c(
  modified_peptide = "modified_peptide",
  protein_group = "protein_name", charge = "StPeterQuant_peptide.charge",
  abundance = "StPeterQuant_peptide.SI"
)

# combine all quantified proteins in one table
# and remove duplicates
quant <- do.call(cbind, quant_out)
quant <- quant[!duplicated(quant[, "protein_name"]), ]
# combine all quantified peptides in one table
# and remove duplicates
quant_pep <- do.call(cbind, quant_pep_out)
quant_pep <- quant_pep[!duplicated(quant_pep[, "modified_peptide"]), ]

# to map protein groups to their info
prot_info <- pep_info <- NULL
for (file in exp_design[, 1]) {
  t_quant <- read.csv(sub(
    ".raw",
    ".pep.interact.pep.prot_stpeter.prot_prot.csv", file
  ))
  t_pep_quant <- read.csv(sub(
    ".raw",
    ".pep.interact.pep.prot_stpeter.prot_pep.csv", file
  ))
  # filter out non-quantified peptides
  t_pep_quant <- t_pep_quant[
    !is.na(t_pep_quant[, "StPeterQuant_peptide.sequence"]),
  ]
  # remove sporadic multiple entries
  # (not sure why they appear, seem to be some modified peptides)
  t_pep_quant <- t_pep_quant[
    !duplicated(apply(t_pep_quant[, c("modified_peptide", "charge")],
      1, paste,
      collapse = ";"
    )),
  ]

  t_pep_rownames <- apply(t_pep_quant[, c("modified_peptide", "charge")],
    1, paste,
    collapse = ";"
  )

  rownames(t_pep_quant) <- t_pep_rownames
  rownames(t_quant) <- t_quant[, "protein_name"]
  t_prot_info <- t_quant[, keep_columns_once, drop = F]
  t_pep_info <- t_pep_quant[, keep_pep_columns_once]
  quant_out[[file]] <- cbind(rownames(t_quant),
    t_quant[, keep_columns_all],
    stringsAsFactors = FALSE
  )
  print(colnames(t_pep_quant))
  quant_pep_out[[file]] <- cbind(rownames(t_pep_quant),
    t_pep_quant[, keep_pep_columns_all],
    stringsAsFactors = FALSE
  )
  prot_info <- rbind(prot_info, t_prot_info)
  pep_info <- rbind(pep_info, t_pep_info)
  colnames(quant_out[[file]]) <- paste(colnames(quant_out[[file]]),
    exp_design[file, 2], exp_design[file, 3],
    sep = "_"
  )
  colnames(quant_pep_out[[file]]) <- paste(colnames(quant_pep_out[[file]]),
    exp_design[file, 2], exp_design[file, 3],
    sep = "_"
  )
}

prot_info <- unique(prot_info)
pep_info <- unique(pep_info)
all_quant <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), quant_out)
all_pep_quant <- Reduce(function(x, y) {
  merge(x, y, by = 1, all = TRUE)
}, quant_pep_out)
all_quant <- cbind(
  prot_info[all_quant[, 1], ],
  all_quant[, 2:ncol(all_quant)]
)
all_pep_quant <- cbind(
  pep_info[all_pep_quant[, 1], ],
  all_pep_quant[, 2:ncol(all_pep_quant)]
)

# reorder columns
ind_cols <- colnames(all_quant)[
  -which(colnames(all_quant) %in% keep_columns_once)
]
all_quant <- all_quant[
  , c(keep_columns_once, ind_cols[order(ind_cols)])
]
ind_cols <- colnames(all_pep_quant)[
  -which(colnames(all_pep_quant) %in% keep_pep_columns_once)
]
all_pep_quant <- all_pep_quant[
  , c(keep_pep_columns_once, ind_cols[order(ind_cols)])
]

write.csv(all_quant, "all_prot_quant_merged.csv", row.names = FALSE)
write.csv(all_pep_quant, "all_pep_quant_merged.csv", row.names = FALSE)


# Now simplifying the tables into new standard format for benchmarking
stand_prot_quant <- all_quant[, grep(
  paste("^", prot_convert_names, collapse = "|", sep = ""),
  colnames(all_quant)
)]
cnames <- colnames(stand_prot_quant)
for (i in 1:length(prot_convert_names)) {
  cnames <- gsub(prot_convert_names[i], names(prot_convert_names)[i], cnames)
}
colnames(stand_prot_quant) <- cnames

stand_pep_quant <- all_pep_quant[, grep(
  paste("^", pep_convert_names, collapse = "|", sep = ""),
  colnames(all_pep_quant)
)]
cnames <- colnames(stand_pep_quant)
for (i in 1:length(pep_convert_names)) {
  cnames <- gsub(pep_convert_names[i], names(pep_convert_names)[i], cnames)
}
colnames(stand_pep_quant) <- cnames

# merging charge states
library(dplyr)
stand_pep_quant <- merge(
  stand_pep_quant %>%
    group_by(modified_peptide) %>%
    summarize_at(c("protein_group", "charge"), paste, collapse = ";"),
  stand_pep_quant %>%
    group_by(modified_peptide) %>%
    summarize_at(grep("^abundance", colnames(stand_pep_quant), value = TRUE),
      sum,
      na.rm = TRUE
    ),
  by = 1
)

for (c in grep("^abundance", colnames(stand_pep_quant))) {
  stand_pep_quant[stand_pep_quant[, c] == 0, c] <- NA
}

write.csv(stand_prot_quant, "stand_prot_quant_merged_pre.csv",
  row.names = FALSE
)
write.csv(stand_pep_quant, "stand_pep_quant_merged_pre.csv",
  row.names = FALSE
)
