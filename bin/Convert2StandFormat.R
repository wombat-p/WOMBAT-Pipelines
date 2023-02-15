library(matrixStats)
library(stringi)

#Reading files
peptides <- read.csv("polystest_pep_res.csv")
proteins <- read.csv("polystest_prot_res.csv")
exp_design <- read.csv("exp_design.txt", sep="\t")
exp_design[,2] <- make.names(exp_design[,2])

if(is.null(exp_design$raw_file))
     exp_design$raw_file <- exp_design$mzdb_file

exp_design$sample_name <- sub("\\.mzDB","", sub("\\.\\/", "", exp_design$raw_file))
colnames(exp_design)[1:2] <- c("raw_file", "exp_condition")
# Create column for (biological) replicate number if not existing already
if (is.null(exp_design$biorep)) {
exp_design$biorep <- 1
for (i in unique(exp_design$exp_condition)) {
  ttt <- exp_design[exp_design$exp_condition == i, "biorep"]
  exp_design[exp_design$exp_condition == i, "biorep"] <- 1:length(ttt)
}
}

write.table(exp_design, "exp_design.txt", sep="\t", row.names=F)

# Converting column names
for (i in 1:nrow(exp_design)) {
colnames(peptides) <- sub(paste0("^psm_count_", exp_design$sample_name[i],"$"), 
                          paste("number_of_psms", exp_design$exp_condition[i], exp_design$biorep[i], sep="_"), 
                          colnames(peptides))
colnames(proteins) <- sub(paste0("^peptides_count_", exp_design$sample_name[i],"$"), 
                          paste("number_of_peptides", exp_design$exp_condition[i], exp_design$biorep[i], sep="_"), 
                          colnames(proteins))
}
colnames(peptides) <- sub("^log\\.ratios\\.", "log_fold_change_", colnames(peptides))
for (s in unique(exp_design$exp_condition)) colnames(peptides) <- sub(paste0("^",s,"\\."), paste0("abundance_",s,"_"), colnames(peptides))
colnames(peptides) <- sub("^FDR\\.PolySTest\\.", "differential_abundance_qvalue_", colnames(peptides))

#Creating modified sequences
modified_peptides <- strsplit(as.character(peptides$modifications), "; ")
modified_peptides <- lapply(modified_peptides, function(x) {
if(any(!is.na(x))) {
  tt <- matrix(unlist(strsplit(x, " \\(")), nrow=2)
  tt[2,] <- sub("\\)","",tt[2,])
  modpos <- NULL
  for (i in 1:ncol(tt)) {
    modpos <- append(modpos, ifelse(grepl("N-term|C-term", tt[2,i]), 0, sub("[A-Z]","",tt[2,i])))
  }
  tt <- rbind(tt, modpos)
} else {
  NA
}
})

modified_sequence <- peptides$sequence
for (i in 1:nrow(peptides)) {
if (!is.na(modified_peptides[[i]][1])) {
  modified_sequence[i] <- stri_sub_replace_all(modified_sequence[i], replacement = paste0("[",modified_peptides[[i]][1,],"]"), 
                                               from=as.numeric(modified_peptides[[i]][3,])+1, 
                       to=as.numeric(modified_peptides[[i]][3,]))
}
}
peptides$modified_sequence <- modified_sequence

stand_peps <- data.frame("modified_peptide"=peptides$modified_sequence, protein_group=peptides$samesets_accessions, 
                       peptides[, grep("^number_of_psms", colnames(peptides)), drop=F],
                       2^peptides[, grep("^abundance", colnames(peptides)), drop=F],
                       peptides[, grep("^log_ratios", colnames(peptides)), drop=F],
                       peptides[, grep("^differential_abundance_qvalue", colnames(peptides)), drop=F])

# deleting charge states with lower intensities to maintain max. 1 (modified) peptide sequence
stand_peps <- stand_peps[order(rowMeans(peptides[, grep("^abundance", colnames(peptides))], na.rm = T), decreasing=T),]
stand_peps <- stand_peps[!duplicated(stand_peps$modified_peptide), ]
stand_peps <- stand_peps[order(stand_peps$protein_group), ]
write.csv(stand_peps, "stand_pep_quant_merged.csv")

# Converting column names
colnames(proteins) <- sub("^log\\.ratios\\.", "log_ratios_", colnames(proteins))
for (s in unique(exp_design$exp_condition)) colnames(proteins) <- sub(paste0("^",s,"\\."), paste0("abundance_",s,"_"), colnames(proteins))
colnames(proteins) <- sub("^FDR\\.PolySTest\\.[X]?", "differential_abundance_qvalue_", colnames(proteins))
stand_prots <- data.frame(protein_group=proteins$samesets_accessions, 
                        proteins[, grep("^number_of_peptides", colnames(proteins)), drop=F],
                        proteins[, grep("^abundance_", colnames(proteins)), drop=F],
                        proteins[, grep("^log_ratios", colnames(proteins)), drop=F],
                        proteins[, grep("^differential_abundance_qvalue_", colnames(proteins)), drop=F]
                        )
write.csv(stand_prots, "stand_prot_quant_merged.csv")

