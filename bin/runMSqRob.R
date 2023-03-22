library(MSnbase)
library(MSqRob)

# reading cmd arguments
args <- commandArgs(trailingOnly = TRUE)
normalization_method <- strsplit(grep('--normalization', args, value = TRUE), split = '=')[[1]][[2]]
min_peptides <- strsplit(grep('--min_peptides', args, value = TRUE), split = '=')[[1]][[2]]
if (!any(normalization_method == c("sum", "median", "mean", "quantiles","none"))) {
  stop("Invalid normalization method, should be one of: sum, median, mean, quantiles, none")
}
if (any(normalization_method == c("mean","median")))
  normalization_method <- paste0("center.", normalization_method)
if (any(normalization_method == c("quantiles")))
  normalization_method <- "quantiles.robust"


options(stringsAsFactors=T)
peptides <- read2MSnSet("q_input.txt",pattern="Intensity_")
#head(exprs(peptides))
exp_annotation <- read.csv("exp_design.txt",sep="\t")
exp_annotation$raw_file <- tools::file_path_sans_ext(exp_annotation$raw_file)
exp_annotation$genotype <- make.names(exp_annotation$exp_condition)
exp_annotation$run <- paste0("Intensity_",exp_annotation$raw_file)
for (c in 1:ncol(exp_annotation)) {
  exp_annotation[,c] <- as.factor(exp_annotation[,c])
}
# Create column for replicate number if not existing already
if (is.null(exp_annotation$biorep)) {
  exp_annotation$biorep <- 1
  for (i in unique(exp_annotation$exp_condition)) {
    ttt <- exp_annotation[exp_annotation$exp_condition == i, "biorep"]
    exp_annotation[exp_annotation$exp_condition == i, "biorep"] <- 1:length(ttt)
  }
}


## Running MSqRob
# Had to change normalization due to error in preprocesscore
# normalization methods: should be one of “sum”, “max”, “center.mean”, “center.median”, “div.mean”, “div.median”, “diff.median”, “quantiles”, “quantiles.robust”, “vsn”
# we stick here with sum, center.mean, center.median and quantiles
peptides <- preprocess_MSnSet(peptides,accession="Protein.Groups",split=",", useful_properties="Sequence", exp_annotation=exp_annotation, normalisation=normalization_method)

# necessary due to change to R 4.x    
fData(peptides)[,"Protein.Groups"] <- as.factor(fData(peptides)[,"Protein.Groups"])
# Set data.frame for experimental design with columns run genotype biorep
proteins <- MSnSet2protdata(peptides, accession="Protein.Groups")  

# filter for proteins with more than min_peptides peptides
keep_prots <- sapply(proteins@data, function(x) length(unique(x$Sequence))) >= min_peptides
# needed to fix problem when accession numbers are actual numbers
proteins[as.character(proteins@accession[which(keep_prots)])]
#proteins <- proteins[which(keep_prots)]
system.time(protLM <- fit.model(proteins, response="quant_value", fixed=c("genotype"),  random=c("run","Sequence"), add.intercept=TRUE))
#create comparisons vs first
contrasts <- NULL
levels <- as.factor(paste0("genotype",make.names(unique(exp_annotation$genotype))))
# levels <- unique(exp_annotation$genotype)
if (length(levels) > 1) {
  for (i in levels[-1]) {
    contrasts <- append(contrasts, paste(i, levels[1], sep="-"))  
  }
} else {
  contrasts <- levels
}
L <- makeContrast(contrasts=contrasts,levels=as.character(levels))
result <- test.protLMcontrast(protLM, L)
result <- prot.p.adjust(result, method="fdr")
write.csv(result, "MSqRobOut.csv")
result <- as.data.frame(result)

## Merging more data from peptide and protein level
all_peptides <- list()
protnames <- NULL
for (type in exp_annotation$raw_file) {
  tin <- read.csv(paste0(type,"_peptides.txt"), sep="\t")
  protnames <- tin[, c("Modified.Sequence","Protein.s.")]
  tin <- tin[,c( "Modified.Sequence","X.Validated.PSMs")]
  colnames(tin) <- paste0(c("modified_peptide", "number_of_psms"), "_", type)
  all_peptides[[type]] <- tin
}
protnames <- unique(protnames)
rownames(protnames) <- protnames[,1]
all_pep <- Reduce(function(x, y) merge(x, y, by=1, all=TRUE), all_peptides)

rownames(all_pep) <- all_pep[,1]
# merging with quant data
stand_pep_quant <- cbind(all_pep[fData(peptides)$Sequence,], protein_group=protnames[fData(peptides)$Sequence, 2], 2^exprs(peptides))
for (r in 1:nrow(exp_annotation)) {
  colnames(stand_pep_quant) <- sub(paste0("^Intensity_", exp_annotation$raw_file[r], "$"), 
                                   paste0("abundance_", exp_annotation$exp_condition[r], "_",
                                          exp_annotation$biorep[r]), colnames(stand_pep_quant))
  colnames(stand_pep_quant) <- sub(paste0("^number_of_psms_", exp_annotation$raw_file[r], "$"), 
                                   paste0("number_of_psms_", exp_annotation$exp_condition[r], "_",
                                          exp_annotation$biorep[r]), colnames(stand_pep_quant))
}
# standardizing colnames
colnames(stand_pep_quant)[1] <- "modified_peptide"

write.csv(stand_pep_quant, "stand_pep_quant_merged.csv", row.names=F)


# Merging data from peptideshaker, flashlfq and msqrob
quant_prots <- read.csv("q_prot.txt", sep="\t")
rownames(quant_prots) <- quant_prots[, "Protein.Groups"]
quant_prots <- quant_prots[, grep("^Intensity",colnames(quant_prots))]
# making zeroes to NA
quant_prots[quant_prots == 0] <- NA
all_proteins <- list()
for (type in exp_annotation$raw_file) {
  tin <- read.csv(paste0(type,"_proteins.txt"), sep="\t")
  tin <- tin[,c("Protein.Group", "X.Unique.Peptides")]
  colnames(tin) <- c("protein_group", paste0("number_of_peptides_", type))
  all_proteins[[type]] <- tin
}
all_prot <- Reduce(function(x, y) merge(x, y, by=1, all=TRUE), all_proteins)
# merging with quant data
rownames(all_prot) <- all_prot[,1]
stand_prot_quant <- cbind(all_prot[rownames(result),], log2(quant_prots[rownames(result), ]), result[,grep("estimate$|qval$|pval$", colnames(result))])
# Exchanging file name based columns names to the ones defined in the experimental design
for (r in 1:nrow(exp_annotation)) {
  colnames(stand_prot_quant) <- sub(paste0("^Intensity_", exp_annotation$raw_file[r], "$"), 
                                    paste0("abundance_", exp_annotation$exp_condition[r], "_",
                                           exp_annotation$biorep[r]), colnames(stand_prot_quant))
  colnames(stand_prot_quant) <- sub(paste0("^number_of_peptides_", exp_annotation$raw_file[r], "$"), 
                                    paste0("number_of_peptides_", exp_annotation$exp_condition[r], "_",
                                           exp_annotation$biorep[r]), colnames(stand_prot_quant))
}

# Change column names for statistics
ttt <- colnames(stand_prot_quant)[grep("estimate$", colnames(stand_prot_quant))] 
ttt <- sub("estimate$", "", ttt)
if (ttt == "")
  ttt <- paste(rev(unique(exp_annotation$exp_condition)), collapse = "_vs_")
colnames(stand_prot_quant)[grep("estimate$", colnames(stand_prot_quant))] <- paste0("log_fold_change_", ttt)
ttt <- colnames(stand_prot_quant)[grep("qval$", colnames(stand_prot_quant))] 
ttt <- sub("qval$", "", ttt)
if (ttt == "")
  ttt <- paste(rev(unique(exp_annotation$exp_condition)), collapse = "_vs_")
colnames(stand_prot_quant)[grep("qval$", colnames(stand_prot_quant))] <- paste0("differential_abundance_qvalue_", ttt)
ttt <- colnames(stand_prot_quant)[grep("pval$", colnames(stand_prot_quant))] 
ttt <- sub("pval$", "", ttt)
if (ttt == "")
  ttt <- paste(rev(unique(exp_annotation$exp_condition)), collapse = "_vs_")
colnames(stand_prot_quant)[grep("pval$", colnames(stand_prot_quant))] <- paste0("differential_abundance_pvalue_", ttt)
write.csv(stand_prot_quant, "stand_prot_quant_merged.csv", row.names=F)



