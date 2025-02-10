#!/usr/bin/env Rscript

#### simple script to read Proline result files and set the parameter values for PolySTest
library(readxl)

# checking flags an input
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("\n**** Script to read Proline output, take the log, to configure the PolySTest parameter file and to write csv-files to be read by PolySTest (both protein and peptide level).****
You need to specify these two files in the given order: experimental design (txt-file) and Excel output file from Proline.
The output are two parameter files pep_param.yml and prot_param_yml that can be read by runPolySTestCLI.R as well as the files peptide_ions_profline.csv
peptides_proline.csv and proteins_proline.csv
       ", call. = FALSE)
}

# reading experimental design
filename <- args[1]
if (!file.exists(filename)) {
  stop("experimental design file does not exist", call. = FALSE)
}
expDesign <- read.csv(filename, sep = "\t")

# Proline output file
filename <- args[2]
if (!file.exists(filename)) {
  stop("Proline output file does not exist", call. = FALSE)
}

# reading peptide ions
peptide_ions <- as.data.frame(read_excel(filename, "Quantified peptide ions"))

# reading protein table
proteins <- as.data.frame(read_excel(filename, "Protein sets"))

## Display experimental design
NumReps <- max(table(expDesign[, 2]))
NumCond <- length(unique(expDesign[, 2]))

cat(paste0("****** Found ", NumCond, " different experimental conditions with up to ", NumReps, " replicates\n"))
knitr::kable(table(expDesign[, 2]))

# create names of samples
sample_names <- basename(expDesign[, 1])
sample_names <- gsub(".mzDB", "", sample_names)
sample_names <- gsub(".raw", "", sample_names)
sample_names <- paste0("abundance_", sample_names)
expDesign[, 1] <- sample_names
expDesign <- expDesign[, 1:2]
colnames(expDesign) <- c("sample", "condition")

# adding replicate numbers
tmpDesign <- NULL
NA_cols <- 0
for (cond in unique(expDesign[, 2])) {
  allreps <- expDesign[expDesign[, 2] == cond, , drop = F]
  tmpDesign <- rbind(tmpDesign, cbind(allreps, replicate = 1:nrow(allreps)))
  if (NumReps > nrow(allreps)) {
    numcols <- NumReps - nrow(allreps)
    tmpDesign <- rbind(tmpDesign, cbind(sample = paste0("NA column ", (NA_cols + 1):(NA_cols + numcols)), condition = cond, replicate = (nrow(allreps) + 1):NumReps))
    NA_cols <- NA_cols + numcols
  }
}
expDesign <- tmpDesign
rownames(expDesign) <- paste(expDesign$condition, expDesign$replicate)

### On peptide level
# set number of columns before input
ColQuant <- min(which(colnames(peptide_ions) %in% expDesign$sample))
print(expDesign)
print(colnames(peptide_ions))
peptides <- peptide_ions[, 1:(ColQuant - 1)]
for (rep in 1:NumReps) {
  for (cond in unique(expDesign$condition)) {
    sample <- expDesign[paste(cond, rep), 1]
    if (any(sample == colnames(peptide_ions))) {
      dat <- log2(peptide_ions[, sample])
      peptides <- cbind(peptides, dat)
    } else {
      peptides <- cbind(peptides, NA)
    }
    colnames(peptides)[ncol(peptides)] <- paste("abundance", cond, rep, sep = "_")
  }
}
write.csv(peptide_ions, "polystest_ions_res.csv", row.names = F)

write.csv(peptides, "peptides_proline.csv", row.names = F)

## Write parameter values for PolySTest, assuming unpaired design
params <- paste0(
  "numreps: ", NumReps, "\nnumcond: ", NumCond, "\npaired: false\nrefcond: 0\nfirstquantcol: ", ColQuant,
  "\nrep_grouped: true\ncsvfile: peptides.csv\ndelim: ','\ndecimal: '.'\nheader: true\noutfile: 'polystest_pep_res.csv'\nthreads: 2\n"
)
writeLines(params, "pep_param.yml")
cat(paste0("Files pep_param.yml , peptide_ions_proline.csv and peptides_proline.txt written\n"))

## now on protein level
ColQuant <- min(which(colnames(proteins) %in% expDesign$sample))
oldproteins <- proteins
proteins <- proteins[, 1:(ColQuant - 1)]
for (rep in 1:NumReps) {
  for (cond in unique(expDesign$condition)) {
    sample <- expDesign[paste(cond, rep), 1]
    if (any(sample == colnames(oldproteins))) {
      dat <- log2(oldproteins[, sample])
      # print(head(dat))
      # dat <- dat - median(dat, na.rm=T)
      proteins <- cbind(proteins, dat)
    } else {
      proteins <- cbind(proteins, NA)
    }
    colnames(proteins)[ncol(proteins)] <- paste(cond, rep)
  }
}
write.csv(proteins, "proteins.csv", row.names = FALSE)

## Write parameter values for PolySTest, assuming unpaired design
params <- paste0(
  "numreps: ", NumReps, "\nnumcond: ", NumCond, "\npaired: false\nrefcond: 0\nfirstquantcol: ", ColQuant,
  "\nrep_grouped: true\ncsvfile: proteins.csv\ndelim: ','\ndecimal: '.'\nheader: true\noutfile: 'polystest_prot_res.csv'\nthreads: 2\n"
)
writeLines(params, "prot_param.yml")
cat(paste0("Files prot_param.yml and proteins_proline.txt written\n"))
