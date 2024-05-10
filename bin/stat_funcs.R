## Functions to read and process quantitative protein and peptide data from WOMBAT-P

library(limma)
library(ROTS)

rotsparam_B <- 1000 # default: 1000 
rotsparam_K <- NULL # default: NULL


## reading cmd arguments
read_args <- function(args) {
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
return(list(normalization_method = normalization_method, min_peptides = min_peptides))
}

# Read experimental design 
read_expdesign <- function(file) {
  exp_design <- read.csv(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
for (col in c("raw_file", "exp_condition", "biorep", "techrep", "fraction")) {
  if (!any(col %in% colnames(exp_design))) {
    stop(paste("Missing column in experimental design file:", col))
  }
}

  if (length(unique(exp_design$biorep)) == 1) {
    exp_design$biorep <- exp_design$techrep
    exp_design$techrep <- 1
  }
  exp_design$biorep <- as.character(exp_design$biorep)
 exp_design$run <- paste0("abundance_", exp_design$exp_condition, "_", exp_design$biorep)

  return(exp_design)
}


# Read protein data
read_proteins <- function(file) {
  D <- read.csv(file)
  ab_cols <- grepl("^abundance_", colnames(D))
  # remove rows with only missing values
  D <- D[rowSums(is.na(D[, ab_cols])) < sum(ab_cols),]
  intensities <- D[, ab_cols]
  rownames(D) <- D$protein_group
  D[, ab_cols] <- intensities
  return(D)
}


# Read peptide data
read_peptides <- function(file) {
  D <- read.csv(file)
  ab_cols <- grepl("^abundance_", colnames(D))
  # remove rows with only missing values
  D <- D[rowSums(is.na(D[, ab_cols])) < sum(ab_cols),]
  intensities <- D[, ab_cols]
  rownames(D) <- D$modified_peptide
  intensities <- log2(intensities)
  D[,ab_cols] <- intensities
  return(D)
}

# Match experimental design with protein and peptide data
check_names <- function(exp_design, data) {
    all_names <- make.names(unique(paste("abundance", exp_design$exp_condition, exp_design$biorep, exp_design$techrep, sep = "_")))
    # each name in the data should be once in all_names
    for (n in all_names) {
        if (!any(n == colnames(data))) {
            stop(paste("Missing column in data:", n))
        } else if (sum(n == colnames(data)) > 1) {
            stop(paste("Column name in data is not unique:", n))
        }            
    }
    return(all_names)
}

# Create Qfeatures object from data
make_qfeatures <- function(data, all_names, exp_annotation, is_peptides = TRUE) {

qfeat <- readQFeatures(
  table = data,
  fnames = ifelse(is_peptides,"modified_peptide","protein_group"),
  ecol = all_names,
  name = ifelse(is_peptides, "peptideRaw", "proteinRaw"),
)

colData(qfeat)$exp_condition[exp_annotation$run] <- exp_annotation$exp_condition
colData(qfeat)$biorep[exp_annotation$run] <- exp_annotation$biorep
colData(qfeat)$exp_condition <- factor(colData(qfeat)$exp_condition)
colData(qfeat)$biorep <- factor(colData(qfeat)$biorep)
qfeat <- zeroIsNA(qfeat, ifelse(is_peptides, "peptideRaw", "proteinRaw"))
if (is_peptides) 
  qfeat <- logTransform(qfeat, base = 2, i = "peptideRaw", name = "peptideLog")
  return(qfeat)
}

## Normalization
# take exp and log for sum only
normalize_qfeat <- function(qfeat, normalization_method, is_peptides = TRUE) {
  assay_in_name <- ifelse(is_peptides, "peptideLog", "proteinRaw")
  assay_name <- ifelse(is_peptides, "peptideNorm", "proteinNorm")
  if (normalization_method == "none") {
    qfeat <- addAssay(qfeat, qfeat[[assay_in_name]],
      name = assay_name
    )
    addAssayLinkOneToOne(qfeat, assay_in_name, assay_name)
  } else {
    if (normalization_method == "sum") {
      ttt <- qfeat[[assay_in_name]]
      assay(ttt) <- 2^assay(ttt)
      qfeat[[assay_in_name]] <- ttt
    }
    qfeat <- normalize(qfeat,
      i = assay_in_name,
      name = assay_name,
      method = normalization_method
    )
    if (normalization_method == "sum") {
      ttt <- qfeat[[assay_name]]
      assay(ttt) <- log2(assay(ttt))
      qfeat[[assay_name]] <- ttt
    }
  }
  return(qfeat)
}

# filter proteins for min_peptides
filter_min_peps <- function(proteins, min_peptides) {
  prots <- proteins[["proteinNorm"]]
nums <- rowData(proteins)[["proteinRaw"]]
nums <- as.matrix(nums[,grep("^number_of_peptides_", names(nums))])
assay(prots)[nums < min_peptides ] <- NA
proteins[["proteinNorm"]] <- prots
return(proteins)
}


# Average over technical replicates
average_technical_replicates <- function(data, exp_design) {
    all_names <- make.names(unique(paste("abundance", exp_design$exp_condition, exp_design$biorep, sep = "_")))
    data_out <- data.frame(row.names=rownames(data))
    for (n in all_names) {
        techreps <- grep(n, colnames(data))
            data_out[,n] <- rowMeans(data[, techreps, drop=F], na.rm = TRUE)
        
    }
    return(data_out)
}



# Remove rows with less than 2 values per group
filter_data <- function(data, groups) {
    keep <- rep(TRUE, nrow(proteins))
for (i in unique(groups)) {
    group <- groups == i
    keep <- keep & rowSums(!is.na(data[, group])) >= 2
}
data_out <- data[keep,]

    return(data_out)
}

# Run ROTS for each comparison
run_rots <- function(data, groups, comparisons) {
RES_tmp <- list()
for (i in 1:ncol(comparisons)) {
    group1 <- comparisons[1, i]
    group2 <- comparisons[2, i]
    data_tmp <- data[, groups %in% c(group1, group2)]
    groups_tmp <- groups[groups %in% c(group1, group2)]
    name <- paste(as.character(group2), as.character(group1), sep = "_vs_")
    RES <- ROTS(data_tmp, groups = as.numeric(groups_tmp), log = TRUE,
        paired = FALSE, B = rotsparam_B, K = rotsparam_K, progress = TRUE)
    data[[paste("differential_abundance_pvalue", name, sep = "_")]] <- RES$pvalue
    data[[paste("differential_abundance_qvalue", name, sep = "_")]] <- RES$FDR
    data[[paste("log_ratios", name, sep = "_")]] <- RES$logfc
}
    return(data)
}

# Create contrasts
create_contrasts <- function(data) {
  conditions <- levels(colData(data)$exp_condition)

  # Create a design matrix
design <- model.matrix(~ 0 + colData(data)$exp_condition)
colnames(design) <- conditions

  # Choose the type of contrast: all-vs-all (TODO, use optional parameter for all-vs-first)
  contrast_formulas <- combn(conditions, 2, FUN = function(x) paste(x[1], "-", x[2]), simplify = FALSE)
  contrast_names <- combn(conditions, 2, FUN = function(x) paste(x[1], "vs", x[2], sep="_"), simplify = FALSE)
  contrasts <- setNames(contrast_formulas, contrast_names)
  contrast_matrix <- makeContrasts(levels = design, contrasts = contrasts)
  rownames(contrast_matrix) <- paste0("(Intercept)exp_condition", rownames(contrast_matrix))
  return(contrast_matrix)
}

## adding the new columns to the data frame
final_assembly <- function(data, contrasts, is_peptides) {
   assay_name <- ifelse(is_peptides, "peptideNorm", "proteinNorm")
  stand_out <- as.data.frame(assay(data, assay_name))
  add_cols <- rowData(data)[[assay_name]]
  if (!is_peptides) {
  stand_out <- cbind(protein_group = add_cols$protein_group, stand_out, add_cols[, grep("^number_of_peptides_", colnames(add_cols))])
  } else {
  stand_out <- cbind(add_cols$modified_peptide, stand_out, add_cols[, grep("^number_of_psms_", colnames(add_cols))])
  }
  for (i in colnames(contrasts)) {
    ttt <- add_cols[, i]
    cname <- sub(" - ", "_vs_", i)
    colnames(ttt) <- paste0(c("log_ratios_", "standard_error_", "degrees_freedom_",
                              "t_", "differential_abundance_pvalue_", "differential_abundance_qvalue_"), cname)
    stand_out <- cbind(stand_out, ttt[,c(1,5,6)])
  }
  return(stand_out)
}

