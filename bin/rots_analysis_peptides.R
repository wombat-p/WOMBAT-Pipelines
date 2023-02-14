
### installation of R packages if necessary

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("ROTS")

library(limma)
library(ROTS)

rotsparam_B <- 1000 # default: 1000
rotsparam_K <- NULL # default: NULL

################################################################################
### Peptides:

### data input
D <- read.csv("stand_pep_quant_merged_pre.csv")

### columns of abundances values
ab_cols <- grepl("^abundance_", colnames(D))

### remove rows with only missing values
D <- D[rowSums(is.na(D[ab_cols])) < sum(ab_cols),]

pep_sequence <- D$modified_peptide
intensities <- D[, ab_cols]

### extract group information
group <- limma::strsplit2(colnames(intensities), "_")[,2]
group <- factor(group)
nr_groups <- length(levels(group))

RES_complete <- NULL

### statistical test with ROTS package for each combination of groups
for (i in 1:(nr_groups-1)) {
  for(j in (i+1):nr_groups) {
    print(paste("calculating ", i, ":", j, " of ", (nr_groups-1), ":", nr_groups))

    i_level <- levels(group)[i]
    j_level <- levels(group)[j]

    ### extract data from relevant groups
    group_tmp <- droplevels(group[group %in% c(i_level,j_level)])
    D_tmp <- intensities[, group %in% c(i_level,j_level)]

    ### ROTS only works when you have at least 2 valid values in each group!
    id_keep <- apply(D_tmp, 1, function(x) {
      group1 <- x[group_tmp == i_level]
      group2 <- x[group_tmp == j_level]

      nr_non_missing <- c(sum(!is.na(group1)), sum(!is.na(group2)))
      if (any(nr_non_missing < 2)) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    })

    ### remove rows where at least one group has < 2 valid values
    D_tmp2 <- D_tmp[id_keep,]
    pep_sequence2 <- pep_sequence[id_keep]

    ### ROTS test
    RES <- ROTS(D_tmp2, groups = as.numeric(group_tmp), log = TRUE,
        paired = FALSE, B = rotsparam_B, K = rotsparam_K, progress = TRUE)

    RES2 <- data.frame(pep_sequence = pep_sequence2, differential_abundance_pvalue = RES$pvalue, differential_abundance_qvalue = RES$FDR, log_ratios = RES$logfc)
    rm(RES)
    colnames(RES2)[2:4] <- paste0(colnames(RES2)[2:4], "_", i_level, "_vs_", j_level)

    if (is.null(RES_complete)) {
      RES_complete <- RES2
    } else {
      RES_complete <- merge(RES_complete, RES2, by = "pep_sequence", all = TRUE)
    }

    # clean a bit, the peptide data uses up a lot of RAM
    rm(RES2)
    gc()
  }
}

# put all data together
RES_complete <- merge(D, RES_complete, by.x = "modified_peptide", by.y = "pep_sequence", all = TRUE)
colnames(RES_complete)[1] <- "modified_peptide"

write.csv(RES_complete, "stand_pep_quant_merged.csv", row.names = FALSE)
