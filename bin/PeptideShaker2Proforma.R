library(matrixStats)
library(stringi)

# Reading files
ions <- read.delim("peptideshaker_filtered_out.txt", sep = "\t", check.names = F)
peptides <- read.delim("peptideshaker_peptides_out.txt", sep = "\t", check.names = F)
proteins <- read.delim("peptideshaker_proteins_out.txt", sep = "\t", check.names = F)
ptmmaptable <- read.delim("ptm_mapping.txt", sep = "\t", check.names = F)

ptmmaptable <- ptmmaptable[, c("unimod_title", "searchgui_name")]
# Take only the first entry for each searchgui_name
ptmmaptable <- ptmmaptable[!duplicated(ptmmaptable$searchgui_name), ]
ptmmaptable <- ptmmaptable[!is.na(ptmmaptable$searchgui_name), ]
ptmmapping <- ptmmaptable[, "unimod_title"]
names(ptmmapping) <- ptmmaptable$searchgui_name


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
        x <- gsub(" $", "", x)
        # Split by " of "
        x <- ptmmapping[x]
        # print(x)o
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
    return(mods)
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
  print(modified_sequence)

  return(modified_sequence)
}

peptides$"Modified Sequence" <- modify_sequence(peptides$"Variable Modifications", peptides$"Fixed Modifications", peptides$Sequence)
ions$"Modified Sequence" <- modify_sequence(ions$"Variable Modifications", ions$"Fixed Modifications", ions$Sequence)

# Not needed: Reduce protein accessions from long format (e.g. "sp|P12345|A1BG_HUMAN;sp|P12346|A1BG_HUMAN") to a string of only the accession numbers
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

write.table(peptides, "peptides_proforma.txt", row.names = F, sep = "\t", quote = F)
write.table(ions, "psms_proforma.txt", row.names = F, sep = "\t", quote = F)
write.table(proteins, "proteins_proforma.txt", row.names = F, sep = "\t", quote = F)
