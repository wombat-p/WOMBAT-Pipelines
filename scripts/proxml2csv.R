######################################################################
# AUTH  Jeff Jones, extended by Veit Schwaemmle
# DATE  2022.07.23
# DESC  parse protXML post TPP: ProteinProphet
#
######################################################################

rm(list=ls())
require(XML)
library(progress)
library(stringi)
options(warn=-1)

help_text <- "
 NAME
    pepxml2csv.R

 SYNOPSIS
    pepxml2csv.R --xml=<path_pepxml> --fdr=0.1 --npep=2

 DESCRIPTION
    extract peptide id data from PeptideProphet pepXML files for downstream analysis

 COMMAND LINE

    --xml <path_pepxml>

    --csv <path_csv> (optional: path_pepxml.csv)

    --fdr false discovery rate value between 0-1 (default: 0.1 ~ 10% FDR) 

    --npep value >= 1 (default: 2) assumes user avoids 1-hit proteins

 EXAMPLE

    Rscript proxml2csv.R --xml=<path_to.proXML>

"

###############################################################################
# USER INPUT
path_xml                      <- NULL
path_csv                      <- NULL
fdr_cutf                      <- 0.1
min_npep                      <- 2

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--xml", arg) ) path_xml <- arg_value
    if( grepl("--csv", arg) ) path_csv <- arg_value
    if( grepl("--fdr", arg) ) fdr_cutf <- arg_value
    if( grepl("--npep", arg) ) min_npep <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

###############################################################################
# INPUT VALIDATION
message <- NULL
if(is.null(path_xml)) message <- stop("ERROR\n", "  no mzXML file declared\n")
if(!grepl("pro[t.]*xml$", tolower(path_xml))) message <- paste0(message, "  xml file (--xml) not a supported format\n")
if(is.null(path_csv))
    path_csv = paste0(path_xml, ".csv")
if(!grepl(".csv$", path_csv)) message <- paste0(message, "  csv file (--csv) not a supported format\n")

if(!is.null(message)) stop("ERROR\n", message)

cat("pepXML to CSV started\n")
cat(" xml file:                        ", path_xml, "\n")
cat(" csv file:                        ", path_csv, "\n")
cat(" fdr cut off:                     ", fdr_cutf, "\n")
cat(" min n peptides:                  ", min_npep, "\n")

#
# Read in the data
#
cat(" reading xml file ...")
data <- xmlParse(path_xml)
file_xml <- sub('.*/', '', path_xml)

#
# convert to a master list
#
xml_data <- xmlToList(data)
cat("\n")


#
# probe for model performance values
#
roc_dat <- xml_data$protein_summary_header$program_details$proteinprophet_details

if(is.null(roc_dat)) stop("ERROR\n", "  xml does not contain FDR statistics\n")

w_dp <- which(names(roc_dat) == 'protein_summary_data_filter')
w_ep <- which(names(roc_dat) == 'error_point')

roc_df <- do.call(rbind, lapply(lapply(roc_dat[w_dp], unlist), "[",
                                unique(unlist(c(sapply(roc_dat[w_dp],names))))))
row.names(roc_df) <- 1:dim(roc_df)[1]
roc_df <- as.data.frame(roc_df)
roc_df$analysts <- 'roc_data_point'
colnames(roc_df)[c(1,3,4,5)] <- c('min_prob','error','num_corr','num_incorr')

err_df <- do.call(rbind, lapply(lapply(roc_dat[w_ep], unlist), "[",
                                unique(unlist(c(sapply(roc_dat[w_ep],names))))))
row.names(err_df) <- 1:dim(err_df)[1]
err_df <- as.data.frame(err_df)
err_df$analysts <- 'error_point'
err_df$sensitivity <- NA

pdf <- rbind(roc_df, err_df)
row.names(pdf) <- 1:dim(pdf)[1]

pdf$min_prob <- as.numeric(as.character(pdf$min_prob))
pdf$sensitivity <- as.numeric(as.character(pdf$sensitivity))
pdf$error <- as.numeric(as.character(pdf$error))
pdf$num_corr <- as.numeric(as.character(pdf$num_corr))
pdf$num_incorr <- as.numeric(as.character(pdf$num_incorr))

w_pcf <- which(pdf$error == max(pdf[pdf$error <= fdr_cutf &
                                        pdf$analysts == 'error_point',]$error) &
                   pdf$analysts == 'error_point')
min_prob <- pdf[w_pcf,]$min_prob
cat("\n")
print(pdf[w_pcf,], row.names = FALSE)
cat("\n")

#
# General definitions
#
substitute_names <- c(SI="Spectral index", SIn="Normalized spectral index", counts="number of quantified peptides", 
                      NSAF="Normalized spectral abundance factor")


#
# traverse the master list
#
w_protg <- which(names(xml_data) == 'protein_group')
pb <- progress_bar$new(
    format = " extracting data [:bar] :percent eta: :eta",
    total = length(w_protg), clear = FALSE, width= 60)
hits <- list()
hit_i <- 0
for ( prot_i in w_protg ){
    
    pb$tick()
    
    prot <- xml_data[prot_i]
    
    
    #
    # grab the protein attributes
    #
    prot_att <- append(unlist(prot$protein_group$.attrs), file_xml)
    names(prot_att)[length(prot_att)] <- 'fileName'
    w_proteins <- which(names(prot$protein_group) == 'protein')
    
    for(result in prot$protein_group[w_proteins]){
        
        #
        # grab the result attributes
        # 
        
        result_att <- unlist(result$.attrs)

        w_peptides <- which(names(result) == 'peptide')
        
        #
        # add sequence+charge to find same peptides in stpeter fields
        #
        names(w_peptides) <- sapply(w_peptides, function(x) paste0(unlist(result[x]$peptide$.attr)[c("peptide_sequence","charge")], collapse=""))
        
        # 
        # get results from StPeter analysis (if available)
        # 
        
        a_peptides <- a_prot <- NULL
        if (!is.null(result[['analysis_result']])) {
            a_prot <- result[['analysis_result']]$StPeterQuant
            a_peptides <- which(names(a_prot) == 'StPeterQuant_peptide')
            names(a_peptides) <- sapply(a_peptides, function(x) paste0(a_prot[x]$StPeterQuant_peptide[c("sequence","charge")], collapse=""))
            a_prot$.attrs <- a_prot$.attrs[names(substitute_names)]
            names(a_prot$.attrs) <- substitute_names
        }
        
        

        for(peptide in result[w_peptides]){
            
            peptide_vals <- unlist(peptide$.attrs)
            modified_peptide <- peptide_vals['peptide_sequence']
            names(modified_peptide) <- NULL
            if (!is.null(peptide$modification_info)) {
                w_mods <- which(names(peptide$modification_info) == 'mod_aminoacid_mass')
                mod_info <- NULL
                for (mod in 1:length(w_mods)) {
                    mod_info <- rbind(mod_info, unlist(peptide$modification_info[mod]))
                }
                if (length(unlist(mod_info)) > 1) {
                modified_peptide <- unlist(stri_sub_replace_all(peptide_vals['peptide_sequence'], from=as.numeric(mod_info[,1])+1,
                                     to=as.numeric(mod_info[,1]), 
                                     replacement=paste0("[",round(as.numeric(mod_info[,2])),"]")))
                }
            }
            peptide_vals <- c(peptide_vals, modified_peptide=modified_peptide)


            #
            # append all the values to our new list
            #
            
            hit_i <- hit_i + 1
            
            hits[[hit_i]] <- list()
            hits[[hit_i]] <- append(hits[[hit_i]], as.list(prot_att))
            hits[[hit_i]] <- append(hits[[hit_i]], as.list(result_att))
            hits[[hit_i]] <- append(hits[[hit_i]], as.list(peptide_vals))
 
            if (!is.null(a_peptides)) {
                hits[[hit_i]] <- append(hits[[hit_i]], as.list(unlist(a_prot$.attrs)))
                pep_pos <- which(names(a_peptides) == paste0(peptide_vals[c("modified_peptide","charge")],collapse=""))
                if (length(pep_pos) > 0) {
                    hits[[hit_i]] <- append(hits[[hit_i]], as.list(unlist(a_prot[a_peptides[pep_pos]])))
                }
                
            }
        }
        
        
    }
}


#
# use lapply to converge the hits list into a matrix with named columns
#
data <- do.call(rbind, lapply(lapply(hits, unlist), "[",
                              unique(unlist(c(sapply(hits,names))))))

#
# convert the matrix to a data frame
#
data <- as.data.frame(data)
colnames(data) <- unique(unlist(c(sapply(hits,names))))
data$probability <- as.numeric(as.character(data$probability))
data$total_number_distinct_peptides <- as.numeric(as.character(data$total_number_distinct_peptides))
#data$percent_coverage <- as.numeric(as.character(data$percent_coverage))

#
# filter data by cutoff min_prob and rank
#
data <- data[data$probability >= min_prob & data$total_number_distinct_peptides >= min_npep & !is.na(data$group_number),]

pdata <- unique(data[,c('group_number','probability','protein_name',
                        'n_indistinguishable_proteins',#'percent_coverage',
                        'total_number_peptides','total_number_distinct_peptides',
                        substitute_names)])

#pdata <- pdata[order(-pdata$percent_coverage),]

#print(head(pdata[,c(3,5,7)]), row.names = FALSE)

write.csv(data, path_csv, row.names = FALSE)
cat("\nFINISHED:\n")
cat(" data written to                  ", path_csv, "\n")

path_csv <- sub(".csv", "_summary.csv", path_csv)
write.csv(pdata, path_csv, row.names = FALSE)
cat("\nFINISHED:\n")
cat(" protein summary written to       ", path_csv, "\n")

path_csv <- sub(".csv", "_performance.csv", path_csv)
write.csv(pdf, path_csv, row.names = FALSE)
cat(" performance roc written to       ", path_csv, "\n")