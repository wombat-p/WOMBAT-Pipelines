import groovy.json.JsonOutput

  process CONVERT_PROFORMA {
    label 'process_low'
    label 'process_single_thread'
    publishDir "${params.outdir}", mode:'copy'
    conda (params.enable_conda ? "conda-forge::notyetavailable" : null)
    if (workflow.containerEngine == 'singularity'|| workflow.containerEngine == 'apptainer') {
        container "docker://wombatp/maxquant-pipeline:v0.2"
    } else {
        container "wombatp/maxquant-pipeline:v0.2"
    }


  input:
   path peptideshaker_filtered_out
   path peptideshaker_peptides_out
   path peptideshaker_proteins_out
 
  output:
   path "psm_proforma.txt"    , emit: peptideshaker_proforma_filtered
   path "pep_proforma.txt"    , emit: peptideshaker_proforma_peptides
   path "proteins_proforma.txt"    , emit: peptideshaker_proforma_proteins
  
  script:
  """
  cp "${peptideshaker_out}" peptideshaker_out.txt
  Rscript $baseDir/bin/PeptideShaker2Proforma.R
  """
}


