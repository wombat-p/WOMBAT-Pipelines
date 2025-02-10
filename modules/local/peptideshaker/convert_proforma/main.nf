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
   path peptideshaker_out
 
  output:
   path "psm_proforma.txt"    , emit: peptideshaker_proforma_file_filtered
  
  script:
  """
  cp "${peptideshaker_out}" peptideshaker_out.txt
  Rscript $baseDir/bin/PeptideShaker2Proforma.R
  """
}


