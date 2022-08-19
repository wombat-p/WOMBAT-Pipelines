process ROTS {
label 'process_medium'
  label 'process_single_thread'
  conda (params.enable_conda ? "bioconda::bioconductor-rots::1.22.0" : null)
if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "quay.io/biocontainers/bioconductor-rots:1.22.0--r41h619a076_1"
} else {
        container "quay.io/biocontainers/bioconductor-rots:1.22.0--r41h619a076_1"
}
  
publishDir "${params.outdir}/rots", mode:'copy'
 
 
  when:
  params.run_statistics

  input:
  path protein_quants
  path peptide_quants

  output:
  path "stand_prot_quant_merged.csv", includeInputs: true, emit: protein_quants_rots
  path "stand_pep_quant_merged.csv", includeInputs: true, emit: peptide_quants_rots

  script:
  """
  R CMD BATCH $baseDir/scripts/rots_analysis_proteins.R
  R CMD BATCH $baseDir/scripts/rots_analysis_peptides.R
  """
}    
