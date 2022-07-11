process CONVERT_POLYSTEST {
label 'process_low'

conda (params.enable_conda ? "bioconda::polystest-1.1" : null)
if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://quay.io/biocontainers/polystest:1.1--hdfd78af_2"
} else {
        container "quay.io/biocontainers/polystest:1.1--hdfd78af_2"
}
  
  publishDir "${params.outdir}/polystest", mode:'copy'
  
  input:
  path exp_design 
  path pep_quant 
  path prot_quant
  
  output:
  path "stand_prot_quant_merged.csv", emit: stdprotquant
  path "stand_pep_quant_merged.csv", emit: stdpepquant
  path "exp_design.txt", emit: exp_design
  
  when:
  params.run_statistics
  
  script:
 """
  cp "${exp_design}" exp_design.txt
  Rscript $baseDir/scripts/Convert2StandFormat.R
  """


}