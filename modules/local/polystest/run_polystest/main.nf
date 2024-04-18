process POLYSTEST {
  label 'process_high'
  // TODO: REVERT TO OFFICIAL PACKAGE
  conda (params.enable_conda ? "bioconda::polystest-1.5.2" : null)
  if (workflow.containerEngine == 'singularity'|| workflow.containerEngine == 'apptainer') {
        container "docker://quay.io/biocontainers/polystest:1.5.01--hdfd78af_0"
  } else {
// TODO        container "quay.io/biocontainers/polystest:1.5.01--hdfd78af_0"
        container "veitveit/polystest"
  }
  
  publishDir "${params.outdir}/polystest", mode:'copy'
  
  when:
  parameters.run_statistics

  input:
  path exp_design
  path peptides
  path proteins
  path pep_param
  path prot_param
  val parameters
  
  output:
    path "stand_prot_quant_merged.csv", emit: stdprotquant
    path "stand_pep_quant_merged.csv", emit: stdpepquant
    path "exp_design.txt", emit: exp_design
  
  script:
  """
  
  sed -i "s/threads: 2/threads: ${task.cpus}/g" "${pep_param}"
  sed -i "s/threads: 2/threads: ${task.cpus}/g""${prot_param}"
  
  runPolySTestCLI.R pep_param.yml
  runPolySTestCLI.R prot_param.yml

  # Convert to standard format
  cp "${exp_design}" exp_design.txt
  Rscript $baseDir/bin/Convert2StandFormat.R


  """

}
