process CONVERT_PROLINE {
label 'process_low'

conda (params.enable_conda ? "bioconda::polystest-1.5.1" : null)
if (workflow.containerEngine == 'singularity'|| workflow.containerEngine == 'apptainer') {
//        container "docker://quay.io/biocontainers/polystest:1.5.1--hdfd78af_0"
        container "docker://veitveit/polystest"
} else {
//        container "quay.io/biocontainers/polystest:1.5.1--hdfd78af_0"
        container "veitveit/polystest"
}
  
  publishDir "${params.outdir}/polystest", mode:'copy'
  
  input:
  path exp_design 
  path proline_res
  
  output:
  path "petide_ions_proline.csv", emit: proline_ions
  path "peptides_proline.csv", emit: proline_peptides
  path "proteins_proline.csv", emit: proline_proteins
  path "pep_param.yml", emit: pep_param
  path "prot_param.yml", emit: prot_param
  path "exp_design.txt", emit: exp_design
  
  when:
  params.run_statistics
  
  script:
"""
  convertProline=\$(which runPolySTestCLI.R)
  
  echo \$convertProline
  convertProline=\$(dirname \$convertProline)
  
  echo \$convertProline
  Rscript \${convertProline}/convertFromProline.R ${exp_design} ${proline_res}
"""


}
