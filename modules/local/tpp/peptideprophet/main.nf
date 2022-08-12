process PEPTIDEPROPHET {
label 'process_medium'

conda (params.enable_conda ? "bioconda::tpp-5.0.0-pl5.22.0_0" : null)
if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://spctools/tpp:version6.1"
} else {
        container "spctools/tpp:version6.1"
        containerOptions = '-u $(id -u):$(id -g)'
}
  
  publishDir "${params.outdir}/tpp", mode:'copy'
  
  input:
  path pepxml_file
  path fasta
  path mzml
  val parameters
  
  output:
  path("${pepxml_file.baseName}.interact.pep.xml"), emit: peptideprophet
  
  script:
  enzymemap = ["Trypsin": "", "Trypsin/P": "", "Lys_C": "-eN", "Lys_N": "-eL", "Arg_C": "-eN", "Asp_N": "-eA", "CNBr": "-eM", "Glu_C": "-eG", "PepsinA": "-eN", "Chymotrypsin": "-eC", "Unspecified": "-eN"]
  enzyme = enzymemap[parameters.enzyme]
  
  """
  cp ${pepxml_file} t
  cp t ${pepxml_file}
  xinteract -N"${pepxml_file.baseName}.interact.pep.xml" -p"${parameters.ident_fdr_psm}" ${enzyme} -l"${parameters.min_peptide_length}" -THREADS=${task.cpus} -PPM -O -D"${fasta}" "${pepxml_file}"
  """

  

  }    
