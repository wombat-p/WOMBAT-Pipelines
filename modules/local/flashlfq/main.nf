process FLASHLFQ {
  label 'process_high'
  conda (params.enable_conda ? "bioconda::flashlfq-1.1.1" : null)
if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://quay.io/biocontainers/flashlfq:1.1.1--hdfd78af_1"
} else {
        container "quay.io/biocontainers/flashlfq:1.1.1--hdfd78af_1"
}
  
publishDir "${params.outdir}/flashlfq", mode:'copy'
 
  input:
  path peptideshaker_out
  path mzmlfiles
  val parameters
  
  output:
  path "QuantifiedPeaks.tsv", emit: flashlfq_peaks
  path "QuantifiedPeptides.tsv", emit: flashlfq_peptides
  path "QuantifiedProteins.tsv", emit: flashlfq_proteins
  
  script:
  """
  first_line=""
  for file in *.txt
  do
    echo \$file
    tail -n +2 "\$file" >> tlfq_ident.tabular
    first_line=\$(head -n1 "\$file")
  done
  echo "\$first_line" | cat - tlfq_ident.tabular > lfq_ident.tabular
  FlashLFQ --idt "lfq_ident.tabular" --rep "./" --out ./ --mbr ${parameters.match_between_runs} --ppm ${parameters.precursor_mass_tolerance} --thr ${task.cpus}
  """
 
}    
