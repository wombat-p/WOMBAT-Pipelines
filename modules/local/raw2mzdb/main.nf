process RAW2MZDB {
    label 'process_low'
    label 'process_single_thread'
    conda (params.enable_conda ? "conda-forge::python-3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://wombatp/proline-pipeline:dev"
    } else {
        container "wombatp/proline-pipeline:dev"
    }

  publishDir "${params.outdir}/mzdb", mode:'copy'

  input:
  path rawfile
  
  output:
  path "${rawfile.baseName}.mzDB" , emit: mzdbs
  
  script:
  """
  ls -la
  thermo2mzdb "${rawfile}"
  mv "${rawfile}.mzDB" "${rawfile.baseName}.mzDB" 
  """
}
