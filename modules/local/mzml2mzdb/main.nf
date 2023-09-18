process MZML2MZDB {    
    label 'process_low'
    label 'process_single_thread'
    conda (params.enable_conda ? "conda-forge::python-3.8.3" : null)
    if (workflow.containerEngine == 'singularity'|| workflow.containerEngine == 'apptainer') {
        container "docker://wombatp/proline-pipeline:v0.21"
    } else {
        container "wombatp/proline-pipeline:v0.21"
    }

  publishDir "${params.outdir}/mzdb", mode:'copy'

  input:
  path mzmlfile
  
  output:
  path "${mzmlfile.baseName}.mzDB" , emit: mzdbs
  
  script:
  """
     mzML2mzDB "${mzmlfile}"
     mv "${mzmlfile}.mzDB" "${mzmlfile.baseName}.mzDB"
  """
}
