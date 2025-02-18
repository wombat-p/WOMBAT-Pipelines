process MAXQUANT_LFQ {
//    tag "$meta.id"
    publishDir "${params.outdir}/maxquant"
    
    label 'process_high'
    conda (params.enable_conda ? "bioconda::maxquant==" + params.MAXQUANT.version : null)
    if (workflow.containerEngine == 'singularity'|| workflow.containerEngine == 'apptainer') {
        container "docker://" + params.MAXQUANT.docker_container
    } else {
        container params.MAXQUANT.docker_container
    }

    input:
    path fasta
    path paramfile
    path raw

    output:
//    tuple val(meta), path("*.txt"), emit: maxquant_txt
    path "*.txt"	, emit: maxquant_txt
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 //   def args = task.ext.args ?: ''
  //  def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            maxquant: \$(maxquant --version 2>&1 > /dev/null | cut -f2 -d\" \")
    END_VERSIONS
    sed \"s_<numThreads>.*_<numThreads>$task.cpus</numThreads>_\" ${paramfile} > mqpar_changed.xml
    sed -i \"s|PLACEHOLDER|\$PWD/|g\" mqpar_changed.xml
    mkdir temp
    chmod -R a+rw *
    maxquant mqpar_changed.xml
    mv combined/txt/*.txt .
    mv combined/proc/*unningTimes.txt runningTimes.txt
    """
}
