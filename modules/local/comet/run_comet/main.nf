process RUN_COMET {
    label 'process_high'

    conda params.enable_conda ? "bioconda::comet-ms-2021010-h87f3376_1" : null
    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
        ? 'docker://quay.io/biocontainers/comet-ms:2021010--h87f3376_1'
        : 'quay.io/biocontainers/comet-ms:2021010--h87f3376_1'}"
    containerOptions "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
        ? null
        : '-u $(id -u):$(id -g)'}"

    publishDir "${params.outdir}/comet", mode: 'copy'

    input:
    path mzmlfile
    path fasta
    path comet_param

    output:
    path "${mzmlfile.baseName}.pep.xml", emit: comet

    script:
    """
    # set to the allowed number of CPUs for this process
    sed -i 's/^num_threads.*/'"num_threads = ${task.cpus}"/ ${comet_param}

    comet -P"${comet_param}" -N"${mzmlfile.baseName}" -D"${fasta}" "${mzmlfile}"
    """
}
