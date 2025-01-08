process RAW2MZML {
    label 'process_intermediate'
    
    conda params.enable_conda ? "bioconda::thermorawfileparser==1.4.5--h05cac1d_1" : null

    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
        ? 'docker://quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1'
        : 'quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1'}"
    
    publishDir "${params.outdir}/mzml", mode:'copy'
    
    input:
    path rawfile
    
    output:
    path "${rawfile.baseName}.mzML", includeInputs:true , emit: mzmls
    
    script:
    """
    # Check if the file is a mzML file
    if [[ "${rawfile}" =~ .*\\.(mzML|mzml)\$ ]]
    then
        # check if same file
        if [[ "${rawfile}" != "${rawfile.baseName}.mzML" ]]
        then
            cp "${rawfile}" "${rawfile.baseName}.mzML"
        fi
    else
        thermorawfileparser -i "${rawfile}" -b "${rawfile.baseName}.mzML" -f 2
    fi
    """
}
