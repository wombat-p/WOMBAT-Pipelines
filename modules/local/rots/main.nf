process ROTS {
    label 'process_medium'
    label 'process_single_thread'
    
    conda params.enable_conda ? "bioconda::bioconductor-rots::1.22.0" : null
    
    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
        ? 'docker://wombatp/transproteomic-pipeline:0.24'
        : 'wombatp/transproteomic-pipeline:0.24'}"
    
    publishDir "${params.outdir}/rots", mode:'copy'
    
    input:
    path protein_quants
    path peptide_quants
    val parameters
    
    output:
    path "stand_prot_quant_merged.csv", includeInputs: true, emit: protein_quants_rots
    path "stand_pep_quant_merged.csv", includeInputs: true, emit: peptide_quants_rots
    
    script:
    """
    rots_analysis_proteins.R
    rots_analysis_peptides.R
    """
}    
