process PREPARE_SEARCHGUI {
  label 'process_low'
  label 'process_single_thread'
    conda (params.enable_conda ? "bioconda::searchgui-4.0.41" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://quay.io/biocontainers/searchgui:4.0.41--h779adbc_1"
    } else {
        container "quay.io/biocontainers/searchgui:4.0.41--h779adbc_1"
    }
  publishDir "${params.outdir}/searchgui_params", mode:'copy'
  
  input:
  path fasta_decoy
  val parameters
  
  output:
  path "searchgui.par", emit: searchgui_param
  
  script:
  """
  searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI -prec_tol ${params.precursor_mass_tolerance} \\
    -frag_tol ${params.fragment_mass_tolerance} -enzyme "${parameters["enzyme"]}" -mc ${parameters["allowed_miscleavages"]} -max_isotope ${parameters["isotope_error_range"]} \\
    -fixed_mods "${parameters["fixed_mods"]}" -variable_mods "${parameters["variable_mods"]}" -min_charge ${params.min_charge} -max_charge ${params.max_charge} \\
    -fi ${parameters["fions"]} -ri ${parameters["rions"]} -xtandem_quick_acetyl 0 -xtandem_quick_pyro 0 -peptide_fdr ${parameters["ident_fdr_peptide"]}\\
    -protein_fdr ${parameters["ident_fdr_protein"]} -psm_fdr ${parameters["ident_fdr_psm"]} -psm_fdr ${parameters["max_mods"]}\\
    # only available as search engine specific parameters
    -myrimatch_num_ptms ${parameters["max_mods"] -ms_amanda_max_mod ${parameters["max_mods"] -msgf_num_ptms ${parameters["max_mods"] -meta_morpheus_max_mods_for_peptide \\
    ${parameters["max_mods"] -directag_max_var_mods ${parameters["max_mods"]} -comet_num_ptms ${parameters["max_mods"] -tide_max_ptms ${parameters["max_mods"]}  \\
    -db ${fasta_decoy}  -out searchgui.par
  """    
}
  
