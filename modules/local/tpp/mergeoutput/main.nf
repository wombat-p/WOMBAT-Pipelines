process MERGEOUTPUT {
    label 'process_medium'

    conda params.enable_conda ? "conda-forge::r-dplyr-1.0.9-r40ha35a809_0" : null

    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
        ? 'docker://rocker/tidyverse:4.2.1'
        : 'rocker/tidyverse:4.2.1'}"
    
    publishDir "${params.outdir}/tpp", mode:'copy'
    
    input:
    path csv_files
    val raw_files
    path exp_design_file
    
    output:
    path "all_prot_quant_merged.csv" , emit: allprotquant
    path "all_pep_quant_merged.csv" , emit: allpepquant
    path "stand_prot_quant_merged_pre.csv" , emit: stdprotquant_qc
    path "stand_pep_quant_merged_pre.csv" , emit: stdpepquant_qc
    path "exp_design.txt" , emit: expdesign
    
    script:
    if (exp_design_file.getName() == "none") {
        // no experimental design file provided
        def expdesign_text = "raw_file\texp_condition"

        raw_files.each{ raw ->
            expdesign_text += "\n${raw.getName()}\texp_condition"
        }
        """
        touch exp_design.txt  
        echo "${expdesign_text}" >> exp_design.txt

        MergeTPPOutput.R
        """
    } else {
        """
        if [[ "${exp_design_file}" != "exp_design.txt" ]] 
        then
            cp "${exp_design_file}" exp_design.txt
        fi
        
        MergeTPPOutput.R
        """
    }
}
