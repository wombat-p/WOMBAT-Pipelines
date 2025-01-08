//
// Get assets like general parameter files
//

//
// Run Transproteomic Pipeline (TPP) with ROTS statistical analysis
//

include { RAW2MZML }       from '../../modules/local/raw2mzml/main'  
include { WRITE_CONFIG }   from '../../modules/local/comet/write_config/main'
include { RUN_COMET }      from '../../modules/local/comet/run_comet/main'
include { PEPTIDEPROPHET } from '../../modules/local/tpp/peptideprophet/main'
include { PROTEINPROPHET } from '../../modules/local/tpp/proteinprophet/main'
include { STPETER }        from '../../modules/local/tpp/stpeter/main'
include { PROTXML2CSV }    from '../../modules/local/tpp/protxml2csv/main'
include { MERGEOUTPUT }    from '../../modules/local/tpp/mergeoutput/main'
include { ROTS }           from '../../modules/local/rots/main'

workflow TPP {
    take:
    fasta // fasta file
    raws // raw files
    parameters // map of parameters
    exp_design // experimental design file
    ptm_mapping // map to convert from unimod to searchgui

    main:
    mzmls               = RAW2MZML(raws)
    comet_params        = WRITE_CONFIG(parameters, ptm_mapping)
    comet_results       = RUN_COMET(mzmls, fasta, comet_params)
    pepprophet_results  = PEPTIDEPROPHET(comet_results, fasta, parameters)
    protprophet_results = PROTEINPROPHET(pepprophet_results, parameters)
    stpeter_results     = STPETER(protprophet_results.proteinprophet_xml, parameters, mzmls.collect(), fasta)
    protquants          = PROTXML2CSV(stpeter_results, parameters)
    merged_quants       = MERGEOUTPUT(protquants.collect(), raws.collect(), exp_design)

    rots_results_protein = null
    rots_results_peptide = null
    if (parameters.run_statistics) {
        rots_results = ROTS(merged_quants.stdprotquant_qc, merged_quants.stdpepquant_qc, parameters )
        rots_results_protein = rots_results.protein_quants_rots
        rots_results_peptide = rots_results.peptide_quants_rots
    }
    
    emit:
    expdesign = merged_quants.expdesign
    rots_results_protein
    rots_results_peptide
}
