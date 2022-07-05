//
// Run maxquant and normalyzer
//

include { RAW2MZDB }                 from '../../modules/local/raw2mzdb/main'  
include { MZDB2MGF }                      from '../../modules/local/mzdb2mgf/main'
include { CREATE_DECOY_DATABASE }                      from '../../modules/local/searchgui/create_decoy_database/main'
include { PREPARE_SEARCHGUI }                      from '../../modules/local/searchgui/prepare_searchgui/main'
/*include { RUN_SEARCHGUI }                     from '../../modules/local/searchgui/run_searchgui/main'
include { CONFIG_PROLINE }                    from '../../modules/local/proline/config_proline/main'
include { EXP_DESIGN_PROLINE}                   from '../../modules/local/proline/exp_design_proline/main'
include { RUN_PROLINE }                      from '../../modules/local/proline/run_proline/main'
include { POLYSTEST }                     from '../../modules/local/polytest/main'
 */
workflow PROLINE {
    take:
    fasta // fasta file
    raws // raw files
    parameters // map of parameters 
    exp_design // experimental design file


    main:
    RAW2MZDB ( raws )
    MZDB2MGF ( RAW2MZDB.out )
    CREATE_DECOY_DATABASE ( fasta , parameters["add_decoys"])
    //PREPARE_SEARCHGUI ( parameters, fasta )
    /*RUN_SEARCHGUI ( MZDB2MGF.out, PREPARE_SEARCHGUI.out )
    CONFIG_PROLINE ( RUN_SEARCHGUI.out.collect() )
    EXP_DESIGN_PROLINE ( RUN_SEARCHGUI.out.collect() , exp_design )
    RUN_PROLINE ( RUN_SEARCHGUI.out.collect(), RAW2MZDB.out.collect(), CONFIG_PROLINE.out,  EXP_DESIGN_PROLINE.out )
    POLYSTEST ( EXP_DESIGN_PROLINE.out, RUN_PROLINE.out.collect() )

    emit:
    POLYSTEST.out.exp_design
    POLYSTEST.out.std_peps
    POLYSTEST.out.std_prots */
}
