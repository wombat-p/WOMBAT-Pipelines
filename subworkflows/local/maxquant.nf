//
// Run maxquant and normalyzer
//

include { CONVERT_MAXQUANT }                 from '../../modules/local/sdrfpipelines/convert_maxquant/main'  
include { MAXQUANT_LFQ }                      from '../../modules/nf-core/modules/maxquant/lfq/main'
include { NORMALYZERDE }                  from '../../modules/local/normalyzerde/main'

workflow MAXQUANT {
    take:
    sdrf_local // sdrf with parameter values 
    fasta // fasta file
    raws // raw files


    main:
    CONVERT_MAXQUANT (sdrf_local, fasta)
    MAXQUANT_LFQ ( fasta, CONVERT_MAXQUANT.out.maxquantpar, raws )
    NORMALYZERDE (MAXQUANT_LFQ.out.maxquant_txt, CONVERT_MAXQUANT.out.exp_design, 
                  sdrf_local.splitCsv(sep: "\t", header: true).map{row -> row["comment[normalization method]"]}.unique())

    emit:
    NORMALYZERDE.out.exp_design
    NORMALYZERDE.out.std_prots
    NORMALYZERDE.out.std_peps
 
}
