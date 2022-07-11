// loading library for json and yaml
import groovy.json.JsonOutput
import org.yaml.snakeyaml.Yaml
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowWombat.initialise(params, log)

// Check input path parameters to see if they exist
//def checkPathParamList = [ params.fasta, params.params, params.raws, params.mzmls, params.sdrf, params.exp_design]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Input fasta file not specified!' }
if (params.sdrf) { ch_sdrf = file(params.sdrf) } else { ch_sdrf = file("no_sdrf") }
if (params.exp_design) { ch_exp_design = file(params.exp_design) } else { ch_exp_design = file("no_exp_design") }
if (params.raws) { ch_raws = Channel.fromPath(params.raws) }  else { ch_raws = [file("no_raws")] }
if (params.mzmls) { ch_mzmls = Channel.fromPath(params.mzmls) }  else { ch_mzmls = [file("no_mzmls")] }
if (!params.raws && !params.mzmls && !params.sdrf) { exit 1, 'Neither raw files, mzml files nor sdrf file provided!' }
if (params.parameters) { ch_params = file(params.parameters) } else { ch_params = file("no_params") }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

ch_sdrfmapping = file("https://raw.githubusercontent.com/bigbio/proteomics-metadata-standard/master/sdrf-proteomics/assets/param2sdrf.yml", checkIfExists: true)
ch_ptm_mapping = Channel.fromPath("assets/unimod2searchgui_mapping.tsv").splitCsv(header: true, sep:"\t", quote:'\"').map{ row -> [("$row.unimod_title of $row.residue".toString()): row.searchgui_name] }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { TRANSPROTEOMICS } from '../subworkflows/local/transproteomics'
include { PROLINE } from '../subworkflows/local/proline'
include { MAXQUANT } from '../subworkflows/local/maxquant'
//include { COMPOMICS } from '../subworkflows/local/compomics'
include { PREPARE_FILES } from '../modules/local/prepare_files/main'
include { CALCBENCHMARKS } from '../modules/local/calcbenchmarks/main'
include { SDRFMERGE } from '../modules/local/sdrfpipelines/sdrfmerge/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
//include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
//def multiqc_report = []

workflow WOMBAT {

    ch_versions = Channel.empty()

    //
    // MODULE: Prepare input files
    //
    PREPARE_FILES (ch_sdrf, ch_params, ch_exp_design, ch_raws.collect(), ch_mzmls.collect(), ch_sdrfmapping)    


    //
    // MODULE: Create sdrf with data analysis parameters
    //
    SDRFMERGE (PREPARE_FILES.out.sdrf_local, PREPARE_FILES.out.params, ch_sdrfmapping)
 
    // 
    // Reading parameter from created parameter yaml file
    // 
    ch_parameters = PREPARE_FILES.out.params.map{ new Yaml().load(it)["params"] }


//    CUSTOM_DUMPSOFTWAREVERSIONS (
//        ch_versions.unique().collectFile(name: 'collated_versions.yml')
//    )

    //
    // SUBWORKFLOW 1:
    //
    // Maxquant-based
    if (params.workflow.contains("all") || params.workflow.contains("maxquant")) {
       MAXQUANT (SDRFMERGE.out.params, ch_fasta, PREPARE_FILES.out.raws.collect())

    //
    // MODULE: calculate benchmarks
    //
    CALCBENCHMARKS ( JsonOutput.prettyPrint(JsonOutput.toJson(params)), MAXQUANT.out[0], MAXQUANT.out[1], MAXQUANT.out[2], ch_fasta )
    }

    //
    // SUBWORKFLOW 1:
    //
    // Proline-based
    if (params.workflow.contains("all") || params.workflow.contains("proline")) {
        PROLINE (ch_fasta, ch_raws, ch_parameters, PREPARE_FILES.out.exp_design, ch_ptm_mapping)

    //
    // MODULE: calculate benchmarks
    //
    CALCBENCHMARKS ( JsonOutput.prettyPrint(JsonOutput.toJson(params)), PROLINE.out[0], PROLINE.out[1], PROLINE.out[2], ch_fasta )
    }



    workflow_summary    = WorkflowWombat.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    //ch_multiqc_files = Channel.empty()
    //ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    //MULTIQC (
    //    ch_multiqc_files.collect()
    //)
    //multiqc_report = MULTIQC.out.report.toList()
    //ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
