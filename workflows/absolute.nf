include { INPUT_CHECK }         from '../subworkflows/local/input_check'
include { PHASE_I }         from '../subworkflows/local/absolute'
include { PHASE_II }         from '../subworkflows/local/absolute'

workflow ABSOLUTE_WORKFLOW {
    if (!params.samplesheet) {
        exit 1, 'Samplesheet not specified. Please, provide a --samplesheet=/path/to/samplesheet.csv !'
    }
    file(params.samplesheet, checkIfExists: true)

    if (!params.outdir) {
        exit 1, 'Output directory not specified. Please, provide a --outdir=/path/to/outdir !'
    }

    if (!params.phase) {
        exit 1, 'Phase not specified. Please specify --phase 1 or 2.'
    }

    INPUT_CHECK( file(params.samplesheet) )

    if ("${params.phase}" == "1") {
        PHASE_I( INPUT_CHECK.out.sample_map )
    } else if ("${params.phase}" == "2") {
        PHASE_II( INPUT_CHECK.out.sample_map )
    } else {
        error "Invalid phase: choose --phase 1 or --phase 2"
    }
}
