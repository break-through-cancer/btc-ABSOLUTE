nextflow.enable.dsl = 2

include { ABSOLUTE_WORKFLOW } from './workflows/absolute.nf'

workflow {
    ABSOLUTE_WORKFLOW()
} 