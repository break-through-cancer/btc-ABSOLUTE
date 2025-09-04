include { PREPROCESS_FUNCOTATOR; PREPROCESS_ICONICC } from '../../modules/local/preprocess'

workflow PREPROCESS {
    take:
    sample_map

    main:
    sample_map
        .map { sample, maf, _segfile, _processed_counts ->
            [sample, maf]
        }
        .set { funcotator_files }

    sample_map
        .map { sample, _maf, segfile, processed_counts ->
            [sample, segfile, processed_counts]
        }
        .set { iconicc_files }

    PREPROCESS_FUNCOTATOR(
        funcotator_files,
        params.model_colnames
    )

    PREPROCESS_ICONICC(
        iconicc_files,
        params.capseg_rscript
    )
}