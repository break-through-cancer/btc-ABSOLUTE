include { ABSOLUTE_RUN_I; ABSOLUTE_FORCECALL; ABSOLUTE_RUN_II } from '../../modules/local/absolute'

workflow PHASE_I {
    take:
    sample_map

    main:
    sample_map
        .map { sample, seg, indel, snp, _purity, _ploidy, _rdata ->
            [sample, seg, indel, snp]
        }
        .set { phase1_inputs }

    ABSOLUTE_RUN_I(
        phase1_inputs, params.ssnv_skew
    )
}

workflow PHASE_II {
    take:
    sample_map

    main:
    sample_map
        .map { sample, seg, indel, snp, purity, ploidy, _rdata ->
            [sample, seg, indel, snp, purity, ploidy]
        }
        .set { phase2_inputs }

    sample_map
        .map { sample, _seg, _indel, _snp, _purity, _ploidy, rdata ->
            [sample, rdata]
        }
        .set { forcecall_inputs }

    ABSOLUTE_RUN_II(
        phase2_inputs, params.ssnv_skew
    )

    ABSOLUTE_FORCECALL(
        forcecall_inputs
    )
}