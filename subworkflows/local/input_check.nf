//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet

    main:

    // 1. Run samplesheet_check
    SAMPLESHEET_CHECK( samplesheet, params.phase )
        .samplesheet_utf8
        .set { samplesheet_utf8 }

    // 2. Parse samplesheet
    samplesheet_utf8
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def sample =  row.sample
            def seg_path   = file(row.seg_path)
            def indel_path = file(row.indel_path)
            def snp_path   = file(row.snp_path)
            def purity    = row.purity ? row.purity.toDouble() : null
            def ploidy    = row.ploidy ? row.ploidy.toDouble() : null
            def rdata_path = (row.rdata_path && row.rdata_path.trim()) ? file(row.rdata_path) : null
            return [sample, seg_path, indel_path, snp_path, purity, ploidy, rdata_path]
        }
        .set { sample_map }

    emit:
    sample_map          //input to sample-level analysis
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
