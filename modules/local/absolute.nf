process ABSOLUTE_RUN_I {
    container "ghcr.io/btc-absolute:latest"

    input:
    tuple val(sample),
        path(seg_data),
        path(indel_data),
        path(snp_data)

    output:
    path "${sample}", emit: result_dir

    script:
    """
    Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_cli_start.R \
        --seg_dat_fn ${seg_data} \
        --maf_fn ${snp_data} \
        --indelmaf_fn ${indel_data} \
        --sample_name ${sample} \
        --results_dir . \
        --ssnv_skew 0.99 \
        --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5 \

    # Create a directory for this sample's results
    mkdir ${sample}

    # Move all relevant output files into that directory
    mv ${sample}.ABSOLUTE_* ${sample}/
    mv ${sample}.PP-* ${sample}/
    mv ${sample}.ABSOLUTE.RData ${sample}/
    """
}

process ABSOLUTE_RUN_II {
    container "ghcr.io/btc-absolute:latest"

    input:
    tuple val(sample),
        path(seg_data),
        path(indel_data),
        path(snp_data),
        val(purity),
        val(ploidy)

    output:
    path "${sample}", emit: result_dir

    script:
    """
    Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_cli_start.R \
        --seg_dat_fn ${seg_data} \
        --maf_fn ${snp_data} \
        --indelmaf_fn ${indel_data} \
        --sample_name ${sample} \
        --results_dir . \
        --ssnv_skew 0.99 \
        --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5 \
        --force_alpha ${purity} \
        --force_tau ${ploidy}
    
    # Create a directory for this sample's results
    mkdir ${sample}

    # Move all relevant output files into that directory
    mv ${sample}.ABSOLUTE_* ${sample}/
    mv ${sample}.PP-* ${sample}/
    mv ${sample}.ABSOLUTE.RData ${sample}/
    """
}

process ABSOLUTE_FORCECALL{
    container "ghcr.io/btc-absolute:latest"
    // publishDir params.outdir
    input:
    tuple val(sample),
        path(rdata)

    output:
    path "${sample}", emit: result_dir

    script:
    """
    Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_extract_cli_start.R \
        --solution_num 1 \
        --analyst_id force_called \
        --rdata_modes_fn ${rdata} \
        --sample_name ${sample} \
        --results_dir . \
        --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5
        
    mkdir -p ${sample}

    mv reviewed/samples \
        reviewed/SEG_MAF \
        reviewed/${sample}.aggregate_MAF.Rds \
        reviewed/${sample}.force_called.ABSOLUTE.table.txt \
        reviewed/${sample}.MAF_list.Rds \
        ${sample}/
    """
}