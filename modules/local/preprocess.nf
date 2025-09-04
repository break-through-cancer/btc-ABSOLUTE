process PREPROCESS_FUNCOTATOR {
    container "ghcr.io/break-through-cancer/btc-absolute:latest"
    label 'process_low'

    input:
    tuple val(sample),
        path(maf)
    path model_colnames
    output:
    path "${sample}.snp", emit: indel_data
    path "${sample}.indel", emit: snp_data

    script:
    """
    Rscript - <<EOF
    library('tools')

    maf <- read.delim('${maf}',comment.char = '#')
    model_cols <- tolower(readLines('${model_colnames}'))

    maf_colnames <- tolower(colnames(maf))
    retain_columns <- maf_colnames %in% model_cols
    maf <- maf[,retain_columns]

    names(maf)[names(maf) == "Start_Position"] <- "Start_position"
    names(maf)[names(maf) == "End_Position"] <- "End_position"
    maf[['Chromosome']] <- gsub("^chr", "", maf[['Chromosome']])

    snp_cols <- c('SNP', 'DNP', 'TNP', 'MNP')
    snp_maf <- maf[maf[['Variant_Type']] %in% snp_cols,]
    write.table(snp_maf,'${sample}.snp',sep = '\t',col.names = T,row.names = F,quote = F)

    indel_cols <- c('INS', 'DEL')
    indel_maf <- maf[maf[['Variant_Type']] %in% indel_cols,]
    write.table(indel_maf,'${sample}.indel',sep = '\t',col.names = T,row.names = F,quote = F)
    
    EOF
    """
}

process PREPROCESS_ICONICC {
    container "wchukwu/r-docker:latest"
    label 'process_low'

    input:
    tuple val(sample),
        path(segfile),
        path(processed_counts)
    path (capseg_rscript)

    output:
    path "${sample}.capseg.txt", emit: seg_data

    script:
    """
    Rscript '${capseg_rscript}' \
        -s '${segfile}' \
        -c '${processed_counts}' \
        -i '${sample}'
    """
}