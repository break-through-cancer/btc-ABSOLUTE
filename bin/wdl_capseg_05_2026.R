#################################################################
##                      READ IN LIBRARIES                      ##
#################################################################
library(data.table)
library(parallel)
library(optparse)
library(tidyverse)
library(dplyr)

##################################################################
##                DEFINE INPUT OPTIONS AND FLAGS                ##
##################################################################
option_list <- list(
  make_option(c("-s", "--segfile"), type = "character", default = NA,
              help = "Path to iconicc segfile output"),
  make_option(c("-c", "--processed_cts"), type = "character", default = NA,
              help = "Path to processed allelic counts file"),
  make_option(c("-i", "--participant_id"), type = "character", default = NA,
              help = "Unique sample identiifer")
  
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

iconicc_segments <- fread(opt$segfile)
processed_counts <- fread(opt$processed_cts)
sample_id <- opt$participant_id

canonical_contigs <- c(paste0("chr", 1:22), "chrX", "chrY")

#' Assign segment ID to each bin based on chromosomal overlap.
#'
#' Intersects each bin in the processed counts table with the segment table and returns the
#' overlapping segment ID. Used to map bins to their corresponding segments for
#' calculating allele-specific copy number estimates.
#'
#' @param row_index Integer index of the bin row to process.
#' @param source_dt Data table of processed counts with CONTIG, START, END columns.
#' @param segment_table Data table of segments with chrom, loc.start, loc.end, SegmentID columns.
#' @param segment_source Character prefix for the output column name.
#'
#' @return Updated row with segment ID assignment.
assign_segment_id <- function(row_index, source_dt, segment_table, segment_source) {
  current_row <- source_dt[row_index, ]
  current_chrom <- current_row$CONTIG
  chrom_segment_info <- subset(segment_table, chrom == current_chrom)
  segment_id <- chrom_segment_info[
    which((as.numeric(current_row$START) >= loc.start) &
            (as.numeric(current_row$END) <= (loc.end + 24999))),
  ]$SegmentID
  if (identical(integer(0), segment_id)) {
    segment_id <- NA
  }
  return(current_row[, paste0(segment_source, "_segment") := segment_id])
}


#' Convert ICONICC segmentation output to cancer allele-specific copy number (CAPSEG) format.
#'
#' Takes an ICONICC segment file and processed allelic counts, calculates cancer allele-specific
#' copy number estimates (tau), allele frequencies (f), and associated standard errors. Outputs
#' a CAPSEG file compatible with downstream CNV calling and visualization tools.
#'
#' @param segment_table Data table of final ICONICC segments.
#' @param processed_data Data table of processed allelic counts.
#'
#' @return Data table in CAPSEG format with one row per segment.
convert_iconicc_to_capseg <- function(segment_table, processed_data) {
  estimate_minor_fraction_from_afmin <- function(afmin_values, trial_counts) {
    # Peak-finding:
    # 1) Keep only bins with valid AFMIN and usable SNP support.
    # 2) Fit a weighted kernel density on AFMIN in [0, 0.5].
    # 3) Return the AFMIN mode (x at max density) as the segment-level minor fraction.
    valid_idx <- which(!is.na(afmin_values) & !is.na(trial_counts) & trial_counts > 0)

    if (length(valid_idx) < 2) {
      return(NA_real_)
    }

    # Weight by SNP support so bins with more evidence contribute more to the peak.
    afmin_values <- afmin_values[valid_idx]
    trial_counts <- as.numeric(trial_counts[valid_idx])
    density_weights <- trial_counts / sum(trial_counts)

    density_fit <- stats::density(
      afmin_values,
      weights = density_weights,
      from = 0,
      to = 0.5,
      n = 1024,
      cut = 0
    )

    density_fit$x[which.max(density_fit$y)]
  }

  segment_table[, SegmentID := seq_len(nrow(segment_table))]
  cancer_allele_segments <- segment_table

  processed_data <- processed_data[order(factor(CONTIG, levels = canonical_contigs), as.numeric(START))]
  processed_data <- rbindlist(mclapply(
    seq_len(nrow(processed_data)),
    assign_segment_id,
    processed_data,
    segment_table,
    segment_source = "final",
    mc.cores = 7
  ))

  processed_data[, copy_number_linear := (2^log2_tangent)]
  cancer_allele_segments[, tau := 2 * (2^seg.mean)]

  copy_number_sd_by_segment <- aggregate(
    processed_data$copy_number_linear,
    by = list(processed_data$final_segment),
    FUN = function(x) {
      return(sd(x, na.rm = TRUE))
    }
  )
  copy_number_sd_by_segment <- as.data.table(copy_number_sd_by_segment)
  colnames(copy_number_sd_by_segment) <- c("SegmentID", "copy_number_sd")
  cancer_allele_segments <- merge(cancer_allele_segments, copy_number_sd_by_segment, by = "SegmentID")

  processed_data[, AFMIN := ifelse(snp_count < 5, NA, AFMIN)]
  processed_data[, AFMAX := ifelse(snp_count < 5, NA, AFMAX)]
  processed_data[, AF1 := ifelse(snp_count < 5, NA, AF1)]
  processed_data[, AF2 := ifelse(snp_count < 5, NA, AF2)]

  # Segment-level f estimation via AFMIN peak-finding
  # For each segment, collect AFMIN values from high-confidence bins.
  # Use the AFMIN mean as a fallback for sparse segments.
  # If enough bins are available, replace fallback with weighted KDE mode.
  # Snap near-balanced values to exactly 0.5 to reduce numerical jitter.
  min_bins_for_density <- 8
  # threshold for rounding to balanced f
  # balance_snap_tolerance <- 0.00
  minor_af_stats <- processed_data[, {
    # Per-segment AFMIN subset after earlier SNP-count filtering.
    afmin_values <- AFMIN[!is.na(AFMIN)]
    bins_with_sufficient_snps <- sum(snp_count > 5, na.rm = TRUE)

    afmin_mean <- if (length(afmin_values) > 0) mean(afmin_values, na.rm = TRUE) else NA_real_
    afmin_sd <- if (length(afmin_values) > 1) sd(afmin_values, na.rm = TRUE) else NA_real_
    # Fallback for short/sparse segments where density mode is unstable.
    fitted_minor_fraction <- afmin_mean

    if (length(afmin_values) >= min_bins_for_density) {
      # Rebuild aligned weights from SNP counts for the non-missing AFMIN bins.
      trial_counts <- as.integer(round(snp_count[!is.na(AFMIN)]))
      valid_idx <- which(!is.na(trial_counts) & trial_counts > 0L)
      trial_counts <- trial_counts[valid_idx]
      afmin_test_values <- afmin_values[valid_idx]

      # Main estimator: AFMIN peak from weighted KDE.
      fitted_minor_fraction <- estimate_minor_fraction_from_afmin(afmin_test_values, trial_counts)
    }

    f_value <- fitted_minor_fraction

    # Snap near-balanced segments to exactly 0.5 to reduce numerical jitter.
    # if (!is.na(f_value) && abs(f_value - 0.5) <= balance_snap_tolerance) {
    #   f_value <- 0.5
    # }

    list(
      allele_freq_min_mean = afmin_mean,
      allele_freq_min_sd = afmin_sd,
      bins_snp_count_gt5 = bins_with_sufficient_snps,
      f = f_value
    )
  }, by = final_segment]
  setnames(minor_af_stats, "final_segment", "SegmentID")

  major_af_stats <- aggregate(
    processed_data$AFMAX,
    by = list(processed_data$final_segment),
    FUN = function(x) {
      return(sd(x, na.rm = TRUE))
    }
  )
  major_af_stats <- as.data.table(major_af_stats)
  colnames(major_af_stats) <- c("SegmentID", "allele_freq_maj_sd")

  cancer_allele_segments <- Reduce(merge, list(cancer_allele_segments, minor_af_stats, major_af_stats))

  cancer_allele_segments[, minor_copy_estimate := f * tau]
  cancer_allele_segments[, major_copy_estimate := (1 - f) * tau]
  cancer_allele_segments[, minor_copy_sd := sqrt((copy_number_sd^2) + (allele_freq_min_sd^2))]
  cancer_allele_segments[, major_copy_sd := sqrt((copy_number_sd^2) + (allele_freq_maj_sd^2))]

  snp_count_by_segment <- aggregate(
    processed_data[snp_count > 5, snp_count],
    by = list(processed_data[snp_count > 5, final_segment]),
    FUN = sum
  )
  snp_count_by_segment <- as.data.table(snp_count_by_segment)
  colnames(snp_count_by_segment) <- c("SegmentID", "n_hets")
  cancer_allele_segments <- merge(cancer_allele_segments, snp_count_by_segment, by = "SegmentID")

  cancer_allele_segments[, loh_label := ifelse(
    (tau >= 1.8 & tau <= 2.2) & allele_freq_min_mean <= 0.05,
    0,
    2
  )]

  cancer_allele_segments[, length := loc.end - loc.start]

  capseg_output <- cancer_allele_segments[, c(
    "chrom", "loc.start", "loc.end", "num.mark", "length", "n_hets",
    "f", "tau", "copy_number_sd", "minor_copy_estimate", "minor_copy_sd",
    "major_copy_estimate", "major_copy_sd", "loh_label"
  )]

  colnames(capseg_output) <- c(
    "Chromosome", "Start.bp", "End.bp", "n_probes", "length", "n_hets",
    "f", "tau", "sigma.tau", "mu.minor", "sigma.minor",
    "mu.major", "sigma.major", "SegLabelCNLOH"
  )

  capseg_output[, Chromosome := gsub("chr", "", Chromosome)]
  capseg_output <- capseg_output %>% drop_na(tau)

  write.table(
    capseg_output,
    file = paste0(sample_id, ".capseg.txt"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )

  return(capseg_output)
}

# execute main conversion pipeline.
convert_iconicc_to_capseg(iconicc_segments, processed_counts)