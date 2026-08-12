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

  # tau = 2 * (copy ratio). seg.mean may be in log2 scale (diploid ~ 0, deletions negative) or
  # in linear copy-ratio scale (diploid ~ 1, strictly positive). Older ICONICC output was linear;
  # newer output is log2. Detect the scale and compute tau accordingly (a linear value fed to
  # 2*2^seg.mean would blow up amplicons to astronomical tau).
  seg_mean_vals <- cancer_allele_segments$seg.mean
  seg_mean_is_log2 <- (min(seg_mean_vals, na.rm = TRUE) < -0.1) ||
    (median(seg_mean_vals, na.rm = TRUE) < 0.5)
  message("seg.mean scale detected: ", if (seg_mean_is_log2) "log2" else "linear")
  if (seg_mean_is_log2) {
    cancer_allele_segments[, tau := 2 * (2^seg.mean)]
  } else {
    cancer_allele_segments[, tau := 2 * seg.mean]
  }

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

  # Segment-level f estimation: noise-deconvolved minor-allele deviation + depth-aware balance snap.
  #
  # Per-bin AFMIN = min(AF1, AF2) is biased BELOW 0.5 even for genuinely balanced segments, because
  # taking the per-bin minimum of two noisy ~0.5 fractions systematically selects the downward
  # fluctuation. A mean/mode of AFMIN therefore never collapses balanced segments to 0.5.
  #
  # AF1 is the MEAN over `snp_count` het-SNP fractions (ICONICC AC_calc), so its sampling variance
  # is sigma_i^2 = v / snp_count_i, where v is a per-sample per-SNP allelic overdispersion. For a
  # bin, dev = 0.5 - AFMIN and E[dev^2] = (0.5 - f)^2 + sigma_i^2. Deconvolve:
  #   (0.5 - f) = sqrt(max(0, mean(dev^2) - mean(sigma_i^2)))
  # v is SELF-CALIBRATED from the data: for balanced bins snp_count*dev^2 ~ v, and real imbalance
  # only inflates it, so a low quantile of snp_count*dev^2 isolates the balanced-noise floor.
  # Balanced segment -> excess ~ 0 -> f = 0.5; real LOH -> recovers the true minor fraction;
  # depth-aware via snp_count. A consistency test then snaps to f = 0.5 any
  # segment whose per-bin excess deviation is not significantly positive across bins (see T_SNAP).
  # Because AF2 = 1 - AF1 exactly, dev = 0.5 - AFMIN = |AF1 - 0.5|, so for a balanced bin
  # snp_count * dev^2 ~ v * chi-square(1) (right-skewed). The q-quantile of that product equals
  # v * qchisq(q, 1); dividing the observed quantile by qchisq(q, 1) recovers the scale v. Using the
  # raw quantile (without this correction) severely underestimates v and under-collapses balanced
  # segments. A quantile at/below the median keeps the estimate in the balanced-dominated region.
  Q_FLOOR <- 0.5     # tunable: quantile used for the balanced-noise calibration (<=0.5 recommended)
  T_SNAP  <- as.numeric(Sys.getenv("T_SNAP", "3"))
                     # tunable: consistency-test threshold. Snap f -> 0.5 when the per-bin excess of
                     # squared minor-deviation over the balanced-noise floor is NOT significantly
                     # positive: t = mean(excess) / (sd(excess)/sqrt(n_bins)) < T_SNAP. Larger =>
                     # snaps more segments to balanced. Replaces the old Z_SNAP/TOL magnitude floor,
                     # which collapsed mild but real (low-purity) imbalance because it tested the
                     # deviation's *size* rather than its *consistency* across bins.
  hs <- processed_data[!is.na(AFMIN) & !is.na(snp_count) & snp_count > 0]
  v_hat <- as.numeric(quantile(hs$snp_count * (0.5 - hs$AFMIN)^2, Q_FLOOR, na.rm = TRUE)) /
    qchisq(Q_FLOOR, df = 1)
  message(sprintf("self-calibrated per-SNP allelic variance v = %.5f (floor q = %.2f)", v_hat, Q_FLOOR))

  minor_af_stats <- processed_data[, {
    ok <- !is.na(AFMIN) & !is.na(snp_count) & snp_count > 0
    afmin_values <- AFMIN[ok]
    sc <- as.numeric(snp_count[ok])
    n_bins <- length(afmin_values)
    bins_with_sufficient_snps <- sum(snp_count > 5, na.rm = TRUE)

    afmin_mean <- if (n_bins > 0) mean(afmin_values) else NA_real_
    afmin_sd   <- if (n_bins > 1) sd(afmin_values) else NA_real_

    if (n_bins < 2) {
      f_value <- afmin_mean
    } else {
      dev2  <- (0.5 - afmin_values)^2       # per-bin squared deviation
      sig2  <- v_hat / sc                   # per-bin AF variance (self-calibrated, depth-aware)
      delta <- sqrt(max(0, mean(dev2) - mean(sig2)))   # noise-deconvolved minor deviation
      excess <- dev2 - sig2                            # per-bin excess over the balanced-noise floor
      # Snap on the CONSISTENCY of the excess across bins, not its magnitude. Under a genuinely
      # balanced segment E[excess] = 0 with bin-to-bin scatter sd(excess); a real imbalance (even a
      # mild, low-purity one) shifts every bin's excess positive, so mean(excess) sits many empirical
      # SEs above 0. delta and the t-statistic are coupled (t ~ delta^2 * sqrt(n_bins) / sd(excess)),
      # so a tiny deviation cannot score a large t: the statistic folds effect-size and consistency
      # into one test, keeping a small-but-consistent shift (real LOH) while snapping a small-but-
      # noisy one (min/max artifact) -- unlike the old absolute-magnitude TOL floor.
      se_excess <- if (n_bins > 1) sd(excess) / sqrt(n_bins) else Inf
      t_stat    <- if (se_excess > 0) mean(excess) / se_excess else 0
      f_value   <- if (t_stat < T_SNAP) 0.5 else (0.5 - delta)
    }

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

  # NOTE (2026-08): dropped the legacy SegLabelCNLOH column. ABSOLUTE never reads it (grep-confirmed across
  # getzlab/ABSOLUTE v1.5; the allelic reader read.delim->AllelicMakeSegObj indexes seg cols BY NAME, with no
  # presence-check and no column-count/positional assumption), and the canonical generators (getzlab/ASCAT-parser,
  # aaronmck/CapSeg AllelicCapseg) both OMIT it. Our old 0/2 rule was also malformed vs the canonical 0/1
  # (false/true CNLOH) spec and thresholded on the switch-biased raw AFMIN. Output now matches ascat-parser's acs.
  cancer_allele_segments[, length := loc.end - loc.start]

  capseg_output <- cancer_allele_segments[, c(
    "chrom", "loc.start", "loc.end", "num.mark", "length", "n_hets",
    "f", "tau", "copy_number_sd", "minor_copy_estimate", "minor_copy_sd",
    "major_copy_estimate", "major_copy_sd"
  )]

  colnames(capseg_output) <- c(
    "Chromosome", "Start.bp", "End.bp", "n_probes", "length", "n_hets",
    "f", "tau", "sigma.tau", "mu.minor", "sigma.minor",
    "mu.major", "sigma.major"
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