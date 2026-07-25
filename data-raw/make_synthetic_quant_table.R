# ---------------------------------------------------------------------------
# Generate a SYNTHETIC MRMkit-style quant table.
#
# Purpose: document the expected input format for the scripts in analysis/ so
# that code can be read and reshaping checked without any real data.
#
# This file is repository tooling, not part of the thesis analysis. Every value
# it produces is drawn from a random number generator. It contains no real
# measurements and no real participants.
# ---------------------------------------------------------------------------

set.seed(20260726)

n_info_col   <- 3    # filename, type, batch
n_analyte    <- 40
n_metric_row <- 10   # row 7 = CoV, row 9 = D-ratio, row 10 = S/N
n_subject    <- 12
n_timepoint  <- 3
n_bqc        <- 4    # pooled QC samples

analytes <- paste0("analyte_", sprintf("%03d", seq_len(n_analyte)))

# --- metric rows -----------------------------------------------------------
# Deliberately spread across the gate boundaries so that filtering has an
# effect that is visible when the example is run.
metric <- matrix(NA_character_, nrow = n_metric_row, ncol = n_analyte,
                 dimnames = list(NULL, analytes))
for (r in seq_len(n_metric_row)) {
  metric[r, ] <- sprintf("%.2f", runif(n_analyte, 0, 100))
}
metric[7, ]  <- sprintf("%.2f", c(runif(n_analyte * 0.7, 2, 30),      # CoV pass
                                 runif(n_analyte - n_analyte * 0.7, 30, 90)))
metric[9, ]  <- sprintf("%.2f", c(runif(n_analyte * 0.6, 5, 50),      # D-ratio pass
                                 runif(n_analyte - n_analyte * 0.6, 50, 180)))
metric[10, ] <- sprintf("%.2f", c(runif(n_analyte * 0.85, 5, 400),    # S/N pass
                                  runif(n_analyte - n_analyte * 0.85, 0.5, 5)))

metric_info <- data.frame(
  filename = paste0("metric_row_", seq_len(n_metric_row)),
  type     = c("mean", "sd", "median", "mad", "min", "max",
               "CoV", "blank_ratio", "D_ratio", "SN"),
  batch    = NA_character_,
  stringsAsFactors = FALSE
)

# --- sample rows -----------------------------------------------------------
subjects   <- sprintf("SUBJ%03d", seq_len(n_subject))
timepoints <- c(0, 60, 120)

samples <- expand.grid(subject = subjects, time = timepoints,
                       stringsAsFactors = FALSE)
samples <- samples[order(samples$subject, samples$time), ]

sample_info <- data.frame(
  filename = paste0("SYNTH_", samples$subject, "_T", samples$time, ".mzML"),
  type     = "sample",
  batch    = rep(1:2, length.out = nrow(samples)),   # two alternating columns
  stringsAsFactors = FALSE
)

bqc_info <- data.frame(
  filename = paste0("SYNTH_BQC_", sprintf("%02d", seq_len(n_bqc)), ".mzML"),
  type     = "BQC",
  batch    = rep(1:2, length.out = n_bqc),
  stringsAsFactors = FALSE
)

n_row  <- nrow(sample_info) + nrow(bqc_info)
values <- matrix(sprintf("%.1f", rlnorm(n_row * n_analyte, meanlog = 8, sdlog = 1)),
                 nrow = n_row, ncol = n_analyte, dimnames = list(NULL, analytes))

# --- assemble --------------------------------------------------------------
out <- rbind(
  cbind(metric_info, metric),
  cbind(rbind(sample_info, bqc_info), values)
)

dir.create("data", showWarnings = FALSE)
write.table(out, file = "data/synthetic_quant_table.txt",
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

cat("wrote data/synthetic_quant_table.txt:",
    nrow(out), "rows x", ncol(out), "columns\n")
cat("  metric rows :", n_metric_row, "(CoV row 7, D-ratio row 9, S/N row 10)\n")
cat("  sample rows :", nrow(sample_info), "\n")
cat("  pooled QC   :", nrow(bqc_info), "\n")
cat("  analytes    :", n_analyte, "\n")
