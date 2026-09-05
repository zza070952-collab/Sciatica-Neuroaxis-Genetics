#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: phase8_drg_recurrence.R DRG_EXPRESSION_TSV OUTPUT_DIR")
workbook <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv_stable <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

genes <- c("TXNL1", "MAPK3", "FGFR3", "PDPR", "GFPT1")
cell_types <- c("Nociceptor", "Mechanoreceptor", "Schwann cell")
primary_fraction <- 0.05
primary_min_cells <- 20L

x <- read.delim(workbook, check.names = FALSE) %>%
  filter(gene %in% genes, cell_type %in% cell_types) %>%
  mutate(
    unique_donor_key = biological_replicate_key,
    n_expressing_approx = round(n_cells * fraction_expressing),
    assessable_primary = n_cells >= primary_min_cells,
    positive_primary = assessable_primary & fraction_expressing >= primary_fraction
  )

if (n_distinct(x$sample_preparation) != 6L || n_distinct(x$unique_donor_key) != 5L) {
  stop("Expected 6 preparations and 5 unique donors in GSE168243")
}

donor_status <- x %>%
  group_by(gene, cell_type, unique_donor_key) %>%
  summarise(
    donor_assessable = any(assessable_primary),
    donor_positive = any(positive_primary),
    preparations = paste(sort(unique(sample_preparation)), collapse = ","),
    .groups = "drop"
  )

summary <- x %>%
  group_by(gene, cell_type) %>%
  summarise(
    dataset_preparations = 6L,
    dataset_unique_donors = 5L,
    observed_preparations = n_distinct(sample_preparation),
    assessable_preparations = sum(assessable_primary),
    expression_positive_preparations = sum(positive_primary),
    .groups = "drop"
  ) %>%
  left_join(
    donor_status %>%
      group_by(gene, cell_type) %>%
      summarise(
        observed_unique_donors = n_distinct(unique_donor_key),
        assessable_unique_donors = sum(donor_assessable),
        expression_positive_unique_donors = sum(donor_positive),
        .groups = "drop"
      ),
    by = c("gene", "cell_type")
  ) %>%
  mutate(
    detection_fraction_threshold = primary_fraction,
    minimum_nuclei_per_preparation_cell_class = primary_min_cells,
    preparation_recurrence = paste0(expression_positive_preparations, "/", assessable_preparations, " assessable preparations"),
    donor_recurrence = paste0(expression_positive_unique_donors, "/", assessable_unique_donors, " assessable unique donors"),
    interpretation = case_when(
      assessable_preparations == 0 ~ "UNAVAILABLE_IN_ASSESSABLE_PREPARATIONS",
      expression_positive_preparations == assessable_preparations ~ "RECURRENT_ACROSS_ASSESSABLE_PREPARATIONS",
      expression_positive_preparations == 0 ~ "NOT_DETECTED_AT_LOCKED_THRESHOLD",
      TRUE ~ "PARTIAL_OR_SPARSE_RECURRENCE"
    )
  ) %>%
  arrange(factor(gene, levels = genes), factor(cell_type, levels = cell_types))

before_after <- summary %>%
  transmute(
    gene, cell_type,
    before_dataset_coverage_preparations = observed_preparations,
    before_dataset_coverage_unique_donors = observed_unique_donors,
    after_positive_preparations = expression_positive_preparations,
    after_assessable_preparations = assessable_preparations,
    after_positive_unique_donors = expression_positive_unique_donors,
    after_assessable_unique_donors = assessable_unique_donors,
    correction = "Dataset coverage is no longer used as expression-positive recurrence"
  )

sensitivity_grid <- tidyr::crossing(
  detection_fraction_threshold = c(0, 0.01, 0.05, 0.10),
  minimum_nuclei = c(1L, 10L, 20L)
)
sensitivity <- bind_rows(lapply(seq_len(nrow(sensitivity_grid)), function(i) {
  thr <- sensitivity_grid$detection_fraction_threshold[[i]]
  min_n <- sensitivity_grid$minimum_nuclei[[i]]
  y <- x %>%
    mutate(
      assessable = n_cells >= min_n,
      positive = assessable & if (thr == 0) fraction_expressing > 0 else fraction_expressing >= thr
    )
  ds <- y %>%
    group_by(gene, cell_type, unique_donor_key) %>%
    summarise(donor_assessable = any(assessable), donor_positive = any(positive), .groups = "drop")
  y %>%
    group_by(gene, cell_type) %>%
    summarise(
      assessable_preparations = sum(assessable),
      positive_preparations = sum(positive),
      .groups = "drop"
    ) %>%
    left_join(
      ds %>% group_by(gene, cell_type) %>% summarise(
        assessable_unique_donors = sum(donor_assessable),
        positive_unique_donors = sum(donor_positive),
        .groups = "drop"
      ),
      by = c("gene", "cell_type")
    ) %>%
    mutate(detection_fraction_threshold = thr, minimum_nuclei = min_n)
}))

detail <- x %>%
  select(
    gene, cell_type, donor, sample_preparation, unique_donor_key, n_cells,
    mean_log_expression, fraction_expressing, n_expressing_approx,
    assessable_primary, positive_primary
  ) %>%
  arrange(factor(gene, levels = genes), factor(cell_type, levels = cell_types), sample_preparation)

write_tsv_stable(summary, file.path(out_dir, "DRG_expression_positive_recurrence.tsv"))
write_tsv_stable(before_after, file.path(out_dir, "DRG_recurrence_before_after.tsv"))
write_tsv_stable(sensitivity, file.path(out_dir, "DRG_recurrence_threshold_sensitivity.tsv"))
write_tsv_stable(detail, file.path(out_dir, "DRG_preparation_level_detection.tsv"))

definition <- c(
  "# DRG expression-positive recurrence definition",
  "",
  "The Phase 3 files did not state a numerical expression-positive recurrence threshold. They reported fraction_expressing but counted every available donor row as donor coverage, even when fraction_expressing was zero or near zero. Phase 8 therefore separates dataset coverage from expression-positive recurrence.",
  "",
  "## Locked primary rule for this correction",
  "",
  paste0("- A preparation-by-cell-class record was assessable when it contained at least ", primary_min_cells, " nuclei."),
  paste0("- An assessable preparation was expression-positive when fraction_expressing was at least ", primary_fraction, "."),
  "- hDRG3 and hDRG5 were counted as separate preparations but as one unique donor.",
  "- A unique donor was positive when at least one assessable preparation from that donor was positive; both preparation-level and donor-level recurrence are reported.",
  "- Records below the minimum cell count were unavailable for recurrence inference rather than coded negative.",
  "",
  "## Boundary",
  "",
  "This uniformly applied descriptive gate was locked for the Phase 8 correction and was not used to re-rank candidates. Threshold sensitivity at 0%, 1%, 5% and 10% fraction expressing and minimum-cell thresholds of 1, 10 and 20 is reported separately. Normal adult DRG reference expression does not constitute case-control differential expression."
)
writeLines(definition, file.path(out_dir, "DRG_threshold_definition.md"), useBytes = TRUE)
