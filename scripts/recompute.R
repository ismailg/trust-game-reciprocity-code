#!/usr/bin/env Rscript
# Run from the private working directory created by run_analysis.py.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source("scripts/reproduction_helpers.R")
release_mode <- Sys.getenv("RTG_RUN_MODE", "fresh")
release_B <- as.integer(Sys.getenv("RTG_BOOTSTRAP_DRAWS", "10000"))
if (!release_mode %in% c("fresh", "reuse-fits")) stop("Unknown computation mode")
if (is.na(release_B) || release_B < 8L || release_B %% 8L != 0L) {
  stop("Bootstrap count must be a positive multiple of eight")
}
Sys.setenv(VTC_HMM_FULL_RECOMPUTE = "1",
           VTC_TRUSTEE_BOOTSTRAP_DRAWS = release_B,
           VTC_TRUSTEE_SENSITIVITY_DRAWS = release_B)
if (release_B != 10000L) Sys.setenv(VTC_TRUSTEE_ALLOW_NONCANONICAL = "1")
release_chunks <- read_analysis_chunks()
release_deferred <- Filter(function(x) x$label == "trustee-long-run-sensitivity-results", release_chunks)[[1]]
release_started <- Sys.time()
dir.create("output", showWarnings = FALSE)
for (release_chunk in release_chunks) {
  release_label <- release_chunk$label
  if (release_label %in% c("release-state", "trustee-long-run-sensitivity-results")) next
  message("ANALYSIS STAGE: ", release_label)
  release_code <- release_chunk$code
  if (release_label == "investor-standardized-search-foundations") {
    # The validated investor search used treatment coding; the trustee search
    # used sum coding. Set this explicitly instead of inheriting session options.
    release_previous_contrasts <- getOption("contrasts")
    options(contrasts = c("contr.treatment", "contr.poly"))
    eval(parse(text = release_code), envir = .GlobalEnv)
    options(contrasts = release_previous_contrasts)
    next
  }
  if (release_mode == "reuse-fits" && release_label == "investor-standardized-search-run") {
    release_code <- strsplit(release_code, "if (investor_full_recompute) {", fixed = TRUE)[[1]][[1]]
    eval(parse(text = release_code), envir = .GlobalEnv)
    investor_standardized_search_result <- readRDS("private_fits/investor.rds")
    validate_saved_search_inputs(investor_standardized_search_result,
      list(simple = investor_simple_templates, group = investor_group_templates),
      investor_hmm_dat$investment, order_mod_truncdiscgaus)
    next
  }
  if (release_mode == "reuse-fits" && release_label == "trustee-standardized-model-search") {
    # Define the original functions, then authenticate saved input arrays before reuse.
    release_prefix <- strsplit(release_code, "if (trustee_full_recompute) {", fixed = TRUE)[[1]][[1]]
    eval(parse(text = release_prefix), envir = .GlobalEnv)
    release_search <- readRDS("private_fits/trustee.rds")
    validate_saved_search_inputs(release_search,
      list(simple = simple_HMMs, group = Ctrst_HMMs), trustee_dat$return_0, order_mod_vtdgaus)
    release_search$data_sha256 <- digest::digest(file = "Data/full_RTG_data.csv", algo = "sha256")
    dir.create("results/HMM/recomputed_from_analysis_main", recursive = TRUE, showWarnings = FALSE)
    saveRDS(release_search, "results/HMM/recomputed_from_analysis_main/standardized_search_result.rds")
    Sys.setenv(VTC_TRUSTEE_SEARCH_REUSE = "1")
  }
  if (release_label == "trustee-corrected-five-state-artifact") {
    release_code <- sub("trustee_corrected_result <- trustee_run_corrected_five_state_analysis(",
                        "trustee_corrected_result <- release_pool_primary(", release_code, fixed = TRUE)
  }
  if (release_label == "trustee-corrected-state-count-sensitivity") {
    release_code <- sub("      trustee_run_state_count_sensitivity(",
                        "      release_pool_state_count(", release_code, fixed = TRUE)
  }
  eval(parse(text = release_code), envir = .GlobalEnv)
  if (release_label == "trustee-corrected-five-state-recompute-functions") {
    # The separately fitted decoded-transition regressions used treatment
    # coding. Keep the manuscript's ordinary confidence-interval reporting.
    release_decoded_fit_original <- trustee_safe_decoded_glmm
    trustee_safe_decoded_glmm <- function(...) {
      args <- list(...)
      with_treatment_contrasts(function() do.call(release_decoded_fit_original, args))
    }
  }
  if (release_label == "trustee-corrected-five-state-artifact") {
    trustee_bootstrap_sampling_description <- paste(
      "We obtained 10,000 successful samples in eight independently seeded",
      "batches of 1,250 samples (seeds 11–18), pooling the samples before",
      "calculating confidence intervals and tests.")
  }
  if (release_label == "trustee-corrected-state-count-sensitivity") {
    trustee_sensitivity_sampling_description <- paste(
      "Analyses with the neighbouring four- and six-state trustee models",
      "assessed sensitivity to state count, using 10,000 successful bootstrap",
      "samples per model. For each model, eight independently seeded batches",
      "of 1,250 samples were pooled before calculating confidence intervals",
      "and applying the original comparison-family corrections.")
  }
}

message("ANALYSIS STAGE: demographic sensitivity analyses")
source("scripts/hmm_cache_utils.R")
source("scripts/trustee_model_helpers.R")
source("scripts/compute_covariates.R")
diagnoses <- c("CC", "BPD", "AD")
investment_grid <- c(5, 10, 15)
release_previous_contrasts <- getOption("contrasts")
options(contrasts = c("contr.treatment", "contr.poly"))
release_covariates <- compute_covariates(trustee_selected_model,
  prepare_rtg_data(".")$trustee_dat %>% arrange(subject_ID, round),
  "results/HMM/corrected_submission/covariate_sensitivity_2026-07-29_v2")
options(contrasts = release_previous_contrasts)
source("scripts/demographic_functions.R")
source("scripts/prepare_demographics.R")
source("scripts/sample_demographics.R")
trustee_long_run_sens_result <- compute_demographic_cooperation(trustee_selected_model, release_B)
trustee_long_run_sens_dir <- "results/HMM/long_run_covariate_2026-09-30"
dir.create(trustee_long_run_sens_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(trustee_long_run_sens_result$counts,
                 file.path(trustee_long_run_sens_dir, "sample_counts.csv"))
# Populate reporting objects from the newly computed result, rather than
# attempting to authenticate a historical archive that is not an input here.
release_deferred_code <- substring(release_deferred$code,
  regexpr("trustee_long_run_sens_order <-", release_deferred$code, fixed = TRUE)[[1]])
if (release_B == 10000L) {
  eval(parse(text = release_deferred_code), envir = .GlobalEnv)
}
saveRDS(list(primary = trustee_corrected_result,
             state_counts = trustee_state_count_sensitivity,
             demographics = trustee_long_run_sens_result,
             covariates = release_covariates), "output/numerical_results.rds")
release_code_paths <- sort(c(list.files(pattern = "\\.(R|Rmd|py)$"),
  list.files("scripts", pattern = "\\.R$", full.names = TRUE)))
release_provenance <- list(mode = release_mode, bootstrap_draws = release_B,
  code_sha256 = setNames(vapply(release_code_paths,
    function(p) digest::digest(file = p, algo = "sha256"), character(1)), release_code_paths),
  inputs = setNames(vapply(c("Data/full_RTG_data.csv", "Data/demographics.csv"),
    function(p) digest::digest(file = p, algo = "sha256"), character(1)),
    c("Data/full_RTG_data.csv", "Data/demographics.csv")),
  R_version = R.version.string, RNG_kind = RNGkind(),
  finished_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(Sys.time(), release_started, units = "secs")))
jsonlite::write_json(release_provenance, "output/run_provenance.json", pretty = TRUE, auto_unbox = TRUE)
save(list = ls(envir = .GlobalEnv), file = "output/analysis_state.RData", envir = .GlobalEnv)
message("RECOMPUTATION COMPLETE; all generated files contain restricted analysis material.")
