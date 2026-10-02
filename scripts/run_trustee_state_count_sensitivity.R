#!/usr/bin/env Rscript

# Six- and seven-state sensitivity check for the corrected trustee analysis.
# States with mean returns below 40% form the lower-return block; states at or
# above 40% form the cooperative block. This avoids forcing state labels from
# the five-state solution onto models with a different number of states.

Sys.setenv(TRUSTEE_ROBUSTNESS_LIBRARY_ONLY = "1")
source(file.path("scripts", "run_trustee_state_count_robustness.R"))
Sys.unsetenv("TRUSTEE_ROBUSTNESS_LIBRARY_ONLY")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[length(hit)]])
}

state_counts <- as.integer(strsplit(get_arg("--states", "6,7"), ",", fixed = TRUE)[[1]])
bootstrap_draws <- as.integer(get_arg("--bootstrap", "1000"))
seed <- as.integer(get_arg("--seed", "11"))
selection_artifact_rds <- get_arg(
  "--selection-artifact",
  file.path(
    "results", "HMM", "corrected_submission",
    "model_selection_final_2026-07-21_v2", "validated_five_state_model.rds"
  )
)
output_dir <- get_arg(
  "--output-dir",
  file.path("results", "HMM", "corrected_submission", "state_count_sensitivity")
)
if (any(!state_counts %in% 6:7)) stop("--states must contain 6 and/or 7")
if (!is.finite(bootstrap_draws) || bootstrap_draws < 0) stop("--bootstrap must be >= 0")
if (!file.exists(selection_artifact_rds)) {
  stop("Canonical selection artifact is missing: ", selection_artifact_rds)
}
if (dir.exists(output_dir) && length(list.files(
  output_dir, all.files = TRUE, no.. = TRUE
)) > 0L) {
  stop("Output directory must be new or empty: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")
atomic_write_csv <- function(data, path) {
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  readr::write_csv(data, temporary, na = "")
  if (!file.rename(temporary, path)) stop("Atomic rename failed for ", path)
  invisible(path)
}
atomic_save_rds <- function(object, path) {
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Atomic rename failed for ", path)
  invisible(path)
}

selection_artifact <- readRDS(selection_artifact_rds)
required_selection_fields <- c(
  "selected_state_count", "selected_fit_label", "gate", "model_comparison",
  "k6_search_dir", "k7_search_dir", "provenance"
)
if (!is.list(selection_artifact) ||
    !all(required_selection_fields %in% names(selection_artifact)) ||
    !identical(as.integer(selection_artifact$selected_state_count), 5L) ||
    !identical(as.character(selection_artifact$selected_fit_label), "perturbed_warm_02") ||
    !is.data.frame(selection_artifact$gate) || nrow(selection_artifact$gate) == 0L ||
    anyDuplicated(selection_artifact$gate$check) || !all(selection_artifact$gate$passed)) {
  stop("Selection artifact is not the approved fully gated result")
}

model_path_for_state <- function(state_count) {
  search_dir <- selection_artifact[[paste0("k", state_count, "_search_dir")]]
  model_path <- file.path(search_dir, "fits", "best_ordered_fit.rds")
  if (!file.exists(model_path)) stop("Validated sensitivity model is missing: ", model_path)
  expected_hash <- selection_artifact$provenance[[paste0(
    "k", state_count, "_best_fit_sha256"
  )]]
  if (is.null(expected_hash) || !identical(sha256_file(model_path), expected_hash)) {
    stop(state_count, "-state sensitivity-model hash does not match selection artifact")
  }
  model_path
}

rename_broad_metrics <- function(data) {
  data %>% mutate(metric = recode(
    metric,
    entry = "entry_from_lower_return",
    exploitative_exit = "exit_to_lower_return"
  ))
}

summarise_draws_finite <- function(draws, point) {
  levels <- draws %>%
    group_by(diagnosis, investment, metric) %>%
    summarise(
      bootstrap_mean = mean(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    left_join(point %>% rename(point_estimate = value), by = c("diagnosis", "investment", "metric"))

  contrast_draws <- draws %>%
    dplyr::select(iteration, diagnosis, investment, metric, value) %>%
    pivot_wider(names_from = diagnosis, values_from = value) %>%
    mutate(
      `BPD-CC` = BPD - CC,
      `AD-CC` = AD - CC,
      `AD-BPD` = AD - BPD
    ) %>%
    dplyr::select(iteration, investment, metric, `BPD-CC`, `AD-CC`, `AD-BPD`) %>%
    pivot_longer(
      cols = c(`BPD-CC`, `AD-CC`, `AD-BPD`),
      names_to = "contrast",
      values_to = "difference"
    )

  point_contrasts <- point %>%
    pivot_wider(names_from = diagnosis, values_from = value) %>%
    mutate(
      `BPD-CC` = BPD - CC,
      `AD-CC` = AD - CC,
      `AD-BPD` = AD - BPD
    ) %>%
    dplyr::select(investment, metric, `BPD-CC`, `AD-CC`, `AD-BPD`) %>%
    pivot_longer(
      cols = c(`BPD-CC`, `AD-CC`, `AD-BPD`),
      names_to = "contrast",
      values_to = "point_difference"
    )

  contrasts <- contrast_draws %>%
    group_by(investment, metric, contrast) %>%
    summarise(
      bootstrap_mean = mean(difference),
      lower = quantile(difference, 0.025),
      upper = quantile(difference, 0.975),
      p_two_sided = min(
        1,
        2 * (min(sum(difference <= 0), sum(difference >= 0)) + 1) /
          (dplyr::n() + 1)
      ),
      .groups = "drop"
    ) %>%
    mutate(
      correction_family = if_else(
        metric %in% c("retention", "entry_from_lower_return"),
        "retention_and_entry",
        metric
      )
    ) %>%
    group_by(correction_family) %>%
    mutate(p_holm_family = p.adjust(p_two_sided, method = "holm")) %>%
    ungroup() %>%
    left_join(point_contrasts, by = c("investment", "metric", "contrast"))

  list(levels = levels, contrasts = contrasts, draws = draws)
}

cooperative_occupancy_analysis <- function(model, cooperative, state_count) {
  smoothing <- as.matrix(posterior(model, type = "smoothing"))
  if (nrow(smoothing) != nrow(trustee_data_base) || ncol(smoothing) != state_count) {
    stop(state_count, "-state smoothed-posterior dimensions are invalid")
  }
  participant <- trustee_data_base %>%
    transmute(
      subject_ID,
      group,
      cooperative_probability = rowSums(smoothing[, cooperative, drop = FALSE])
    ) %>%
    group_by(subject_ID, group) %>%
    summarise(cooperative_occupancy = mean(cooperative_probability), .groups = "drop")
  summary <- participant %>%
    group_by(group) %>%
    summarise(
      participants = n(),
      mean_occupancy = mean(cooperative_occupancy),
      standard_error = sd(cooperative_occupancy) / sqrt(n()),
      .groups = "drop"
    )
  fitted <- lm(cooperative_occupancy ~ group, data = participant)
  pair_table <- as.data.frame(summary(
    contrast(emmeans(fitted, ~ group), "pairwise", adjust = "holm"),
    infer = c(TRUE, TRUE)
  ))
  pairs <- as_tibble(pair_table) %>% transmute(
    contrast,
    difference = estimate,
    lower = lower.CL,
    upper = upper.CL,
    p_holm = p.value
  )
  list(participant = participant, summary = summary, pairs = pairs)
}

all_results <- list()
for (state_count in state_counts) {
  model_path <- model_path_for_state(state_count)
  model <- order_mod_vtdgaus(readRDS(model_path))
  comparison_row <- selection_artifact$model_comparison %>%
    filter(.data$state_count == .env$state_count)
  convergence_message <- paste(model@message, collapse = " ")
  if (nrow(comparison_row) != 1L || nstates(model) != state_count ||
      !convergence_message %in% c(
        "Log likelihood converged to within tol. (relative change)",
        "Log likelihood converged to within tol. (absolute change)"
      ) ||
      !isTRUE(all.equal(
        as.numeric(logLik(model)),
        as.numeric(comparison_row$log_likelihood),
        tolerance = 1e-8
      )) ||
      !identical(
        as.integer(attr(logLik(model), "df")),
        as.integer(comparison_row$parameters)
      )) {
    stop(state_count, "-state model does not match the final selection table")
  }
  means <- expected_state_means(model)
  lower_return <- which(means < 0.40)
  cooperative <- which(means >= 0.40)
  if (length(lower_return) == 0L || length(cooperative) == 0L) {
    stop("The 40% threshold does not divide the ", state_count, "-state model")
  }

  state_dir <- file.path(output_dir, paste0("states", state_count))
  dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
  dictionary <- tibble(
    state_count,
    state = seq_len(state_count),
    expected_return_fraction = means,
    broad_category = if_else(state %in% cooperative, "cooperative", "lower_return")
  )
  atomic_write_csv(dictionary, file.path(state_dir, "state_dictionary.csv"))

  occupancy <- cooperative_occupancy_analysis(model, cooperative, state_count)
  atomic_write_csv(
    occupancy$participant,
    file.path(state_dir, "probability_weighted_cooperative_occupancy.csv")
  )
  atomic_write_csv(
    occupancy$summary,
    file.path(state_dir, "probability_weighted_cooperative_occupancy_summary.csv")
  )
  atomic_write_csv(
    occupancy$pairs,
    file.path(state_dir, "probability_weighted_cooperative_occupancy_pairwise.csv")
  )

  point <- direct_transition_metric_grid(model, cooperative, lower_return) %>%
    rename_broad_metrics() %>%
    mutate(state_count = state_count)
  atomic_write_csv(point, file.path(state_dir, "transition_metric_point_estimates.csv"))

  message("Fitting corrected transition summaries for uncertainty at ", state_count, " states")
  xi <- expected_transitions_xi(model, trustee_data_base)
  frames <- build_multinom_frames(model, trustee_data_base, xi)
  blocks <- lapply(frames, fit_block_multinom)
  if (any(vapply(blocks, is.null, logical(1)))) stop("A transition block failed")
  transition_summary_point <- transition_metric_grid(blocks, cooperative, lower_return) %>%
    rename_broad_metrics() %>%
    mutate(state_count = state_count)
  point_agreement <- point %>%
    rename(direct_fitted_model = value) %>%
    left_join(
      transition_summary_point %>% rename(transition_summary_model = value),
      by = c("state_count", "diagnosis", "investment", "metric")
    ) %>%
    mutate(absolute_difference = abs(direct_fitted_model - transition_summary_model))
  if (nrow(point_agreement) != nrow(point) ||
      any(!is.finite(point_agreement$absolute_difference)) ||
      max(point_agreement$absolute_difference) > 0.005) {
    stop("Transition-summary model does not adequately reproduce the fitted HMM")
  }
  atomic_write_csv(
    transition_summary_point,
    file.path(state_dir, "transition_summary_model_point_estimates.csv")
  )
  atomic_write_csv(
    point_agreement,
    file.path(state_dir, "direct_vs_transition_summary_agreement.csv")
  )

  bootstrap <- NULL
  bootstrap_status <- NULL
  if (bootstrap_draws > 0L) {
    set.seed(seed + state_count)
    storage <- vector("list", bootstrap_draws)
    progress_every <- max(1L, floor(bootstrap_draws / 20L))
    successful_draws <- 0L
    attempted_draws <- 0L
    maximum_attempts <- max(bootstrap_draws * 5L, bootstrap_draws + 100L)
    while (successful_draws < bootstrap_draws && attempted_draws < maximum_attempts) {
      attempted_draws <- attempted_draws + 1L
      simulated_blocks <- lapply(blocks, simulate_refit_block)
      if (any(vapply(simulated_blocks, is.null, logical(1)))) next
      candidate <- tryCatch(
        transition_metric_grid(simulated_blocks, cooperative, lower_return),
        error = function(condition) NULL
      )
      if (is.null(candidate) || any(!is.finite(candidate$value))) next
      successful_draws <- successful_draws + 1L
      storage[[successful_draws]] <- candidate %>%
        rename_broad_metrics() %>%
        mutate(iteration = successful_draws)
      if (successful_draws %% progress_every == 0L || successful_draws == bootstrap_draws) {
        message(
          "  ", state_count, "-state bootstrap ", successful_draws, "/",
          bootstrap_draws, " successful (", attempted_draws, " attempted)"
        )
      }
    }
    bootstrap_status <- tibble(
      requested_successful_draws = bootstrap_draws,
      successful_draws = successful_draws,
      attempted_draws = attempted_draws,
      rejected_draws = attempted_draws - successful_draws
    )
    atomic_write_csv(bootstrap_status, file.path(state_dir, "bootstrap_run_status.csv"))
    if (successful_draws != bootstrap_draws) {
      writeLines(
        c(
          "bootstrap_status: FAILED",
          paste0("requested_successful_draws: ", bootstrap_draws),
          paste0("successful_draws: ", successful_draws),
          paste0("attempted_draws: ", attempted_draws),
          paste0("rejected_draws: ", attempted_draws - successful_draws)
        ),
        file.path(state_dir, "FAILED_BOOTSTRAP.txt")
      )
      stop("Could not obtain the requested number of successful bootstrap draws")
    }
    bootstrap <- summarise_draws_finite(
      bind_rows(storage),
      point %>% dplyr::select(-state_count)
    )
    expected_family_sizes <- c(
      retention_and_entry = 18L,
      exit_to_lower_return = 9L,
      stationary_cooperation = 9L
    )
    observed_family_sizes <- table(factor(
      bootstrap$contrasts$correction_family,
      levels = names(expected_family_sizes)
    ))
    if (nrow(bootstrap$contrasts) != sum(expected_family_sizes) ||
        anyNA(bootstrap$contrasts$correction_family) ||
        any(!is.finite(bootstrap$contrasts$p_two_sided)) ||
        any(!is.finite(bootstrap$contrasts$p_holm_family)) ||
        !identical(as.integer(observed_family_sizes), unname(expected_family_sizes))) {
      stop("Sensitivity correction families are incomplete or malformed")
    }
    atomic_write_csv(
      bootstrap$levels,
      file.path(state_dir, sprintf("bootstrap_levels_B%d.csv", bootstrap_draws))
    )
    atomic_write_csv(
      bootstrap$contrasts,
      file.path(state_dir, sprintf("bootstrap_contrasts_B%d.csv", bootstrap_draws))
    )
  }

  result <- list(
    state_count = state_count,
    model_source = normalizePath(model_path),
    model_sha256 = sha256_file(model_path),
    state_dictionary = dictionary,
    lower_return_states = lower_return,
    cooperative_states = cooperative,
    cooperative_occupancy = occupancy,
    point = point,
    transition_summary_point = transition_summary_point,
    direct_vs_transition_summary_agreement = point_agreement,
    bootstrap = bootstrap,
    bootstrap_status = bootstrap_status,
    provenance = list(
      analysis_script_sha256 = sha256_file(
        file.path("scripts", "run_trustee_state_count_sensitivity.R")
      ),
      sourced_robustness_script_sha256 = sha256_file(
        file.path("scripts", "run_trustee_state_count_robustness.R")
      ),
      helper_sha256 = sha256_file(file.path("scripts", "hmm_cache_utils.R")),
      behavioral_data_sha256 = sha256_file(file.path("Data", "full_RTG_data.csv")),
      selection_artifact_sha256 = sha256_file(selection_artifact_rds),
      model_sha256 = sha256_file(model_path),
      requested_successful_draws = bootstrap_draws,
      successful_draws = if (is.null(bootstrap_status)) 0L else bootstrap_status$successful_draws[[1]],
      attempted_draws = if (is.null(bootstrap_status)) 0L else bootstrap_status$attempted_draws[[1]],
      rejected_draws = if (is.null(bootstrap_status)) 0L else bootstrap_status$rejected_draws[[1]]
    ),
    transition_orientation = "depmix_xi_to_from_transposed_to_from_to",
    seed = seed + state_count,
    generated_at = Sys.time()
  )
  atomic_save_rds(
    result,
    file.path(state_dir, sprintf("analysis_result_B%d.rds", bootstrap_draws))
  )
  all_results[[as.character(state_count)]] <- result
}

atomic_write_csv(
  bind_rows(lapply(all_results, `[[`, "state_dictionary")),
  file.path(output_dir, "combined_state_dictionary.csv")
)
atomic_write_csv(
  bind_rows(lapply(all_results, `[[`, "point")),
  file.path(output_dir, "combined_point_estimates.csv")
)
atomic_write_csv(
  bind_rows(lapply(all_results, function(result) {
    result$cooperative_occupancy$summary %>% mutate(state_count = result$state_count)
  })),
  file.path(output_dir, "combined_cooperative_occupancy_summary.csv")
)
atomic_write_csv(
  bind_rows(lapply(all_results, function(result) {
    result$cooperative_occupancy$pairs %>% mutate(state_count = result$state_count)
  })),
  file.path(output_dir, "combined_cooperative_occupancy_pairwise.csv")
)
if (bootstrap_draws > 0L) {
  atomic_write_csv(
    bind_rows(lapply(all_results, function(result) {
      result$bootstrap$contrasts %>% mutate(state_count = result$state_count)
    })),
    file.path(output_dir, sprintf("combined_bootstrap_contrasts_B%d.csv", bootstrap_draws))
  )
}

run_manifest <- list(
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  R_version = R.version.string,
  platform = R.version$platform,
  depmixS4_version = as.character(utils::packageVersion("depmixS4")),
  analysis_script_sha256 = sha256_file(
    file.path("scripts", "run_trustee_state_count_sensitivity.R")
  ),
  sourced_robustness_script_sha256 = sha256_file(
    file.path("scripts", "run_trustee_state_count_robustness.R")
  ),
  helper_sha256 = sha256_file(file.path("scripts", "hmm_cache_utils.R")),
  behavioral_data_sha256 = sha256_file(file.path("Data", "full_RTG_data.csv")),
  selection_artifact_path = normalizePath(selection_artifact_rds),
  selection_artifact_sha256 = sha256_file(selection_artifact_rds),
  state_counts = paste(state_counts, collapse = ","),
  bootstrap_draws_per_state = bootstrap_draws,
  base_seed = seed,
  total_successful_draws = sum(vapply(all_results, function(result) {
    if (is.null(result$bootstrap_status)) 0 else result$bootstrap_status$successful_draws[[1]]
  }, numeric(1))),
  total_attempted_draws = sum(vapply(all_results, function(result) {
    if (is.null(result$bootstrap_status)) 0 else result$bootstrap_status$attempted_draws[[1]]
  }, numeric(1))),
  total_rejected_draws = sum(vapply(all_results, function(result) {
    if (is.null(result$bootstrap_status)) 0 else result$bootstrap_status$rejected_draws[[1]]
  }, numeric(1))),
  transition_orientation = "depmix_xi_to_from_transposed_to_from_to"
)
atomic_write_csv(
  enframe(unlist(run_manifest), name = "field", value = "value"),
  file.path(output_dir, "RUN_MANIFEST.csv")
)
message("State-count sensitivity analysis complete")
