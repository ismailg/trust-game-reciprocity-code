#!/usr/bin/env Rscript

# Five-state trustee analysis with State 1 clearly low return, State 2 near one
# third of the multiplied investment, and States 3-5 cooperative. This is
# standalone and does not source, edit, or render the manuscript.

Sys.setenv(TRUSTEE_ROBUSTNESS_LIBRARY_ONLY = "1")
source(file.path("scripts", "run_trustee_state_count_robustness.R"))
Sys.unsetenv("TRUSTEE_ROBUSTNESS_LIBRARY_ONLY")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[length(hit)]])
}
bootstrap_draws <- as.integer(get_arg("--bootstrap", "1000"))
seed <- as.integer(get_arg("--seed", "11"))
model_rds <- get_arg("--model-rds", "")
technical_checks_csv <- get_arg("--technical-checks", "")
reference_result_rds <- get_arg("--reference-result-rds", "")
output_dir <- get_arg(
  "--output-dir",
  file.path("results", "HMM", "state_count_robustness", "states5_state2_neutral")
)
if (!is.finite(bootstrap_draws) || bootstrap_draws < 0) stop("--bootstrap must be >= 0")
if (!nzchar(model_rds)) stop("--model-rds must identify the approved validated artifact")
if (!file.exists(model_rds)) stop("Validated model file does not exist: ", model_rds)
if (!nzchar(technical_checks_csv)) {
  technical_checks_csv <- file.path(dirname(model_rds), "technical_checks.csv")
}
if (!file.exists(technical_checks_csv)) {
  stop("Technical-check file does not exist: ", technical_checks_csv)
}
if (dir.exists(output_dir) && length(list.files(
  output_dir, all.files = TRUE, no.. = TRUE
)) > 0L) {
  stop("Output directory must be new or empty: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

groups <- list(
  clearly_low_return = 1L,
  near_one_third = 2L,
  cooperative = 3:5
)

weighted_probability <- function(transition_matrix, stationary, from, to) {
  denominator <- sum(stationary[from])
  numerator <- sum(
    stationary[from] * rowSums(transition_matrix[from, to, drop = FALSE])
  )
  if (denominator > 0) numerator / denominator else NA_real_
}

neutral_metric_values <- function(transition_matrix) {
  stationary <- stationary_distribution(transition_matrix)
  c(
    retention = weighted_probability(
      transition_matrix, stationary, groups$cooperative, groups$cooperative
    ),
    entry_from_clearly_low_return = weighted_probability(
      transition_matrix, stationary, groups$clearly_low_return, groups$cooperative
    ),
    entry_from_near_one_third = weighted_probability(
      transition_matrix, stationary, groups$near_one_third, groups$cooperative
    ),
    exit_to_clearly_low_return = weighted_probability(
      transition_matrix, stationary, groups$cooperative, groups$clearly_low_return
    ),
    exit_to_near_one_third = weighted_probability(
      transition_matrix, stationary, groups$cooperative, groups$near_one_third
    ),
    stationary_cooperation = sum(stationary[groups$cooperative]),
    stationary_clearly_low_return = sum(stationary[groups$clearly_low_return]),
    stationary_near_one_third = sum(stationary[groups$near_one_third])
  )
}

neutral_metric_grid <- function(blocks) {
  purrr::map_dfr(diagnoses, function(diagnosis) {
    purrr::map_dfr(investment_grid, function(investment) {
      transition_matrix <- compose_transition_matrix(blocks, investment, diagnosis)
      values <- neutral_metric_values(transition_matrix)
      tibble(
        diagnosis = diagnosis,
        investment = investment,
        metric = names(values),
        value = as.numeric(values)
      )
    })
  })
}

neutral_correction_family <- function(metric) {
  case_when(
    metric %in% c(
      "retention", "entry_from_clearly_low_return", "entry_from_near_one_third"
    ) ~ "retention_and_entry",
    metric %in% c(
      "exit_to_clearly_low_return", "exit_to_near_one_third"
    ) ~ "exit_destination",
    metric == "stationary_cooperation" ~ "stationary_cooperation",
    TRUE ~ NA_character_
  )
}

summarise_neutral_draws <- function(draws, point) {
  levels <- draws %>%
    group_by(diagnosis, investment, metric) %>%
    summarise(
      bootstrap_mean = mean(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    left_join(
      point %>% rename(point_estimate = value),
      by = c("diagnosis", "investment", "metric")
    )

  contrast_draws <- draws %>%
    filter(metric %in% c(
      "retention",
      "entry_from_clearly_low_return",
      "entry_from_near_one_third",
      "exit_to_clearly_low_return",
      "exit_to_near_one_third",
      "stationary_cooperation"
    )) %>%
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
    filter(metric %in% c(
      "retention",
      "entry_from_clearly_low_return",
      "entry_from_near_one_third",
      "exit_to_clearly_low_return",
      "exit_to_near_one_third",
      "stationary_cooperation"
    )) %>%
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
    mutate(correction_family = neutral_correction_family(metric)) %>%
    group_by(correction_family) %>%
    mutate(p_holm_family = p.adjust(p_two_sided, method = "holm")) %>%
    ungroup() %>%
    left_join(point_contrasts, by = c("investment", "metric", "contrast"))
  list(levels = levels, contrasts = contrasts, draws = draws)
}

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

prepare_decoded_data <- function(model) {
  decoded <- trustee_data_base
  decoded$posterior_state <- posterior(model, type = "local")
  decoded
}

decoded_neutral_analyses <- function(data) {
  transitions <- data %>%
    arrange(subject_ID, round) %>%
    group_by(subject_ID) %>%
    mutate(
      next_state = lead(posterior_state),
      next_investment_analysis = lead(investment),
      current_gap = return - fair_return
    ) %>%
    ungroup() %>%
    filter(!is.na(next_state), !is.na(next_investment_analysis))

  investment_mean <- mean(transitions$next_investment_analysis)
  investment_sd <- sd(transitions$next_investment_analysis)
  gap_mean <- mean(transitions$current_gap, na.rm = TRUE)
  gap_sd <- sd(transitions$current_gap, na.rm = TRUE)
  scaled_grid <- (investment_grid - investment_mean) / investment_sd

  common_mutate <- function(input) {
    output <- input %>% mutate(
      diag = factor(group, levels = diagnoses),
      inv_sc = (next_investment_analysis - investment_mean) / investment_sd,
      gap_sc = if_else(is.na(current_gap), NA_real_, (current_gap - gap_mean) / gap_sd)
    )
    attr(output, "investment_scaled_grid") <- scaled_grid
    output
  }

  retention <- transitions %>%
    filter(posterior_state %in% groups$cooperative) %>%
    mutate(stay_cooperative = as.integer(next_state %in% groups$cooperative)) %>%
    common_mutate()
  clearly_low_entry <- transitions %>%
    filter(posterior_state == groups$clearly_low_return) %>%
    mutate(enter_cooperative = as.integer(next_state %in% groups$cooperative)) %>%
    common_mutate()
  near_one_third_entry <- transitions %>%
    filter(posterior_state == groups$near_one_third) %>%
    mutate(enter_cooperative = as.integer(next_state %in% groups$cooperative)) %>%
    common_mutate()
  clearly_low_exit <- transitions %>%
    filter(posterior_state %in% groups$cooperative) %>%
    mutate(
      destination_event = as.integer(next_state == groups$clearly_low_return),
      start_state = factor(posterior_state)
    ) %>%
    common_mutate()
  near_one_third_exit <- transitions %>%
    filter(posterior_state %in% groups$cooperative) %>%
    mutate(
      destination_event = as.integer(next_state == groups$near_one_third),
      start_state = factor(posterior_state)
    ) %>%
    common_mutate()

  analyses <- list(
    safe_glmm_analysis(retention, "stay_cooperative", "retention", FALSE),
    safe_glmm_analysis(
      clearly_low_entry, "enter_cooperative", "entry_from_clearly_low_return", FALSE
    ),
    safe_glmm_analysis(
      near_one_third_entry, "enter_cooperative", "entry_from_near_one_third", FALSE
    ),
    safe_glmm_analysis(
      clearly_low_exit, "destination_event", "exit_to_clearly_low_return", TRUE
    ),
    safe_glmm_analysis(
      near_one_third_exit, "destination_event", "exit_to_near_one_third", TRUE
    )
  )
  list(
    status = bind_rows(lapply(analyses, `[[`, "status")),
    predictions = bind_rows(lapply(analyses, `[[`, "predictions")) %>%
      mutate(investment = investment_grid[match(inv_sc, scaled_grid)]),
    contrasts = bind_rows(lapply(analyses, `[[`, "contrasts")) %>%
      mutate(investment = investment_grid[match(inv_sc, scaled_grid)]),
    coefficients = bind_rows(lapply(analyses, `[[`, "coefficients"))
  )
}

message("Preparing approved five-state model with State 2 near one third")
validated <- readRDS(model_rds)
required_artifact_fields <- c(
  "best_model", "selected_state_count", "selected_fit_label", "gate", "provenance"
)
if (!is.list(validated) || !all(required_artifact_fields %in% names(validated)) ||
    !identical(as.integer(validated$selected_state_count), 5L) ||
    !identical(as.character(validated$selected_fit_label), "perturbed_warm_02") ||
    !is.data.frame(validated$gate) || nrow(validated$gate) == 0L ||
    anyDuplicated(validated$gate$check) || !all(validated$gate$passed)) {
  stop("--model-rds is not the approved, fully gated five-state artifact")
}
technical_checks <- readr::read_csv(technical_checks_csv, show_col_types = FALSE)
expected_technical_checks <- c(
  "sample_has_894_sequences",
  "selected_fit_reports_convergence",
  "sample_has_8940_trials",
  "all_sequences_have_10_rounds",
  "selection_gate_passed",
  "model_table_has_one_named_row_per_state_count",
  "five_states_minimizes_trial_BIC",
  "five_states_minimizes_sequence_BIC",
  "states_are_ordered_after_fitting",
  "state_mapping_is_one_exploitative_one_near_break_even_three_cooperative",
  "posterior_has_expected_shape",
  "no_negligible_state_occupancy",
  "emissions_normalize_for_investments_zero_to_twenty",
  "stored_transition_probabilities_sum_to_one_by_from_state",
  "direct_transition_predictions_match_package_probabilities",
  "expected_transition_orientation_matches_posterior_from_state",
  "five_state_group_transition_test_is_valid_and_below_0.05"
)
if (!all(c("check", "passed", "detail") %in% names(technical_checks)) ||
    !setequal(technical_checks$check, expected_technical_checks) ||
    nrow(technical_checks) != length(expected_technical_checks) ||
    anyDuplicated(technical_checks$check) ||
    anyNA(technical_checks$passed) || !all(technical_checks$passed)) {
  stop("The canonical technical checks are missing, duplicated, or failed")
}
if (!identical(
  normalizePath(dirname(model_rds), mustWork = TRUE),
  normalizePath(dirname(technical_checks_csv), mustWork = TRUE)
)) {
  stop("Technical checks must be adjacent to the validated model artifact")
}
selected <- list(model = validated$best_model, source = normalizePath(model_rds))
if (nstates(selected$model) != 5L) stop("Validated artifact does not contain five states")
model <- order_mod_vtdgaus(selected$model)
means <- expected_state_means(model)
if (!all(diff(means) > 0) ||
    means[1] > 0.30 || means[2] < 0.30 || means[2] >= 0.40 ||
    means[2] - means[1] < 0.05 ||
    !identical(groups$cooperative, which(means >= 0.40))) {
  stop(paste(
    "Five-state means do not support one clearly low-return state,",
    "one near-one-third state, and three cooperative states"
  ))
}

state_dictionary <- tibble(
  state_count = 5L,
  state = 1:5,
  expected_return_fraction = means,
  category = c("clearly_low_return", "near_one_third", rep("cooperative", 3))
)
atomic_write_csv(state_dictionary, file.path(output_dir, "state_dictionary.csv"))

emission_grid <- purrr::map_dfr(seq_len(5L), function(state) {
  return_fraction <- seq(0, 1, length.out = 61L)
  parameters <- getpars(model@response[[state]][[1]])
  tibble(
    state = state,
    return_fraction = return_fraction,
    probability = as.numeric(dens(vtdgaus(
      return_fraction,
      yield = rep(60, length(return_fraction)),
      pstart = parameters
    )))
  )
})
atomic_write_csv(emission_grid, file.path(output_dir, "state_emission_profiles.csv"))

smoothed_probabilities <- as_tibble(
  posterior(model, type = "smoothing"),
  .name_repair = ~ paste0("state_", seq_along(.x))
)
probability_weighted_state_occupancy <- bind_cols(
  trustee_data_base %>% dplyr::select(subject_ID, group, type),
  smoothed_probabilities
) %>%
  pivot_longer(
    starts_with("state_"),
    names_to = "state_label",
    values_to = "posterior_probability"
  ) %>%
  mutate(state = as.integer(sub("state_", "", state_label))) %>%
  group_by(subject_ID, group, type, state) %>%
  summarise(occupancy = mean(posterior_probability), .groups = "drop")
probability_weighted_category_occupancy <- probability_weighted_state_occupancy %>%
  mutate(category = case_when(
    state == 1L ~ "clearly_low_return",
    state == 2L ~ "near_one_third",
    TRUE ~ "cooperative"
  )) %>%
  group_by(subject_ID, group, type, category) %>%
  summarise(occupancy = sum(occupancy), .groups = "drop")
probability_weighted_category_summary <- probability_weighted_category_occupancy %>%
  group_by(group, category) %>%
  summarise(
    participants = n(),
    mean_occupancy = mean(occupancy),
    standard_error = sd(occupancy) / sqrt(n()),
    .groups = "drop"
  )
probability_weighted_category_pairs <- probability_weighted_category_occupancy %>%
  group_split(category) %>%
  purrr::map_dfr(function(category_data) {
    fitted <- lm(occupancy ~ group, data = category_data)
    as_tibble(as.data.frame(summary(
      contrast(emmeans(fitted, ~ group), "pairwise", adjust = "holm"),
      infer = c(TRUE, TRUE)
    ))) %>%
      transmute(
        category = unique(category_data$category),
        contrast,
        difference = estimate,
        lower = lower.CL,
        upper = upper.CL,
        p_holm = p.value
      )
  })
atomic_write_csv(
  probability_weighted_state_occupancy,
  file.path(output_dir, "probability_weighted_state_occupancy.csv")
)
atomic_write_csv(
  probability_weighted_category_occupancy,
  file.path(output_dir, "probability_weighted_category_occupancy.csv")
)
atomic_write_csv(
  probability_weighted_category_summary,
  file.path(output_dir, "probability_weighted_category_occupancy_summary.csv")
)
atomic_write_csv(
  probability_weighted_category_pairs,
  file.path(output_dir, "probability_weighted_category_occupancy_pairwise.csv")
)

decoded_data <- prepare_decoded_data(model)
occupancy <- decoded_occupancy(
  decoded_data,
  list(exploitative = groups$clearly_low_return, cooperative = groups$cooperative)
)
occupancy <- lapply(occupancy, function(table) {
  if (!"category" %in% names(table)) return(table)
  table %>% mutate(
    category = recode(
      as.character(category),
      exploitative = "clearly_low_return",
      intermediate = "near_one_third"
    )
  )
})
decoded <- decoded_neutral_analyses(decoded_data)
atomic_write_csv(occupancy$summary, file.path(output_dir, "decoded_category_occupancy_summary.csv"))
atomic_write_csv(occupancy$tests, file.path(output_dir, "decoded_category_occupancy_tests.csv"))
atomic_write_csv(occupancy$pairs, file.path(output_dir, "decoded_category_occupancy_pairwise.csv"))
atomic_write_csv(decoded$status, file.path(output_dir, "decoded_glmm_status.csv"))
atomic_write_csv(decoded$predictions, file.path(output_dir, "decoded_glmm_predictions.csv"))
atomic_write_csv(decoded$contrasts, file.path(output_dir, "decoded_glmm_pairwise.csv"))
atomic_write_csv(decoded$coefficients, file.path(output_dir, "decoded_glmm_coefficients.csv"))

message("Calculating transition probabilities directly from the fitted model")
point <- purrr::map_dfr(diagnoses, function(diagnosis) {
  purrr::map_dfr(investment_grid, function(investment) {
    transition_matrix <- direct_transition_matrix(model, investment, diagnosis)
    values <- neutral_metric_values(transition_matrix)
    tibble(
      diagnosis = diagnosis,
      investment = investment,
      metric = names(values),
      value = as.numeric(values)
    )
  })
})
atomic_write_csv(point, file.path(output_dir, "transition_metric_point_estimates.csv"))

message("Fitting corrected transition-summary blocks for conditional uncertainty")
xi <- expected_transitions_xi(model, trustee_data_base)
frames <- build_multinom_frames(model, trustee_data_base, xi)
blocks <- lapply(frames, fit_block_multinom)
if (any(vapply(blocks, is.null, logical(1)))) stop("At least one transition block failed")

transition_summary_point <- neutral_metric_grid(blocks)
point_agreement <- point %>%
  rename(direct_fitted_model = value) %>%
  left_join(
    transition_summary_point %>% rename(transition_summary_model = value),
    by = c("diagnosis", "investment", "metric")
  ) %>%
  mutate(absolute_difference = abs(direct_fitted_model - transition_summary_model))
if (nrow(point_agreement) != nrow(point) ||
    any(!is.finite(point_agreement$absolute_difference)) ||
    max(point_agreement$absolute_difference) > 0.005) {
  stop("Transition-summary model does not adequately reproduce the fitted HMM")
}
atomic_write_csv(
  transition_summary_point,
  file.path(output_dir, "transition_summary_model_point_estimates.csv")
)
atomic_write_csv(
  point_agreement,
  file.path(output_dir, "direct_vs_transition_summary_agreement.csv")
)

point_check <- NULL
reference_result <- NULL
if (nzchar(reference_result_rds)) {
  if (!file.exists(reference_result_rds)) {
    stop("Reference result does not exist: ", reference_result_rds)
  }
  reference_result <- readRDS(reference_result_rds)
  reference_point <- reference_result$point %>%
    dplyr::select(diagnosis, investment, metric, value)
  point_invariants <- point %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(
      reference_exit_reconstructed =
        exit_to_clearly_low_return + exit_to_near_one_third,
      reference_entry_reconstructed =
        (stationary_clearly_low_return * entry_from_clearly_low_return +
           stationary_near_one_third * entry_from_near_one_third) /
        (stationary_clearly_low_return + stationary_near_one_third)
    ) %>%
    dplyr::select(
      diagnosis, investment, retention, stationary_cooperation,
      reference_exit_reconstructed, reference_entry_reconstructed
    )
  reference_point_wide <- reference_point %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    dplyr::select(
      diagnosis, investment, retention, stationary_cooperation,
      exploitative_exit, entry
    )
  point_check <- point_invariants %>%
    left_join(
      reference_point_wide,
      by = c("diagnosis", "investment"),
      suffix = c("_new", "_reference")
    ) %>%
    transmute(
      diagnosis,
      investment,
      retention_difference = abs(retention_new - retention_reference),
      stationary_difference = abs(
        stationary_cooperation_new - stationary_cooperation_reference
      ),
      exit_decomposition_difference = abs(
        reference_exit_reconstructed - exploitative_exit
      ),
      entry_decomposition_difference = abs(reference_entry_reconstructed - entry)
    )
  point_difference_columns <- c(
    "retention_difference", "stationary_difference",
    "exit_decomposition_difference", "entry_decomposition_difference"
  )
  if (max(as.matrix(point_check[point_difference_columns])) > 1e-8) {
    stop("Point estimates fail checks against the explicitly supplied reference run")
  }
  atomic_write_csv(
    point_check,
    file.path(output_dir, "reference_run_point_equivalence_check.csv")
  )
}

bootstrap <- NULL
draw_check <- NULL
bootstrap_status <- NULL
if (bootstrap_draws > 0) {
  set.seed(seed)
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
      neutral_metric_grid(simulated_blocks),
      error = function(condition) NULL
    )
    if (is.null(candidate) || any(!is.finite(candidate$value))) next
    successful_draws <- successful_draws + 1L
    storage[[successful_draws]] <- candidate %>% mutate(iteration = successful_draws)
    if (successful_draws %% progress_every == 0 || successful_draws == bootstrap_draws) {
      message(
        "  neutral-state bootstrap ", successful_draws, "/", bootstrap_draws,
        " successful (", attempted_draws, " attempted)"
      )
    }
  }
  bootstrap_status <- tibble(
    requested_successful_draws = bootstrap_draws,
    successful_draws = successful_draws,
    attempted_draws = attempted_draws,
    rejected_draws = attempted_draws - successful_draws
  )
  atomic_write_csv(bootstrap_status, file.path(output_dir, "bootstrap_run_status.csv"))
  if (successful_draws != bootstrap_draws) {
    writeLines(
      c(
        "bootstrap_status: FAILED",
        paste0("requested_successful_draws: ", bootstrap_draws),
        paste0("successful_draws: ", successful_draws),
        paste0("attempted_draws: ", attempted_draws),
        paste0("rejected_draws: ", attempted_draws - successful_draws)
      ),
      file.path(output_dir, "FAILED_BOOTSTRAP.txt")
    )
    stop("Could not obtain the requested number of successful bootstrap draws")
  }
  draws <- bind_rows(storage)
  expected_rows <- bootstrap_draws * length(diagnoses) * length(investment_grid) * 8L
  if (nrow(draws) != expected_rows || any(!is.finite(draws$value))) {
    stop("Neutral-state bootstrap draws are incomplete or non-finite")
  }
  bootstrap <- summarise_neutral_draws(draws, point)
  expected_family_sizes <- c(
    exit_destination = 18L,
    retention_and_entry = 27L,
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
    stop("Bootstrap correction families are incomplete or malformed")
  }

  if (!is.null(reference_result)) {
    reference_draws <- reference_result$bootstrap$draws %>%
      filter(iteration <= bootstrap_draws) %>%
      dplyr::select(
        iteration, diagnosis, investment, metric, reference_value = value
      )
    new_wide <- draws %>%
      pivot_wider(names_from = metric, values_from = value) %>%
      mutate(
        exploitative_exit =
          exit_to_clearly_low_return + exit_to_near_one_third,
        entry =
          (stationary_clearly_low_return * entry_from_clearly_low_return +
             stationary_near_one_third * entry_from_near_one_third) /
          (stationary_clearly_low_return + stationary_near_one_third)
      ) %>%
      dplyr::select(
        iteration, diagnosis, investment, retention, entry,
        exploitative_exit, stationary_cooperation
      ) %>%
      pivot_longer(
        cols = c(retention, entry, exploitative_exit, stationary_cooperation),
        names_to = "metric",
        values_to = "new_value"
      )
    draw_check <- new_wide %>%
      left_join(
        reference_draws,
        by = c("iteration", "diagnosis", "investment", "metric")
      ) %>%
      summarise(
        compared_rows = n(),
        missing_reference_rows = sum(is.na(reference_value)),
        maximum_absolute_difference = max(
          abs(new_value - reference_value),
          na.rm = TRUE
        )
      )
    if (draw_check$compared_rows != bootstrap_draws * 3L * 3L * 4L ||
        draw_check$missing_reference_rows != 0L ||
        draw_check$maximum_absolute_difference > 1e-8) {
      stop("Bootstrap draws fail checks against the explicitly supplied reference run")
    }
  }
  atomic_write_csv(
    bootstrap$levels,
    file.path(output_dir, sprintf("bootstrap_levels_B%d.csv", bootstrap_draws))
  )
  atomic_write_csv(
    bootstrap$contrasts,
    file.path(output_dir, sprintf("bootstrap_contrasts_B%d.csv", bootstrap_draws))
  )
  if (!is.null(draw_check)) {
    atomic_write_csv(
      draw_check,
      file.path(
        output_dir,
        sprintf("reference_run_draw_equivalence_check_B%d.csv", bootstrap_draws)
      )
    )
  }
}

final_bootstrap_status <- if (is.null(bootstrap_status)) {
  tibble(
    requested_successful_draws = bootstrap_draws,
    successful_draws = 0L,
    attempted_draws = 0L,
    rejected_draws = 0L
  )
} else {
  bootstrap_status
}
run_provenance <- list(
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  R_version = R.version.string,
  platform = R.version$platform,
  depmixS4_version = as.character(utils::packageVersion("depmixS4")),
  analysis_script_sha256 = sha256_file(
    file.path("scripts", "run_trustee_state5_neutral_state2.R")
  ),
  sourced_robustness_script_sha256 = sha256_file(
    file.path("scripts", "run_trustee_state_count_robustness.R")
  ),
  helper_sha256 = sha256_file(file.path("scripts", "hmm_cache_utils.R")),
  behavioral_data_sha256 = sha256_file(file.path("Data", "full_RTG_data.csv")),
  demographics_data_sha256 = sha256_file(file.path("Data", "demographics.csv")),
  validated_model_artifact_sha256 = sha256_file(model_rds),
  validated_model_artifact_path = normalizePath(model_rds),
  technical_checks_sha256 = sha256_file(technical_checks_csv),
  technical_checks_path = normalizePath(technical_checks_csv),
  selected_state_count = validated$selected_state_count,
  selected_fit_label = validated$selected_fit_label,
  technical_checks_passed = all(technical_checks$passed),
  technical_check_count = nrow(technical_checks),
  bootstrap_seed = seed,
  requested_successful_draws = final_bootstrap_status$requested_successful_draws[[1]],
  successful_draws = final_bootstrap_status$successful_draws[[1]],
  attempted_draws = final_bootstrap_status$attempted_draws[[1]],
  rejected_draws = final_bootstrap_status$rejected_draws[[1]],
  transition_orientation = "depmix_xi_to_from_transposed_to_from_to"
)
atomic_write_csv(
  enframe(unlist(run_provenance), name = "field", value = "value"),
  file.path(output_dir, "RUN_MANIFEST.csv")
)

result <- list(
  state_count = 5L,
  state_groups = groups,
  state_dictionary = state_dictionary,
  emission_grid = emission_grid,
  probability_weighted_occupancy = list(
    state = probability_weighted_state_occupancy,
    category = probability_weighted_category_occupancy,
    summary = probability_weighted_category_summary,
    pairs = probability_weighted_category_pairs
  ),
  occupancy = occupancy,
  decoded = decoded,
  point = point,
  transition_summary_point = transition_summary_point,
  direct_vs_transition_summary_agreement = point_agreement,
  bootstrap = bootstrap,
  point_equivalence_check = point_check,
  draw_equivalence_check = draw_check,
  bootstrap_status = final_bootstrap_status,
  seed = seed,
  model_source = selected$source,
  technical_checks_source = normalizePath(technical_checks_csv),
  provenance = run_provenance,
  transition_orientation = "depmix_xi_to_from_transposed_to_from_to",
  manuscript_modified_or_rendered = FALSE,
  generated_at = Sys.time()
)
atomic_save_rds(
  result,
  file.path(output_dir, sprintf("analysis_result_B%d.rds", bootstrap_draws))
)
message("Five-state clearly-low / near-one-third analysis complete")
