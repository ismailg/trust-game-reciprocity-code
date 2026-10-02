#!/usr/bin/env Rscript

# Focused checks of whether the decoded five-state trustee findings persist
# (1) among online participants and (2) after adjusting for age and gender.

Sys.setenv(TRUSTEE_ROBUSTNESS_LIBRARY_ONLY = "1")
source(file.path("scripts", "run_trustee_state_count_robustness.R"))
Sys.unsetenv("TRUSTEE_ROBUSTNESS_LIBRARY_ONLY")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[length(hit)]])
}
model_rds <- get_arg(
  "--model-rds",
  file.path(
    "results", "HMM", "corrected_submission",
    "model_selection_final_2026-07-21_v2",
    "validated_five_state_model.rds"
  )
)
technical_checks_csv <- get_arg("--technical-checks", "")
output_dir <- get_arg(
  "--output-dir",
  file.path("results", "HMM", "corrected_submission", "covariate_sensitivity")
)
if (!file.exists(model_rds)) stop("Validated model file does not exist: ", model_rds)
if (!nzchar(technical_checks_csv)) {
  technical_checks_csv <- file.path(dirname(model_rds), "technical_checks.csv")
}
if (!file.exists(technical_checks_csv)) {
  stop("Technical-check file does not exist: ", technical_checks_csv)
}
if (!identical(
  normalizePath(dirname(model_rds), mustWork = TRUE),
  normalizePath(dirname(technical_checks_csv), mustWork = TRUE)
)) {
  stop("Technical checks must be adjacent to the validated model artifact")
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

artifact <- readRDS(model_rds)
if (!is.list(artifact) ||
    !all(c(
      "best_model", "selected_state_count", "selected_fit_label", "gate", "provenance"
    ) %in% names(artifact)) ||
    !identical(as.integer(artifact$selected_state_count), 5L) ||
    !identical(as.character(artifact$selected_fit_label), "perturbed_warm_02") ||
    !is.data.frame(artifact$gate) || nrow(artifact$gate) == 0L ||
    anyDuplicated(artifact$gate$check) || !all(artifact$gate$passed)) {
  stop("Model artifact is not the approved fully gated five-state result")
}
technical_checks <- readr::read_csv(technical_checks_csv, show_col_types = FALSE)
if (!all(c("check", "passed", "detail") %in% names(technical_checks)) ||
    nrow(technical_checks) != 17L || anyDuplicated(technical_checks$check) ||
    anyNA(technical_checks$passed) || !all(technical_checks$passed)) {
  stop("Canonical technical checks are missing, duplicated, or failed")
}
model <- order_mod_vtdgaus(artifact$best_model)
if (nstates(model) != 5L) stop("Covariate sensitivity requires the five-state model")
means <- expected_state_means(model)
if (!all(diff(means) > 0) || means[1] > 0.30 || means[2] < 0.30 ||
    means[2] >= 0.40 || means[2] - means[1] < 0.05 ||
    !identical(which(means >= 0.40), 3:5)) {
  stop("Unexpected clearly-low / near-one-third / cooperative classification")
}

demographics <- read.csv(
  file.path("Data", "demographics.csv"),
  encoding = "UTF-8",
  check.names = FALSE
) %>%
  transmute(
    subject_ID = as.character(ID),
    age = suppressWarnings(as.numeric(Age)),
    gender = recode(
      as.character(Gender),
      `1` = "Male",
      `2` = "Female",
      .default = NA_character_
    )
  ) %>%
  mutate(
    age = if_else(!is.na(age) & age > 0 & age < 100, age, NA_real_),
    gender = factor(gender, levels = c("Male", "Female"))
  )
demographic_conflicts <- demographics %>%
  group_by(subject_ID) %>%
  summarise(
    distinct_age_values = n_distinct(age, na.rm = TRUE),
    distinct_gender_values = n_distinct(gender, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(distinct_age_values > 1L | distinct_gender_values > 1L)
if (nrow(demographic_conflicts) > 0L) {
  stop("Demographic file contains conflicting duplicate participant rows")
}
demographics <- demographics %>% distinct(subject_ID, age, gender)

decoded <- trustee_data_base
decoded$posterior_state <- posterior(model, type = "local")
smoothed_probabilities <- as.matrix(posterior(model, type = "smoothing"))
if (!identical(dim(smoothed_probabilities), c(nrow(decoded), 5L))) {
  stop("Unexpected smoothed-posterior dimensions")
}
decoded$cooperative_probability <- rowSums(smoothed_probabilities[, 3:5, drop = FALSE])
expected_decoded_rows <- nrow(decoded)
decoded <- decoded %>%
  mutate(subject_ID_character = as.character(subject_ID)) %>%
  left_join(
    demographics,
    by = c("subject_ID_character" = "subject_ID"),
    relationship = "many-to-one"
  )
if (nrow(decoded) != expected_decoded_rows) {
  stop("Demographic join changed the number of trustee trials")
}

transitions <- decoded %>%
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
age_mean <- mean(transitions$age, na.rm = TRUE)
age_sd <- sd(transitions$age, na.rm = TRUE)
scaled_investments <- (investment_grid - investment_mean) / investment_sd

transitions <- transitions %>% mutate(
  diag = factor(group, levels = diagnoses),
  inv_sc = (next_investment_analysis - investment_mean) / investment_sd,
  gap_sc = (current_gap - gap_mean) / gap_sd,
  age_sc = (age - age_mean) / age_sd
)

event_data <- list(
  retention = transitions %>%
    filter(posterior_state %in% 3:5) %>%
    mutate(outcome = as.integer(next_state %in% 3:5)),
  entry_from_clearly_low_return = transitions %>%
    filter(posterior_state == 1L) %>%
    mutate(outcome = as.integer(next_state %in% 3:5)),
  entry_from_near_one_third = transitions %>%
    filter(posterior_state == 2L) %>%
    mutate(outcome = as.integer(next_state %in% 3:5)),
  exit_to_clearly_low_return = transitions %>%
    filter(posterior_state %in% 3:5) %>%
    mutate(outcome = as.integer(next_state == 1L), start_state = factor(posterior_state)),
  exit_to_near_one_third = transitions %>%
    filter(posterior_state %in% 3:5) %>%
    mutate(outcome = as.integer(next_state == 2L), start_state = factor(posterior_state))
)

fit_event <- function(data, event, sensitivity, adjusted, include_start_state) {
  analysis_data <- data
  if (startsWith(sensitivity, "online_")) {
    analysis_data <- filter(analysis_data, type == "WEB")
  }
  if (adjusted || sensitivity == "online_complete_case") {
    analysis_data <- filter(analysis_data, !is.na(age_sc), !is.na(gender))
  }
  formula_text <- paste0(
    "outcome ~ diag * inv_sc + gap_sc",
    if (include_start_state) " + start_state" else "",
    if (adjusted) " + age_sc + gender" else "",
    " + (1 | subject_ID)"
  )
  warnings <- character()
  fitted <- tryCatch(
    withCallingHandlers(
      glmer(
        as.formula(formula_text),
        data = analysis_data,
        family = binomial,
        control = glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 1e6),
          check.conv.singular = "ignore"
        )
      ),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) condition
  )
  if (inherits(fitted, "error")) {
    return(list(
      status = tibble(
        event, sensitivity, status = "failed", detail = conditionMessage(fitted),
        observations = nrow(analysis_data), subjects = n_distinct(analysis_data$subject_ID)
      ),
      contrasts = tibble()
    ))
  }
  at_values <- list(inv_sc = scaled_investments, gap_sc = 0)
  if (adjusted) at_values$age_sc <- 0
  marginal <- withCallingHandlers(
    emmeans(
      fitted,
      ~ diag | inv_sc,
      at = at_values,
      weights = "proportional"
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  pair_table <- as.data.frame(withCallingHandlers(
    summary(
      contrast(marginal, "pairwise", by = "inv_sc", adjust = "none"),
      infer = c(TRUE, TRUE),
      adjust = "none"
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  ))
  lower_name <- intersect(c("asymp.LCL", "lower.CL"), names(pair_table))[[1]]
  upper_name <- intersect(c("asymp.UCL", "upper.CL"), names(pair_table))[[1]]
  contrasts <- as_tibble(pair_table) %>% transmute(
    event,
    sensitivity,
    investment = investment_grid[match(inv_sc, scaled_investments)],
    contrast,
    log_odds_difference = estimate,
    standard_error = SE,
    odds_ratio = exp(estimate),
    lower_odds_ratio = exp(.data[[lower_name]]),
    upper_odds_ratio = exp(.data[[upper_name]]),
    p_raw = p.value
  ) %>%
    group_by(investment) %>%
    mutate(
      p_holm_within_investment = p.adjust(p_raw, method = "holm")
    ) %>%
    ungroup()
  details <- c(warnings, fitted@optinfo$conv$lme4$messages)
  list(
    status = tibble(
      event,
      sensitivity,
      status = if (length(details) == 0L) "ok" else "warning",
      detail = paste(unique(details), collapse = " | "),
      observations = nrow(analysis_data),
      subjects = n_distinct(analysis_data$subject_ID),
      event_rate = mean(analysis_data$outcome),
      singular = isSingular(fitted)
    ),
    contrasts = contrasts
  )
}

analyses <- list()
for (event in names(event_data)) {
  include_start_state <- startsWith(event, "exit_")
  analyses[[paste(event, "online", sep = "_")]] <- fit_event(
    event_data[[event]], event, "online_only", FALSE, include_start_state
  )
  analyses[[paste(event, "online_complete_case", sep = "_")]] <- fit_event(
    event_data[[event]], event, "online_complete_case", FALSE, include_start_state
  )
  analyses[[paste(event, "online_adjusted", sep = "_")]] <- fit_event(
    event_data[[event]], event, "online_age_gender_adjusted", TRUE, include_start_state
  )
  analyses[[paste(event, "adjusted", sep = "_")]] <- fit_event(
    event_data[[event]], event, "all_modes_age_gender_adjusted", TRUE, include_start_state
  )
}

status <- bind_rows(lapply(analyses, `[[`, "status"))
contrasts <- bind_rows(lapply(analyses, `[[`, "contrasts"))
atomic_write_csv(status, file.path(output_dir, "decoded_transition_sensitivity_status.csv"))
atomic_write_csv(contrasts, file.path(output_dir, "decoded_transition_sensitivity_contrasts.csv"))

occupancy <- decoded %>%
  group_by(subject_ID, subject_ID_character, group, type, age, gender) %>%
  summarise(
    cooperative_occupancy = mean(cooperative_probability),
    .groups = "drop"
  )

fit_occupancy <- function(data, sensitivity, adjusted) {
  analysis_data <- data
  if (startsWith(sensitivity, "online_")) {
    analysis_data <- filter(analysis_data, type == "WEB")
  }
  if (adjusted || sensitivity == "online_complete_case") {
    analysis_data <- filter(analysis_data, !is.na(age), !is.na(gender))
  }
  analysis_data <- analysis_data %>% mutate(
    diag = factor(group, levels = diagnoses),
    age_sc = (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE)
  )
  formula <- if (adjusted) {
    cooperative_occupancy ~ diag + age_sc + gender
  } else {
    cooperative_occupancy ~ diag
  }
  fitted <- lm(formula, data = analysis_data)
  pair_table <- as.data.frame(summary(
    contrast(
      emmeans(fitted, ~ diag, weights = "proportional"),
      "pairwise",
      adjust = "none"
    ),
    infer = c(TRUE, TRUE),
    adjust = "none"
  ))
  as_tibble(pair_table) %>% transmute(
    sensitivity,
    contrast,
    difference = estimate,
    standard_error = SE,
    df,
    lower = lower.CL,
    upper = upper.CL,
    p_raw = p.value,
    p_holm = p.adjust(p_raw, method = "holm"),
    observations = nrow(analysis_data)
  )
}

occupancy_contrasts <- bind_rows(
  fit_occupancy(occupancy, "online_only", FALSE),
  fit_occupancy(occupancy, "online_complete_case", FALSE),
  fit_occupancy(occupancy, "online_age_gender_adjusted", TRUE),
  fit_occupancy(occupancy, "all_modes_age_gender_adjusted", TRUE)
)
atomic_write_csv(
  occupancy_contrasts,
  file.path(output_dir, "probability_weighted_cooperative_occupancy_sensitivity.csv")
)

run_provenance <- list(
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  R_version = R.version.string,
  platform = R.version$platform,
  analysis_script_sha256 = sha256_file(
    file.path("scripts", "run_trustee_covariate_sensitivity.R")
  ),
  sourced_robustness_script_sha256 = sha256_file(
    file.path("scripts", "run_trustee_state_count_robustness.R")
  ),
  helper_sha256 = sha256_file(file.path("scripts", "hmm_cache_utils.R")),
  behavioral_data_sha256 = sha256_file(file.path("Data", "full_RTG_data.csv")),
  demographics_data_sha256 = sha256_file(file.path("Data", "demographics.csv")),
  validated_model_artifact_path = normalizePath(model_rds),
  validated_model_artifact_sha256 = sha256_file(model_rds),
  technical_checks_path = normalizePath(technical_checks_csv),
  technical_checks_sha256 = sha256_file(technical_checks_csv),
  selected_state_count = artifact$selected_state_count,
  selected_fit_label = artifact$selected_fit_label,
  technical_checks_passed = all(technical_checks$passed)
)
atomic_write_csv(
  enframe(unlist(run_provenance), name = "field", value = "value"),
  file.path(output_dir, "RUN_MANIFEST.csv")
)

atomic_save_rds(
  list(
    model_source = normalizePath(model_rds),
    state_means = means,
    status = status,
    transition_contrasts = contrasts,
    occupancy_contrasts = occupancy_contrasts,
    provenance = run_provenance,
    generated_at = Sys.time()
  ),
  file.path(output_dir, "covariate_sensitivity_result.rds")
)
message("Online and demographic sensitivity checks complete")
