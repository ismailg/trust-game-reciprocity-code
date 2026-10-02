compute_covariates <- function(model, trustee_data_base, output_dir) {
  # Check the fitted input arrays and economic state definitions before fitting
  # second-stage models. The HMM itself used sum coding for diagnostic group.
  same <- function(a, b) isTRUE(all.equal(a, b, tolerance=0, check.attributes=FALSE))
  means <- expected_state_means(model)
  initial <- dplyr::filter(trustee_data_base, round == 1L)
  design <- model.matrix(~ next_investment * (BPD_ctrst + AD_ctrst),
    data=trustee_data_base,
    contrasts.arg=list(BPD_ctrst="contr.sum", AD_ctrst="contr.sum"))
  stopifnot(depmixS4::nstates(model) == 5L, all(diff(means) > 0),
    means[1] <= 0.30, means[2] >= 0.30, means[2] < 0.40,
    means[2] - means[1] >= 0.05, identical(which(means >= 0.40), 3:5),
    same(model@ntimes, rep(10L, nrow(initial))),
    same(model@prior@x, model.matrix(~ investment, data=initial)))
  for (i in 1:5) stopifnot(
    same(as.numeric(model@response[[i]][[1]]@y), trustee_data_base$return_0),
    same(model@response[[i]][[1]]@yield, 3 * trustee_data_base$investment),
    same(model@transition[[i]]@x, design))
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  atomic_write_csv <- function(data, path) readr::write_csv(data, path, na = "")
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


  manifest <- tibble::tibble(field = c("technical_checks_passed", "selected_state_count", "behavioral_data_sha256", "demographics_data_sha256"),
    value = c("TRUE", "5", digest::digest(file="Data/full_RTG_data.csv", algo="sha256"), digest::digest(file="Data/demographics.csv", algo="sha256")))
  manifest <- dplyr::bind_rows(manifest, tibble::tibble(
    field = "validation_basis",
    value = "Recomputed input-alignment and five-state economic-classification checks"))
  readr::write_csv(manifest, file.path(output_dir,"RUN_MANIFEST.csv"))
  invisible(list(occupancy=occupancy_contrasts, transitions=contrasts, status=status))
}
