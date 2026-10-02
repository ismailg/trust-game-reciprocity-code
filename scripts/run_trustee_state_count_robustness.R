#!/usr/bin/env Rscript

# Standalone robustness analysis for the trustee hidden Markov model (HMM).
#
# This script intentionally does not read, render, or modify manuscript files.
# It applies the same transition-summary logic used for the manuscript's
# six-state diagnostic-contrast model to five-, six-, and seven-state
# diagnostic-contrast models. Results are written under results/HMM/ only.

suppressPackageStartupMessages({
  library(tidyverse)
  library(depmixS4)
  library(lme4)
  library(emmeans)
  library(broom)
  library(nnet)
})

source(file.path("scripts", "hmm_cache_utils.R"))
register_hmm_response_classes()

parse_cli <- function(args) {
  get_value <- function(name, default) {
    hit <- grep(paste0("^", name, "="), args, value = TRUE)
    if (length(hit) == 0) return(default)
    sub(paste0("^", name, "="), "", hit[[length(hit)]])
  }
  list(
    states = as.integer(strsplit(get_value("--states", "5,6,7"), ",", fixed = TRUE)[[1]]),
    bootstrap = as.integer(get_value("--bootstrap", "0")),
    seed = as.integer(get_value("--seed", "11")),
    output_dir = get_value(
      "--output-dir",
      file.path("results", "HMM", "state_count_robustness")
    ),
    force = tolower(get_value("--force", "false")) %in% c("true", "1", "yes")
  )
}

config <- parse_cli(commandArgs(trailingOnly = TRUE))
if (length(config$states) == 0 || any(!config$states %in% 5:7)) {
  stop("--states must contain one or more of 5,6,7")
}
if (is.na(config$bootstrap) || config$bootstrap < 0) stop("--bootstrap must be >= 0")
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

message("Trustee state-count robustness analysis")
message("States: ", paste(config$states, collapse = ", "))
message("Bootstrap draws: ", config$bootstrap)
message("Output: ", normalizePath(config$output_dir, mustWork = FALSE))

prepared <- prepare_rtg_data(".")
trustee_data_base <- prepared$trustee_dat %>% arrange(subject_ID, round)

load(file.path("modsData", "fittedSimple.RData"))
load(file.path("modsData", "fittedCtrst.RData"))
load(file.path("modsData", "newCtrst_trustee.RData"))

diagnoses <- c("CC", "BPD", "AD")
investment_grid <- c(5, 10, 15)

choose_contrast_model <- function(k) {
  direct <- fittedCtrst[[k]]
  warm <- newCtrst[[k]]
  if (is.null(direct) || is.null(warm)) stop("Missing contrast-model cache for ", k, " states")
  if (as.numeric(logLik(direct)) > as.numeric(logLik(warm))) {
    list(model = direct, source = "direct")
  } else {
    list(model = warm, source = "warm-started")
  }
}

primary_state_groups <- function(state_means, cooperative_threshold = 0.40) {
  k <- length(state_means)
  coop <- which(state_means >= cooperative_threshold)
  if (length(coop) == 0) {
    stop("No states meet the cooperative mean-return threshold of ", cooperative_threshold)
  }
  list(
    exploitative = 1:2,
    cooperative = coop,
    intermediate = setdiff(seq_len(k), c(1:2, coop))
  )
}

expected_state_means <- function(model) {
  purrr::map_dbl(seq_len(nstates(model)), function(i) {
    pars <- getpars(model@response[[i]][[1]])
    response <- vtdgaus(
      seq(0, 1, length.out = 61),
      pstart = pars,
      yield = rep(60, 61)
    )
    sum(seq(0, 1, length.out = 61) * dens(response))
  })
}

get_nstates_compat <- function(fit) {
  value <- tryCatch(fit@nstates, error = function(e) NULL)
  if (!is.null(value)) return(value)
  stop("Cannot determine number of states for object of class ", paste(class(fit), collapse = "/"))
}

hmm_dens_bundle <- function(fit) {
  dens_mat <- fit@dens
  if (length(dim(dens_mat)) == 3) dens_mat <- dens_mat[, 1, ]
  list(
    dens = as.matrix(dens_mat),
    trDens = fit@trDens,
    init = if (is.null(fit@init)) NULL else as.matrix(fit@init)
  )
}

expected_transitions_xi <- function(fit, data_df, id_col = "subject_ID") {
  k <- get_nstates_compat(fit)
  ids <- data_df[[id_col]]
  id_levels <- unique(as.character(ids))
  seq_index <- split(seq_len(nrow(data_df)), factor(ids, levels = id_levels))
  output <- vector("list", length(seq_index))
  names(output) <- names(seq_index)

  # depmixS4 stores trDens and forwardbackward()$xi as [to_state, from_state],
  # not [from_state, to_state]. Use the package's own forward-backward result
  # and transpose each time slice before building the multinomial responses.
  fb <- depmixS4::forwardbackward(fit, return.all = TRUE)
  xi_to_from <- fb$xi
  if (!identical(dim(xi_to_from), c(nrow(data_df), k, k))) {
    stop(
      "forwardbackward xi shape does not match the analysis data: ",
      paste(dim(xi_to_from), collapse = "x"), " versus ",
      paste(c(nrow(data_df), k, k), collapse = "x")
    )
  }

  for (s in seq_along(seq_index)) {
    ind <- seq_index[[s]]
    sequence_length <- length(ind)
    if (sequence_length < 2) next
    transition_rows <- ind[seq_len(sequence_length - 1)]
    xi_from_to <- aperm(
      xi_to_from[transition_rows, , , drop = FALSE],
      c(1, 3, 2)
    )
    if (any(!is.finite(xi_from_to))) {
      stop("Non-finite expected transition probabilities within subject ", names(seq_index)[s])
    }
    from_marginal <- apply(xi_from_to, c(1, 2), sum)
    expected_from_marginal <- fb$gamma[transition_rows, , drop = FALSE]
    orientation_error <- max(abs(from_marginal - expected_from_marginal))
    if (!is.finite(orientation_error) || orientation_error > 1e-6) {
      stop(
        "Expected-transition orientation check failed within subject ",
        names(seq_index)[s], "; maximum from-state marginal error = ",
        orientation_error
      )
    }
    output[[s]] <- list(idx = ind, Xi = xi_from_to)
  }
  output
}

integerize_counts_const <- function(fractional, count_scale = 2000L, tol = 1e-12) {
  row_sums <- rowSums(fractional)
  keep <- row_sums > tol
  kept <- fractional[keep, , drop = FALSE]
  if (nrow(kept) == 0) return(list(Y_int = kept, keep_idx = integer(0)))
  scaled <- round(kept * count_scale)
  tiny <- which(rowSums(scaled) == 0 & row_sums[keep] > 0)
  for (idx in tiny) scaled[idx, which.max(kept[idx, ])] <- 1L
  storage.mode(scaled) <- "integer"
  list(Y_int = scaled, keep_idx = which(keep))
}

build_multinom_frames <- function(fit, data_df, xi_list, count_scale = 2000L) {
  k <- get_nstates_compat(fit)
  counts <- replicate(k, list(), simplify = FALSE)
  covariates <- replicate(k, list(), simplify = FALSE)

  for (entry in xi_list) {
    if (is.null(entry)) next
    ind <- entry$idx
    xi <- entry$Xi
    if (dim(xi)[1] < 1) next
    cov_block <- data_df[ind[seq_len(length(ind) - 1)], , drop = FALSE] %>%
      transmute(
        next_investment = .data$next_investment,
        BPD = as.integer(.data$group == "BPD"),
        AD = as.integer(.data$group == "AD")
      )
    for (i in seq_len(k)) {
      counts[[i]][[length(counts[[i]]) + 1L]] <- xi[, i, , drop = FALSE]
      covariates[[i]][[length(covariates[[i]]) + 1L]] <- cov_block
    }
  }

  purrr::map(seq_len(k), function(i) {
    if (length(counts[[i]]) == 0) return(NULL)
    fractional <- do.call(rbind, lapply(counts[[i]], matrix, ncol = k))
    predictors <- bind_rows(covariates[[i]])
    integerized <- integerize_counts_const(fractional, count_scale)
    if (nrow(integerized$Y_int) == 0) return(NULL)
    y <- integerized$Y_int
    predictors <- predictors[integerized$keep_idx, , drop = FALSE]
    totals <- rowSums(y)
    keep <- totals > 0
    list(
      Y_int = y[keep, , drop = FALSE],
      predictors = predictors[keep, , drop = FALSE],
      totals = totals[keep],
      K = k
    )
  })
}

ensure_level_rows <- function(data, classes) {
  totals <- data %>%
    group_by(to_state) %>%
    summarise(total = sum(weight), .groups = "drop") %>%
    mutate(to_state = as.character(to_state))
  needed <- union(
    setdiff(as.character(classes), totals$to_state),
    totals$to_state[totals$total <= 0]
  )
  if (length(needed) == 0) return(data)
  dummy <- purrr::map_dfr(needed, function(state) {
    data[1, ] %>% mutate(to_state = factor(state, levels = levels(data$to_state)), weight = 1e-6)
  })
  bind_rows(data, dummy)
}

multinom_block_is_valid <- function(block_fit) {
  if (is.null(block_fit) || is.null(block_fit$fit)) return(FALSE)
  if (is.null(block_fit$fit$convergence) || block_fit$fit$convergence != 0L) return(FALSE)
  coefficients <- coef(block_fit$fit)
  if (any(!is.finite(coefficients))) return(FALSE)
  probe <- tryCatch(
    predict(
      block_fit$fit,
      newdata = expand.grid(
        next_investment = investment_grid,
        BPD = c(0L, 1L, 0L),
        AD = c(0L, 0L, 1L)
      ),
      type = "probs"
    ),
    error = function(condition) condition
  )
  !inherits(probe, "error") && all(is.finite(probe))
}

fit_block_multinom <- function(block_data, decay = 1e-6) {
  if (is.null(block_data) || nrow(block_data$Y_int) == 0) return(NULL)
  y <- block_data$Y_int
  predictors <- block_data$predictors
  k <- block_data$K
  long <- tidyr::expand_grid(
    row_id = seq_len(nrow(y)),
    to_state = factor(seq_len(k), levels = seq_len(k))
  ) %>%
    mutate(
      weight = as.vector(t(y)),
      next_investment = predictors$next_investment[row_id],
      BPD = predictors$BPD[row_id],
      AD = predictors$AD[row_id]
    ) %>%
    ensure_level_rows(seq_len(k)) %>%
    filter(weight > 0)

  fit <- tryCatch(
    nnet::multinom(
      to_state ~ next_investment + BPD + AD + next_investment:BPD + next_investment:AD,
      data = long,
      weights = weight,
      trace = FALSE,
      Hess = FALSE,
      maxit = 200,
      decay = decay,
      MaxNWts = 20000
    ),
    error = function(condition) NULL
  )
  result <- list(
    fit = fit,
    row_predictors = predictors,
    row_totals = block_data$totals,
    K = k
  )
  if (!multinom_block_is_valid(result)) NULL else result
}

predict_row_from_fit <- function(block_fit, investment, diagnosis) {
  k <- block_fit$K
  new_data <- data.frame(
    next_investment = investment,
    BPD = as.integer(diagnosis == "BPD"),
    AD = as.integer(diagnosis == "AD")
  )
  raw <- predict(block_fit$fit, newdata = new_data, type = "probs")
  desired <- as.character(seq_len(k))
  output <- setNames(numeric(k), desired)
  if (is.null(dim(raw))) {
    labels <- names(raw)
    if (is.null(labels)) labels <- desired[seq_along(raw)]
    output[labels] <- as.numeric(raw)
  } else {
    labels <- colnames(raw)
    if (is.null(labels)) labels <- desired[seq_len(ncol(raw))]
    output[labels] <- as.numeric(raw[1, ])
  }
  total <- sum(output)
  if (!is.finite(total) || total <= 0) output[] <- 1 / k else output <- output / total
  output
}

simulate_refit_block <- function(block_fit, count_scale = 2000L, decay = 1e-6) {
  if (is.null(block_fit)) return(NULL)
  desired <- as.character(seq_len(block_fit$K))
  probability <- as.matrix(predict(
    block_fit$fit,
    newdata = block_fit$row_predictors,
    type = "probs"
  ))
  if (is.null(colnames(probability))) colnames(probability) <- desired[seq_len(ncol(probability))]
  missing <- setdiff(desired, colnames(probability))
  if (length(missing) > 0) {
    probability <- cbind(
      probability,
      matrix(0, nrow(probability), length(missing), dimnames = list(NULL, missing))
    )
  }
  probability <- probability[, desired, drop = FALSE]

  simulated <- t(vapply(seq_len(nrow(probability)), function(i) {
    inclusion_probability <- min(block_fit$row_totals[i] / count_scale, 1)
    size <- rbinom(1, 1, inclusion_probability)
    if (size == 0) numeric(length(desired)) else {
      as.numeric(rmultinom(1, 1, probability[i, ]))
    }
  }, numeric(length(desired))))

  long <- tidyr::expand_grid(
    row_id = seq_len(nrow(simulated)),
    to_state = factor(desired, levels = desired)
  ) %>%
    mutate(
      weight = as.vector(t(simulated)),
      next_investment = block_fit$row_predictors$next_investment[row_id],
      BPD = block_fit$row_predictors$BPD[row_id],
      AD = block_fit$row_predictors$AD[row_id]
    ) %>%
    ensure_level_rows(desired) %>%
    filter(weight > 0)

  refit <- tryCatch(
    nnet::multinom(
      to_state ~ next_investment + BPD + AD + next_investment:BPD + next_investment:AD,
      data = long,
      weights = weight,
      trace = FALSE,
      Hess = FALSE,
      maxit = 200,
      decay = decay,
      MaxNWts = 20000
    ),
    error = function(condition) NULL
  )
  result <- list(
    fit = refit,
    row_predictors = block_fit$row_predictors,
    row_totals = block_fit$row_totals,
    K = block_fit$K
  )
  if (!multinom_block_is_valid(result)) NULL else result
}

stationary_distribution <- function(transition_matrix) {
  transition_matrix[!is.finite(transition_matrix)] <- 0
  for (i in seq_len(nrow(transition_matrix))) {
    row_total <- sum(transition_matrix[i, ])
    if (row_total <= 0) {
      transition_matrix[i, ] <- 1 / ncol(transition_matrix)
    } else {
      transition_matrix[i, ] <- transition_matrix[i, ] / row_total
    }
  }
  eig <- eigen(t(transition_matrix))
  vector <- Re(eig$vectors[, which.min(Mod(eig$values - 1))])
  vector / sum(vector)
}

compose_transition_matrix <- function(blocks, investment, diagnosis) {
  k <- length(blocks)
  transition_matrix <- matrix(NA_real_, k, k)
  for (i in seq_len(k)) {
    if (!is.null(blocks[[i]])) {
      transition_matrix[i, ] <- predict_row_from_fit(blocks[[i]], investment, diagnosis)
    }
  }
  transition_matrix
}

direct_transition_matrix <- function(
  model,
  investment,
  diagnosis,
  data_df = trustee_data_base
) {
  k <- nstates(model)
  if (!diagnosis %in% c("CC", "BPD", "AD")) stop("Unknown diagnosis: ", diagnosis)
  output <- matrix(NA_real_, nrow = k, ncol = k)

  for (from_state in seq_len(k)) {
    transition_model <- model@transition[[from_state]]
    design_names <- colnames(transition_model@x)
    design <- setNames(rep(0, length(design_names)), design_names)
    design["(Intercept)"] <- 1
    design["next_investment"] <- investment

    diagnosis_rows <- which(as.character(data_df$group) == diagnosis)
    if (length(diagnosis_rows) == 0L || nrow(transition_model@x) != nrow(data_df)) {
      stop("Cannot recover the fitted diagnosis contrast coding")
    }
    for (main_term in c("BPD_ctrst1", "AD_ctrst1")) {
      fitted_values <- unique(transition_model@x[diagnosis_rows, main_term])
      if (length(fitted_values) != 1L) stop("Diagnosis coding is not constant for ", main_term)
      design[main_term] <- fitted_values
      forward_interaction <- paste0("next_investment:", main_term)
      reverse_interaction <- paste0(main_term, ":next_investment")
      if (forward_interaction %in% design_names) {
        design[forward_interaction] <- investment * fitted_values
      }
      if (reverse_interaction %in% design_names) {
        design[reverse_interaction] <- investment * fitted_values
      }
    }
    coefficients <- transition_model@parameters$coefficients
    coefficients <- coefficients[design_names, , drop = FALSE]
    linear_predictors <- as.numeric(design %*% coefficients)
    stabilized <- linear_predictors - max(linear_predictors)
    output[from_state, ] <- exp(stabilized) / sum(exp(stabilized))
  }

  row_error <- max(abs(rowSums(output) - 1))
  if (!is.finite(row_error) || row_error > 1e-10 || any(!is.finite(output))) {
    stop("Direct fitted-model transition matrix is invalid")
  }
  output
}

metric_values <- function(transition_matrix, cooperative, exploitative) {
  stationary <- stationary_distribution(transition_matrix)
  weighted_probability <- function(from, to) {
    denominator <- sum(stationary[from])
    numerator <- sum(
      stationary[from] * rowSums(transition_matrix[from, to, drop = FALSE])
    )
    if (denominator > 0) numerator / denominator else NA_real_
  }
  c(
    retention = weighted_probability(cooperative, cooperative),
    entry = weighted_probability(exploitative, cooperative),
    exploitative_exit = weighted_probability(cooperative, exploitative),
    stationary_cooperation = sum(stationary[cooperative])
  )
}

transition_metric_grid <- function(blocks, cooperative, exploitative) {
  purrr::map_dfr(diagnoses, function(diagnosis) {
    purrr::map_dfr(investment_grid, function(investment) {
      transition_matrix <- compose_transition_matrix(blocks, investment, diagnosis)
      tibble(
        diagnosis = diagnosis,
        investment = investment,
        metric = names(metric_values(transition_matrix, cooperative, exploitative)),
        value = as.numeric(metric_values(transition_matrix, cooperative, exploitative))
      )
    })
  })
}

direct_transition_metric_grid <- function(model, cooperative, lower_return) {
  purrr::map_dfr(diagnoses, function(diagnosis) {
    purrr::map_dfr(investment_grid, function(investment) {
      transition_matrix <- direct_transition_matrix(model, investment, diagnosis)
      tibble(
        diagnosis = diagnosis,
        investment = investment,
        metric = names(metric_values(transition_matrix, cooperative, lower_return)),
        value = as.numeric(metric_values(transition_matrix, cooperative, lower_return))
      )
    })
  })
}

summarise_bootstrap_draws <- function(draws, point) {
  levels <- draws %>%
    group_by(diagnosis, investment, metric) %>%
    summarise(
      bootstrap_mean = mean(value, na.rm = TRUE),
      lower = quantile(value, 0.025, na.rm = TRUE),
      upper = quantile(value, 0.975, na.rm = TRUE),
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
      bootstrap_mean = mean(difference, na.rm = TRUE),
      lower = quantile(difference, 0.025, na.rm = TRUE),
      upper = quantile(difference, 0.975, na.rm = TRUE),
      p_two_sided = min(
        1,
        2 * min(mean(difference <= 0, na.rm = TRUE), mean(difference >= 0, na.rm = TRUE))
      ),
      .groups = "drop"
    ) %>%
    mutate(
      correction_family = if_else(
        metric %in% c("retention", "entry"),
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

run_transition_bootstrap <- function(blocks, cooperative, exploitative, point, draws, seed,
                                     progress_offset = 0L, progress_total = draws) {
  if (draws <= 0) return(NULL)
  set.seed(seed)
  storage <- vector("list", draws)
  progress_every <- max(1L, floor(draws / 20L))
  for (b in seq_len(draws)) {
    simulated_blocks <- lapply(blocks, simulate_refit_block)
    storage[[b]] <- transition_metric_grid(simulated_blocks, cooperative, exploitative) %>%
      mutate(iteration = b)
    if (b %% progress_every == 0 || b == draws) {
      message("  bootstrap ", b, "/", draws, " (", round(100 * b / draws), "%)")
      message(
        "  overall bootstrap ", progress_offset + b, "/", progress_total,
        " (", round(100 * (progress_offset + b) / progress_total), "%)"
      )
    }
  }
  summarise_bootstrap_draws(bind_rows(storage), point)
}

extract_model_comparison <- function(k, selected, source_name) {
  simple <- fittedSimple[[k]]
  ll_simple <- as.numeric(logLik(simple))
  ll_contrast <- as.numeric(logLik(selected))
  df_simple <- attr(logLik(simple), "df")
  df_contrast <- attr(logLik(selected), "df")
  statistic <- 2 * (ll_contrast - ll_simple)
  df_difference <- df_contrast - df_simple
  tibble(
    state_count = k,
    contrast_cache = source_name,
    simple_log_likelihood = ll_simple,
    contrast_log_likelihood = ll_contrast,
    simple_AIC = AIC(simple),
    contrast_AIC = AIC(selected),
    simple_BIC = BIC(simple),
    contrast_BIC = BIC(selected),
    simple_parameters = df_simple,
    contrast_parameters = df_contrast,
    likelihood_ratio_chisq = statistic,
    likelihood_ratio_df = df_difference,
    likelihood_ratio_p = pchisq(statistic, df_difference, lower.tail = FALSE)
  )
}

categorize_states <- function(state, groups) {
  case_when(
    state %in% groups$exploitative ~ "exploitative",
    state %in% groups$cooperative ~ "cooperative",
    TRUE ~ "intermediate"
  )
}

decoded_occupancy <- function(data, groups) {
  categories <- c("exploitative", "intermediate", "cooperative")
  subject <- data %>%
    mutate(category = factor(categorize_states(posterior_state, groups), levels = categories)) %>%
    count(subject_ID, group, category, name = "trials") %>%
    group_by(subject_ID, group) %>%
    complete(category = factor(categories, levels = categories), fill = list(trials = 0)) %>%
    mutate(proportion = trials / sum(trials)) %>%
    ungroup()

  summary <- subject %>%
    group_by(group, category) %>%
    summarise(
      subjects = n_distinct(subject_ID),
      mean_proportion = mean(proportion),
      sd_proportion = sd(proportion),
      median_proportion = median(proportion),
      .groups = "drop"
    )

  tests <- purrr::map_dfr(categories, function(category_name) {
    current <- subject %>% filter(category == category_name)
    fit <- aov(proportion ~ group, data = current)
    row <- broom::tidy(fit) %>% filter(term == "group")
    tibble(
      category = category_name,
      numerator_df = row$df,
      denominator_df = broom::tidy(fit) %>% filter(term == "Residuals") %>% pull(df),
      F = row$statistic,
      p = row$p.value
    )
  }) %>% mutate(p_holm_across_categories = p.adjust(p, method = "holm"))

  pairs <- purrr::map_dfr(categories, function(category_name) {
    current <- subject %>% filter(category == category_name)
    fit <- aov(proportion ~ group, data = current)
    as.data.frame(emmeans::contrast(emmeans(fit, ~ group), "pairwise", adjust = "holm")) %>%
      as_tibble() %>%
      transmute(
        category = category_name,
        contrast,
        difference = estimate,
        SE,
        df,
        t_ratio = t.ratio,
        p_holm_within_category = p.value
      )
  })
  list(subject = subject, summary = summary, tests = tests, pairs = pairs)
}

decoded_per_state_occupancy <- function(data, k, bootstrap_iterations = 1000L,
                                        seed = 202501L) {
  subject <- data %>%
    count(subject_ID, group, posterior_state, name = "trials") %>%
    group_by(subject_ID, group) %>%
    complete(posterior_state = seq_len(k), fill = list(trials = 0)) %>%
    mutate(proportion = trials / sum(trials)) %>%
    ungroup()
  summary <- subject %>%
    group_by(group, posterior_state) %>%
    summarise(
      subjects = n_distinct(subject_ID),
      mean_proportion = mean(proportion),
      sd_proportion = sd(proportion),
      median_proportion = median(proportion),
      .groups = "drop"
    )
  tests <- purrr::map_dfr(seq_len(k), function(state) {
    current <- subject %>% filter(posterior_state == state)
    fit <- aov(proportion ~ group, data = current)
    row <- broom::tidy(fit) %>% filter(term == "group")
    tibble(
      state = state,
      numerator_df = row$df,
      denominator_df = broom::tidy(fit) %>% filter(term == "Residuals") %>% pull(df),
      F = row$statistic,
      p = row$p.value
    )
  }) %>% mutate(p_holm_across_states = p.adjust(p, method = "holm"))
  pairs <- purrr::map_dfr(seq_len(k), function(state) {
    current <- subject %>% filter(posterior_state == state)
    fit <- aov(proportion ~ group, data = current)
    as.data.frame(emmeans::contrast(emmeans(fit, ~ group), "pairwise", adjust = "holm")) %>%
      as_tibble() %>%
      transmute(
        state = state,
        contrast,
        difference = estimate,
        SE,
        df,
        t_ratio = t.ratio,
        p_holm_within_state = p.value
      )
  })

  set.seed(seed + k)
  bootstrap <- subject %>%
    group_by(group, posterior_state) %>%
    summarise(
      draws = list(replicate(
        bootstrap_iterations,
        mean(sample(proportion, replace = TRUE))
      )),
      .groups = "drop"
    ) %>%
    transmute(
      group,
      state = posterior_state,
      bootstrap_mean = purrr::map_dbl(draws, mean),
      lower = purrr::map_dbl(draws, ~ quantile(.x, 0.025)),
      upper = purrr::map_dbl(draws, ~ quantile(.x, 0.975)),
      bootstrap_iterations = bootstrap_iterations
    )
  list(subject = subject, summary = summary, tests = tests, pairs = pairs, bootstrap = bootstrap)
}

safe_glmm_analysis <- function(data, outcome, event, include_start_state = FALSE) {
  formula_text <- paste0(
    outcome,
    " ~ diag * inv_sc + gap_sc",
    if (include_start_state) " + start_state" else "",
    " + (1 | subject_ID)"
  )
  formula <- as.formula(formula_text)
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      glmer(
        formula,
        data = data,
        family = binomial,
        control = glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 1e6),
          check.conv.singular = "ignore"
        )
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      status = tibble(event = event, status = "failed", detail = conditionMessage(fit)),
      predictions = tibble(),
      contrasts = tibble(),
      coefficients = tibble()
    ))
  }

  at_values <- list(inv_sc = attr(data, "investment_scaled_grid"), gap_sc = 0)
  emm <- emmeans(fit, ~ diag | inv_sc, at = at_values, weights = "proportional")
  pred_raw <- as.data.frame(summary(emm, type = "response", infer = c(TRUE, TRUE)))
  probability_col <- intersect(c("prob", "response"), names(pred_raw))[[1]]
  lower_col <- intersect(c("asymp.LCL", "lower.CL"), names(pred_raw))[[1]]
  upper_col <- intersect(c("asymp.UCL", "upper.CL"), names(pred_raw))[[1]]
  predictions <- as_tibble(pred_raw) %>%
    transmute(
      event = event,
      diagnosis = as.character(diag),
      inv_sc,
      predicted_probability = .data[[probability_col]],
      lower = .data[[lower_col]],
      upper = .data[[upper_col]]
    )

  pair_raw <- as.data.frame(summary(
    contrast(emm, "pairwise", by = "inv_sc", adjust = "holm"),
    infer = c(TRUE, TRUE)
  ))
  pair_lower <- intersect(c("asymp.LCL", "lower.CL"), names(pair_raw))[[1]]
  pair_upper <- intersect(c("asymp.UCL", "upper.CL"), names(pair_raw))[[1]]
  contrasts <- as_tibble(pair_raw) %>%
    transmute(
      event = event,
      inv_sc,
      contrast,
      log_odds_difference = estimate,
      odds_ratio = exp(estimate),
      lower_odds_ratio = exp(.data[[pair_lower]]),
      upper_odds_ratio = exp(.data[[pair_upper]]),
      p_holm_within_investment = p.value
    )

  coefficient_raw <- as.data.frame(summary(fit)$coefficients)
  coefficients <- coefficient_raw %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    transmute(
      event = event,
      term,
      estimate = Estimate,
      standard_error = `Std. Error`,
      z = `z value`,
      p = `Pr(>|z|)`
    )

  convergence <- fit@optinfo$conv$lme4$messages
  details <- c(warnings, convergence)
  list(
    status = tibble(
      event = event,
      status = if (length(details) == 0) "ok" else "warning",
      detail = if (length(details) == 0) "" else paste(unique(details), collapse = " | "),
      observations = nrow(data),
      subjects = n_distinct(data$subject_ID),
      event_rate = mean(data[[outcome]]),
      singular = isSingular(fit)
    ),
    predictions = predictions,
    contrasts = contrasts,
    coefficients = coefficients
  )
}

decoded_transition_analyses <- function(data, groups) {
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

  investment_mean <- mean(transitions$next_investment_analysis, na.rm = TRUE)
  investment_sd <- sd(transitions$next_investment_analysis, na.rm = TRUE)
  if (!is.finite(investment_sd) || investment_sd == 0) investment_sd <- 1
  gap_mean <- mean(transitions$current_gap, na.rm = TRUE)
  gap_sd <- sd(transitions$current_gap, na.rm = TRUE)
  if (!is.finite(gap_sd) || gap_sd == 0) gap_sd <- 1
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

  retention_data <- transitions %>%
    filter(posterior_state %in% groups$cooperative) %>%
    mutate(stay_cooperative = as.integer(next_state %in% groups$cooperative)) %>%
    common_mutate()
  entry_data <- transitions %>%
    filter(posterior_state %in% groups$exploitative) %>%
    mutate(enter_cooperative = as.integer(next_state %in% groups$cooperative)) %>%
    common_mutate()
  exit_data <- transitions %>%
    filter(posterior_state %in% groups$cooperative) %>%
    mutate(
      exploitative_exit = as.integer(next_state %in% groups$exploitative),
      start_state = factor(posterior_state)
    ) %>%
    common_mutate()

  analyses <- list(
    safe_glmm_analysis(retention_data, "stay_cooperative", "retention", FALSE),
    safe_glmm_analysis(entry_data, "enter_cooperative", "entry", FALSE),
    safe_glmm_analysis(exit_data, "exploitative_exit", "exploitative_exit", TRUE)
  )
  list(
    status = bind_rows(lapply(analyses, `[[`, "status")),
    predictions = bind_rows(lapply(analyses, `[[`, "predictions")) %>%
      mutate(investment = investment_grid[match(inv_sc, scaled_grid)]),
    contrasts = bind_rows(lapply(analyses, `[[`, "contrasts")) %>%
      mutate(investment = investment_grid[match(inv_sc, scaled_grid)]),
    coefficients = bind_rows(lapply(analyses, `[[`, "coefficients")),
    support = transitions %>%
      filter(next_investment_analysis %in% investment_grid) %>%
      mutate(
        investment = next_investment_analysis,
        start_category = categorize_states(posterior_state, groups)
      ) %>%
      count(group, investment, start_category, name = "decoded_transitions")
  )
}

write_table <- function(data, path) {
  list_columns <- names(data)[vapply(data, is.list, logical(1))]
  if (length(list_columns) > 0) {
    stop(
      "Refusing to write list columns to ", path, ": ",
      paste(list_columns, collapse = ", ")
    )
  }
  readr::write_csv(data, path, na = "")
  invisible(path)
}

if (!identical(Sys.getenv("TRUSTEE_ROBUSTNESS_LIBRARY_ONLY"), "1")) {
all_results <- list()

for (k in config$states) {
  message("\nPreparing ", k, "-state diagnostic-contrast model")
  state_dir <- file.path(config$output_dir, paste0("states", k))
  dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(state_dir, sprintf("analysis_result_B%d.rds", config$bootstrap))
  selected <- choose_contrast_model(k)
  model <- order_mod_vtdgaus(selected$model)
  means <- expected_state_means(model)
  groups <- primary_state_groups(means, cooperative_threshold = 0.40)
  strict_groups <- list(
    exploitative = 1:2,
    cooperative = which(means >= 0.60),
    intermediate = setdiff(seq_len(k), c(1:2, which(means >= 0.60)))
  )
  if (length(strict_groups$cooperative) == 0) stop("No strict high-return states for ", k, " states")

  if (file.exists(result_path) && !config$force) {
    cached <- tryCatch(readRDS(result_path), error = function(e) NULL)
    cache_valid <- !is.null(cached) &&
      identical(cached$state_count, k) &&
      identical(cached$model_source, selected$source) &&
      identical(cached$transition_orientation, "depmix_xi_to_from_transposed_to_from_to") &&
      identical(cached$config$bootstrap, config$bootstrap) &&
      identical(cached$config$seed, config$seed) &&
      !is.null(cached$state_occupancy) &&
      identical(cached$state_groups$cooperative, groups$cooperative) &&
      identical(cached$state_groups$exploitative, groups$exploitative)
    if (isTRUE(cache_valid)) {
      message("Using validated existing result: ", result_path)
      all_results[[as.character(k)]] <- cached
      next
    }
    message("Ignoring stale or incompatible cached result: ", result_path)
  }

  state_dictionary <- tibble(
    state_count = k,
    state = seq_len(k),
    expected_return_fraction = means,
    primary_category = categorize_states(seq_len(k), groups),
    strict_category = categorize_states(seq_len(k), strict_groups)
  )
  model_comparison <- extract_model_comparison(k, model, selected$source)
  write_table(state_dictionary, file.path(state_dir, "state_dictionary.csv"))
  write_table(model_comparison, file.path(state_dir, "model_comparison.csv"))

  decoded_data <- trustee_data_base
  decoded_data$posterior_state <- posterior(model, type = "local")
  occupancy <- decoded_occupancy(decoded_data, groups)
  strict_occupancy <- decoded_occupancy(decoded_data, strict_groups)
  state_occupancy <- decoded_per_state_occupancy(decoded_data, k)
  decoded <- decoded_transition_analyses(decoded_data, groups)

  write_table(occupancy$subject, file.path(state_dir, "decoded_subject_category_occupancy.csv"))
  write_table(occupancy$summary, file.path(state_dir, "decoded_category_occupancy_summary.csv"))
  write_table(occupancy$tests, file.path(state_dir, "decoded_category_occupancy_tests.csv"))
  write_table(occupancy$pairs, file.path(state_dir, "decoded_category_occupancy_pairwise.csv"))
  write_table(strict_occupancy$summary, file.path(state_dir, "strict_high_return_occupancy_summary.csv"))
  write_table(strict_occupancy$tests, file.path(state_dir, "strict_high_return_occupancy_tests.csv"))
  write_table(strict_occupancy$pairs, file.path(state_dir, "strict_high_return_occupancy_pairwise.csv"))
  write_table(state_occupancy$subject, file.path(state_dir, "decoded_subject_state_occupancy.csv"))
  write_table(state_occupancy$summary, file.path(state_dir, "decoded_state_occupancy_summary.csv"))
  write_table(state_occupancy$tests, file.path(state_dir, "decoded_state_occupancy_tests.csv"))
  write_table(state_occupancy$pairs, file.path(state_dir, "decoded_state_occupancy_pairwise.csv"))
  write_table(state_occupancy$bootstrap, file.path(state_dir, "decoded_state_occupancy_bootstrap.csv"))
  write_table(decoded$status, file.path(state_dir, "decoded_glmm_status.csv"))
  write_table(decoded$predictions, file.path(state_dir, "decoded_glmm_predictions.csv"))
  write_table(decoded$contrasts, file.path(state_dir, "decoded_glmm_pairwise.csv"))
  write_table(decoded$coefficients, file.path(state_dir, "decoded_glmm_coefficients.csv"))
  write_table(decoded$support, file.path(state_dir, "decoded_transition_support.csv"))

  message("Fitting the transition-summary model for ", k, " states")
  xi <- expected_transitions_xi(model, trustee_data_base)
  frames <- build_multinom_frames(model, trustee_data_base, xi)
  blocks <- lapply(frames, fit_block_multinom)
  if (any(vapply(blocks, is.null, logical(1)))) stop("At least one transition block could not be fitted")

  point <- transition_metric_grid(blocks, groups$cooperative, groups$exploitative)
  point_strict <- transition_metric_grid(blocks, strict_groups$cooperative, strict_groups$exploitative)
  point <- point %>% mutate(
    state_count = k,
    classification = "cooperative_mean_return_at_least_0.40"
  )
  point_strict <- point_strict %>% mutate(state_count = k, classification = "strict_expected_return_at_least_0.60")
  write_table(bind_rows(point, point_strict), file.path(state_dir, "transition_metric_point_estimates.csv"))

  bootstrap <- run_transition_bootstrap(
    blocks,
    groups$cooperative,
    groups$exploitative,
    point %>% dplyr::select(diagnosis, investment, metric, value),
    config$bootstrap,
    config$seed,
    progress_offset = (match(k, config$states) - 1L) * config$bootstrap,
    progress_total = length(config$states) * config$bootstrap
  )
  if (!is.null(bootstrap)) {
    bootstrap$levels <- bootstrap$levels %>% mutate(state_count = k)
    bootstrap$contrasts <- bootstrap$contrasts %>% mutate(state_count = k)
    write_table(bootstrap$levels, file.path(state_dir, sprintf("bootstrap_levels_B%d.csv", config$bootstrap)))
    write_table(bootstrap$contrasts, file.path(state_dir, sprintf("bootstrap_contrasts_B%d.csv", config$bootstrap)))
  }

  result <- list(
    state_count = k,
    model_source = selected$source,
    model_comparison = model_comparison,
    state_dictionary = state_dictionary,
    state_groups = groups,
    strict_state_groups = strict_groups,
    occupancy = occupancy,
    strict_occupancy = strict_occupancy,
    state_occupancy = state_occupancy,
    decoded = decoded,
    point = point,
    point_strict = point_strict,
    bootstrap = bootstrap,
    transition_orientation = "depmix_xi_to_from_transposed_to_from_to",
    config = config,
    generated_at = Sys.time()
  )
  saveRDS(result, result_path)
  all_results[[as.character(k)]] <- result
  message("Saved ", k, "-state result")
}

combined_model <- bind_rows(lapply(all_results, `[[`, "model_comparison"))
combined_dictionary <- bind_rows(lapply(all_results, `[[`, "state_dictionary"))
combined_point <- bind_rows(lapply(all_results, function(x) bind_rows(x$point, x$point_strict)))
combined_occupancy <- bind_rows(c(
  unname(lapply(all_results, function(x) {
    x$occupancy$summary %>% mutate(
      state_count = x$state_count,
      classification = "cooperative_mean_return_at_least_0.40"
    )
  })),
  unname(lapply(all_results, function(x) {
    x$strict_occupancy$summary %>% mutate(
      state_count = x$state_count,
      classification = "strict_expected_return_at_least_0.60"
    )
  }))
))
combined_decoded_predictions <- bind_rows(lapply(all_results, function(x) {
  x$decoded$predictions %>% mutate(state_count = x$state_count)
}))
combined_decoded_contrasts <- bind_rows(lapply(all_results, function(x) {
  x$decoded$contrasts %>% mutate(state_count = x$state_count)
}))
combined_bootstrap_levels <- bind_rows(lapply(all_results, function(x) {
  if (is.null(x$bootstrap)) tibble() else x$bootstrap$levels
}))
combined_bootstrap_contrasts <- bind_rows(lapply(all_results, function(x) {
  if (is.null(x$bootstrap)) tibble() else x$bootstrap$contrasts
}))

write_table(combined_model, file.path(config$output_dir, "combined_model_comparison.csv"))
write_table(combined_dictionary, file.path(config$output_dir, "combined_state_dictionary.csv"))
write_table(combined_point, file.path(config$output_dir, "combined_transition_metric_point_estimates.csv"))
write_table(combined_occupancy, file.path(config$output_dir, "combined_decoded_occupancy.csv"))
write_table(combined_decoded_predictions, file.path(config$output_dir, "combined_decoded_glmm_predictions.csv"))
write_table(combined_decoded_contrasts, file.path(config$output_dir, "combined_decoded_glmm_pairwise.csv"))
if (nrow(combined_bootstrap_levels) > 0) {
  write_table(
    combined_bootstrap_levels,
    file.path(config$output_dir, sprintf("combined_bootstrap_levels_B%d.csv", config$bootstrap))
  )
  write_table(
    combined_bootstrap_contrasts,
    file.path(config$output_dir, sprintf("combined_bootstrap_contrasts_B%d.csv", config$bootstrap))
  )
}

manifest <- tibble(
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  requested_states = paste(config$states, collapse = ","),
  bootstrap_draws = config$bootstrap,
  seed = config$seed,
  transition_orientation = "depmix_xi_to_from_transposed_to_from_to",
  manuscript_files_modified = FALSE,
  note = "Standalone trustee HMM state-count robustness analysis"
)
write_table(manifest, file.path(config$output_dir, "run_manifest.csv"))

message("\nAnalysis complete")
}
