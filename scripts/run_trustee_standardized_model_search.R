#!/usr/bin/env Rscript

# Apply one fitting policy to every trustee diagnosis-transition HMM from two
# through seven states. This is an audit runner; the final reviewer-facing code
# is also placed in analysis_main.Rmd.

suppressPackageStartupMessages({
  library(depmixS4)
  library(digest)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

source(file.path("scripts", "hmm_cache_utils.R"))
register_hmm_response_classes()

parse_arguments <- function(arguments) {
  value <- function(name, default) {
    match <- grep(paste0("^", name, "="), arguments, value = TRUE)
    if (!length(match)) return(default)
    sub(paste0("^", name, "="), "", match[[length(match)]])
  }
  list(
    state_counts = as.integer(strsplit(
      value("--states", "2,3,4,5,6,7"), ",", fixed = TRUE
    )[[1]]),
    screen_starts = as.integer(value("--screen-starts", "20")),
    full_starts = as.integer(value("--full-starts", "3")),
    perturbations = as.integer(value("--perturbations", "5")),
    screen_iterations = as.integer(value("--screen-iterations", "20")),
    max_iterations = as.integer(value("--max-iterations", "2000")),
    base_seed = as.integer(value("--seed", "20250708")),
    output_dir = value(
      "--output-dir",
      file.path(
        "results", "HMM", "corrected_submission",
        "standardized_trustee_search_2026-07-21_v2"
      )
    )
  )
}

config <- parse_arguments(commandArgs(trailingOnly = TRUE))
if (anyNA(config$state_counts) || !length(config$state_counts) ||
    any(!config$state_counts %in% 2:7) || anyDuplicated(config$state_counts)) {
  stop("--states must be unique integers drawn from 2,3,4,5,6,7")
}
integer_settings <- unlist(config[c(
  "screen_starts", "full_starts", "perturbations", "screen_iterations",
  "max_iterations", "base_seed"
)])
if (anyNA(integer_settings) || any(integer_settings <= 0L) ||
    config$full_starts > config$screen_starts) {
  stop("Start counts, iteration counts, and seed must be positive and coherent")
}
if (dir.exists(config$output_dir) && length(list.files(
  config$output_dir, all.files = TRUE, no.. = TRUE
)) > 0L) {
  stop("Output directory must be new or empty: ", config$output_dir)
}

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$output_dir, "states"), showWarnings = FALSE)
options(contrasts = c("contr.sum", "contr.poly"))

script_path <- file.path("scripts", "run_trustee_standardized_model_search.R")
helper_path <- file.path("scripts", "hmm_cache_utils.R")
data_path <- file.path("Data", "full_RTG_data.csv")
sha256_file <- function(path) digest(file = path, algo = "sha256")
sha256_object <- function(object) digest(object, algo = "sha256", serialize = TRUE)
atomic_write_csv <- function(data, path) {
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  write_csv(data, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically write ", path)
  }
}

prepared <- prepare_rtg_data(".")
trustee_data <- prepared$trustee_dat %>% arrange(subject_ID, round)
initial_data <- prepared$priordat_trst %>% arrange(subject_ID, round)
if (nrow(trustee_data) != 8940L || nrow(initial_data) != 894L ||
    any(count(trustee_data, subject_ID)$n != 10L)) {
  stop("Expected 894 trustee sequences with ten rounds each")
}

manifest <- tibble(
  field = c(
    "started_utc", "states", "screen_starts", "full_starts",
    "perturbations_per_informed_solution", "screen_iterations",
    "max_iterations", "base_seed", "R_version", "platform",
    "status", "expected_fit_count", "depmixS4_version",
    "initialization_scheme", "screen_initialization", "validity_rule",
    "script_sha256", "helper_sha256", "data_sha256"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    paste(config$state_counts, collapse = ","),
    config$screen_starts, config$full_starts, config$perturbations,
    config$screen_iterations, config$max_iterations, config$base_seed,
    R.version.string, R.version$platform, "running",
    length(config$state_counts) * (
      2L * config$screen_starts + 2L * config$full_starts + 4L +
        2L * config$perturbations
    ),
    as.character(packageVersion("depmixS4")),
    "depmixS4_random_gamma_emission_refit",
    paste(
      "depmixS4 random.start=TRUE: random latent-state membership weights",
      "followed by emission refitting; prior and transition parameters retain",
      "their structured template values"
    ),
    "exact convergence; finite parameters; minimum occupancy 0.005; sigma > 1e-6",
    sha256_file(script_path), sha256_file(helper_path), sha256_file(data_path)
  )
)
atomic_write_csv(manifest, file.path(config$output_dir, "run_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "session_info.txt"))

fits_per_state <- 2L * config$screen_starts + 2L * config$full_starts +
  4L + 2L * config$perturbations
progress_total <- length(config$state_counts) * fits_per_state
progress_completed <- 0L
progress_path <- file.path(config$output_dir, "progress.json")
write_progress <- function(state_count, family, label, status) {
  progress_completed <<- progress_completed + 1L
  payload <- paste0(
    '{"current":', progress_completed,
    ',"total":', progress_total,
    ',"state_count":', state_count,
    ',"family":"', family,
    '","label":"', label,
    '","status":"', status,
    '","updated_utc":"', format(Sys.time(), tz = "UTC", usetz = TRUE), '"}'
  )
  writeLines(payload, progress_path)
  message("FIT_PROGRESS ", progress_completed, "/", progress_total,
          " K=", state_count, " ", family, " ", label, " ", status)
}
writeLines(
  paste0('{"current":0,"total":', progress_total,
         ',"state_count":null,"family":"pending","label":"pending",',
         '"status":"pending"}'),
  progress_path
)

is_converged_fit <- function(model) {
  if (is.null(model)) return(FALSE)
  paste(model@message, collapse = " ") %in% c(
    "Log likelihood converged to within tol. (relative change)",
    "Log likelihood converged to within tol. (absolute change)"
  ) && is.finite(as.numeric(logLik(model)))
}

build_structured_models <- function(state_count) {
  simple_template <- depmix(
    return_0 ~ 1,
    data = trustee_data,
    nstates = state_count,
    transition = ~ next_investment,
    prior = ~ investment,
    initdata = initial_data,
    family = gaussian(),
    ntimes = rep(10, nrow(initial_data))
  )
  group_template <- depmix(
    return_0 ~ 1,
    data = trustee_data,
    nstates = state_count,
    transition = ~ next_investment * (BPD_ctrst + AD_ctrst),
    prior = ~ investment,
    initdata = initial_data,
    family = gaussian(),
    ntimes = rep(10, nrow(initial_data))
  )
  means <- seq_len(state_count) / (state_count + 1)
  standard_deviations <- rep(1 / (state_count + 1), state_count)
  responses <- lapply(seq_len(state_count), function(state) {
    list(vtdgaus(
      y = trustee_data$return_0,
      yield = 3 * trustee_data$investment,
      pstart = c(means[[state]], standard_deviations[[state]])
    ))
  })
  rebuild <- function(template) {
    makeDepmix(
      response = responses,
      transition = template@transition,
      prior = template@prior,
      ntimes = rep(10, nrow(initial_data)),
      homogeneous = FALSE
    )
  }
  list(simple = rebuild(simple_template), group = rebuild(group_template))
}

warm_start_from_simple <- function(simple_fit, group_structure, state_count) {
  original <- getpars(simple_fit)
  prior_covariates <- 2L
  simple_transition_covariates <- 2L
  group_transition_covariates <- 6L
  new_parameters <- original[seq_len(state_count * prior_covariates)]
  for (from_state in seq_len(state_count)) {
    first <- state_count * prior_covariates + 1L +
      (from_state - 1L) * simple_transition_covariates * state_count
    last <- first + simple_transition_covariates * state_count - 1L
    simple_matrix <- matrix(
      original[first:last], ncol = simple_transition_covariates
    )
    expanded <- matrix(
      0, nrow = state_count, ncol = group_transition_covariates
    )
    expanded[, 1:2] <- simple_matrix
    new_parameters <- c(new_parameters, as.numeric(expanded))
  }
  response_first <- state_count * prior_covariates + 1L +
    state_count * simple_transition_covariates * state_count
  new_parameters <- c(new_parameters, original[response_first:length(original)])
  if (length(new_parameters) != length(getpars(group_structure))) {
    stop("Warm-start parameter expansion has the wrong length")
  }
  setpars(group_structure, new_parameters)
}

perturb_structured_start <- function(model, seed) {
  set.seed(seed)
  state_count <- nstates(model)
  prior_parameters <- getpars(model@prior)
  prior_fixed <- getpars(model@prior, which = "fixed")
  prior_parameters[!prior_fixed] <- prior_parameters[!prior_fixed] +
    rnorm(sum(!prior_fixed), sd = 0.02)
  complete_parameters <- prior_parameters
  for (state in seq_len(state_count)) {
    transition_parameters <- getpars(model@transition[[state]])
    transition_fixed <- getpars(model@transition[[state]], which = "fixed")
    transition_parameters[!transition_fixed] <-
      transition_parameters[!transition_fixed] +
      rnorm(sum(!transition_fixed), sd = 0.02)
    complete_parameters <- c(complete_parameters, transition_parameters)
  }
  for (state in seq_len(state_count)) {
    response_parameters <- getpars(model@response[[state]][[1]])
    response_parameters[[1]] <- response_parameters[[1]] + rnorm(1, sd = 0.005)
    response_parameters[[2]] <- response_parameters[[2]] * exp(rnorm(1, sd = 0.02))
    complete_parameters <- c(complete_parameters, response_parameters)
  }
  if (length(complete_parameters) != length(getpars(model))) {
    stop("Perturbed parameter vector has the wrong length")
  }
  output <- setpars(model, complete_parameters)
  cache_error <- max(
    abs(output@init - dens(output@prior)),
    max(vapply(seq_len(state_count), function(state) {
      max(abs(output@trDens[, , state] - dens(output@transition[[state]])))
    }, numeric(1))),
    max(vapply(seq_len(state_count), function(state) {
      max(abs(output@dens[, 1, state] - dens(output@response[[state]][[1]])))
    }, numeric(1)))
  )
  if (!is.finite(cache_error) || cache_error > 1e-12) {
    stop("Perturbed probability-cache mismatch: ", cache_error)
  }
  attr(output, "perturbation_cache_error") <- cache_error
  output
}

fit_candidate <- function(
    model, state_count, family, label, seed, max_iterations,
    random_start = FALSE, save_model = FALSE, state_directory = NULL,
    parent_label = NA_character_) {
  set.seed(seed)
  warnings_seen <- character()
  initial_parameter_sha256 <- sha256_object(getpars(model))
  started <- proc.time()[["elapsed"]]
  fitted <- tryCatch(
    withCallingHandlers(
      fit(
        model,
        emcontrol = em.control(
          maxit = max_iterations,
          tol = 1e-8,
          crit = "relative",
          random.start = random_start
        )
      ),
      warning = function(condition) {
        warnings_seen <<- c(warnings_seen, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) condition
  )
  elapsed <- proc.time()[["elapsed"]] - started
  error_message <- if (inherits(fitted, "error")) conditionMessage(fitted) else ""
  if (inherits(fitted, "error")) fitted <- NULL
  converged <- is_converged_fit(fitted)
  final_parameter_sha256 <- if (is.null(fitted)) {
    NA_character_
  } else {
    sha256_object(getpars(fitted))
  }
  initialization_method <- if (random_start) {
    "depmixS4_random_gamma_emission_refit"
  } else {
    "provided_parameter_start"
  }
  model_path <- NA_character_
  provenance_path <- NA_character_
  if (save_model && !is.null(fitted)) {
    model_path <- file.path(state_directory, "fits", paste0(label, ".rds"))
    provenance_path <- file.path(
      state_directory, "fits", paste0(label, ".provenance.rds")
    )
    saveRDS(fitted, model_path)
    saveRDS(
      list(
        state_count = state_count,
        family = family,
        label = label,
        parent_label = parent_label,
        seed = seed,
        max_iterations = max_iterations,
        random_start = random_start,
        initialization_method = initialization_method,
        classification = "soft",
        tolerance = 1e-8,
        criterion = "relative",
        initial_parameter_sha256 = initial_parameter_sha256,
        data_sha256 = sha256_file(data_path),
        script_sha256 = sha256_file(script_path),
        helper_sha256 = sha256_file(helper_path),
        final_parameter_sha256 = final_parameter_sha256,
        log_likelihood = as.numeric(logLik(fitted)),
        convergence_message = paste(fitted@message, collapse = " ")
      ),
      provenance_path
    )
  }
  status <- if (is.null(fitted)) "failed" else if (converged) "converged" else "nonconverged"
  write_progress(state_count, family, label, status)
  list(
    model = fitted,
    log = tibble(
      state_count = state_count,
      family = family,
      label = label,
      parent_label = parent_label,
      seed = seed,
      max_iterations = max_iterations,
      random_start = random_start,
      initialization_method = initialization_method,
      initial_parameter_sha256 = initial_parameter_sha256,
      final_parameter_sha256 = final_parameter_sha256,
      converged = converged,
      log_likelihood = if (is.null(fitted)) NA_real_ else as.numeric(logLik(fitted)),
      elapsed_seconds = elapsed,
      convergence_message = if (is.null(fitted)) "" else paste(fitted@message, collapse = " "),
      warnings = paste(unique(warnings_seen), collapse = " | "),
      error = error_message,
      model_path = model_path,
      provenance_path = provenance_path
    )
  )
}

run_membership_search <- function(
    template, state_count, family, seed_offset, state_directory) {
  screens <- vector("list", config$screen_starts)
  for (start_id in seq_len(config$screen_starts)) {
    screens[[start_id]] <- fit_candidate(
      template, state_count, family,
      paste0(family, "_membership_screen_", sprintf("%02d", start_id)),
      config$base_seed + state_count * 100000L + seed_offset + start_id,
      config$screen_iterations,
      random_start = TRUE
    )
  }
  successful <- which(vapply(screens, function(item) !is.null(item$model), logical(1)))
  if (length(successful) < config$full_starts) {
    stop("Too few successful random screens for K=", state_count, " ", family)
  }
  likelihoods <- vapply(
    screens[successful], function(item) as.numeric(logLik(item$model)), numeric(1)
  )
  chosen <- successful[order(likelihoods, decreasing = TRUE)[seq_len(config$full_starts)]]
  full <- vector("list", config$full_starts)
  for (rank in seq_len(config$full_starts)) {
    screen_id <- chosen[[rank]]
    full[[rank]] <- fit_candidate(
      screens[[screen_id]]$model,
      state_count,
      family,
      paste0(family, "_membership_full_", sprintf("%02d", rank),
             "_screen_", sprintf("%02d", screen_id)),
      config$base_seed + state_count * 100000L + seed_offset + 1000L + rank,
      config$max_iterations,
      random_start = FALSE,
      save_model = TRUE,
      state_directory = state_directory,
      parent_label = screens[[screen_id]]$log$label[[1]]
    )
  }
  list(
    candidates = full,
    logs = bind_rows(
      lapply(screens, `[[`, "log"),
      lapply(full, `[[`, "log")
    )
  )
}

candidate_diagnostics <- function(model) {
  if (!is_converged_fit(model)) {
    return(tibble(
      valid = FALSE, minimum_occupancy = NA_real_, minimum_sigma = NA_real_,
      maximum_absolute_transition_parameter = NA_real_, reason = "nonconverged"
    ))
  }
  ordered <- order_mod_vtdgaus(model)
  occupancy <- colMeans(forwardbackward(ordered, return.all = TRUE)$gamma)
  sigmas <- vapply(seq_len(nstates(ordered)), function(state) {
    unname(getpars(ordered@response[[state]][[1]])[[2]])
  }, numeric(1))
  transition_parameters <- unlist(lapply(ordered@transition, getpars))
  all_parameters <- getpars(ordered)
  finite <- all(is.finite(c(
    occupancy, sigmas, transition_parameters, all_parameters
  )))
  valid <- finite && min(occupancy) >= 0.005 && min(sigmas) > 1e-6
  tibble(
    valid = valid,
    minimum_occupancy = min(occupancy),
    minimum_sigma = min(sigmas),
    maximum_absolute_transition_parameter = max(abs(transition_parameters)),
    reason = if (valid) "valid" else if (!finite) "nonfinite" else "degenerate"
  )
}

select_best_candidate <- function(candidates, state_count, family, state_directory) {
  usable <- keep(candidates, function(item) !is.null(item$model))
  if (!length(usable)) stop("No fitted candidates for K=", state_count, " ", family)
  comparison <- map_dfr(usable, function(item) {
    bind_cols(
      item$log %>% select(
        state_count, family, label, parent_label, converged, log_likelihood,
        final_parameter_sha256, model_path, provenance_path
      ),
      candidate_diagnostics(item$model)
    )
  })
  valid_rows <- which(comparison$valid)
  if (!length(valid_rows)) stop("No valid candidates for K=", state_count, " ", family)
  best_row <- valid_rows[[which.max(comparison$log_likelihood[valid_rows])]]
  best_label <- comparison$label[[best_row]]
  best_item <- usable[[which(vapply(
    usable, function(item) identical(item$log$label[[1]], best_label), logical(1)
  ))]]
  best_model <- best_item$model
  best_raw_path <- file.path(state_directory, paste0("best_", family, "_model.rds"))
  best_ordered_path <- file.path(
    state_directory, paste0("best_", family, "_model_ordered.rds")
  )
  saveRDS(best_model, best_raw_path)
  best_ordered_model <- order_mod_vtdgaus(best_model)
  saveRDS(best_ordered_model, best_ordered_path)
  comparison <- comparison %>%
    mutate(selected = label == best_label) %>%
    arrange(desc(valid), desc(log_likelihood))
  write_csv(comparison, file.path(state_directory, paste0(family, "_candidates.csv")))
  best_provenance_path <- file.path(
    state_directory, paste0("best_", family, "_model.provenance.rds")
  )
  saveRDS(
    list(
      state_count = state_count,
      family = family,
      selected_label = best_label,
      log_likelihood = as.numeric(logLik(best_model)),
      raw_parameter_sha256 = sha256_object(getpars(best_model)),
      ordered_parameter_sha256 = sha256_object(getpars(best_ordered_model)),
      raw_file_sha256 = sha256_file(best_raw_path),
      ordered_file_sha256 = sha256_file(best_ordered_path),
      selected_candidate = comparison %>% filter(selected),
      script_sha256 = sha256_file(script_path),
      helper_sha256 = sha256_file(helper_path),
      data_sha256 = sha256_file(data_path)
    ),
    best_provenance_path
  )
  list(
    model = best_model,
    ordered_model = best_ordered_model,
    label = best_label,
    comparison = comparison,
    diagnostics = comparison %>% filter(selected) %>% select(-selected)
  )
}

run_state_search <- function(state_count) {
  message("Starting standardized trustee search for K=", state_count)
  state_directory <- file.path(config$output_dir, "states", paste0("states", state_count))
  dir.create(file.path(state_directory, "fits"), recursive = TRUE, showWarnings = FALSE)
  structures <- build_structured_models(state_count)

  simple_random <- run_membership_search(
    structures$simple, state_count, "simple", 10000L, state_directory
  )
  simple_direct <- fit_candidate(
    structures$simple, state_count, "simple", "simple_structured_direct",
    config$base_seed, config$max_iterations, save_model = TRUE,
    state_directory = state_directory, parent_label = "simple_structured_template"
  )
  simple_candidates <- c(simple_random$candidates, list(simple_direct))
  best_simple <- select_best_candidate(
    simple_candidates, state_count, "simple", state_directory
  )

  group_random <- run_membership_search(
    structures$group, state_count, "group", 50000L, state_directory
  )
  group_direct <- fit_candidate(
    structures$group, state_count, "group", "group_structured_direct",
    config$base_seed, config$max_iterations, save_model = TRUE,
    state_directory = state_directory, parent_label = "group_structured_template"
  )
  if (!is_converged_fit(simple_direct$model)) {
    stop("The structured no-diagnosis route did not converge for K=", state_count)
  }
  warm_structure <- warm_start_from_simple(
    simple_direct$model, structures$group, state_count
  )
  group_warm <- fit_candidate(
    warm_structure, state_count, "group", "group_structured_warm",
    config$base_seed, config$max_iterations, save_model = TRUE,
    state_directory = state_directory, parent_label = "simple_structured_direct"
  )
  best_simple_warm_structure <- warm_start_from_simple(
    best_simple$model, structures$group, state_count
  )
  group_best_simple_warm <- fit_candidate(
    best_simple_warm_structure, state_count, "group", "group_best_simple_warm",
    config$base_seed, config$max_iterations, save_model = TRUE,
    state_directory = state_directory, parent_label = best_simple$label
  )

  perturbation_candidates <- list()
  perturbation_logs <- list()
  perturbation_index <- 0L
  for (route in c("direct", "warm")) {
    base_model <- if (route == "direct") structures$group else warm_structure
    route_offset <- if (route == "direct") 1000L else 2000L
    for (perturbation in seq_len(config$perturbations)) {
      perturbation_index <- perturbation_index + 1L
      perturbation_seed <- config$base_seed + route_offset + perturbation
      label <- paste0(
        "group_perturbed_", route, "_", sprintf("%02d", perturbation)
      )
      perturbed <- tryCatch(
        perturb_structured_start(base_model, perturbation_seed),
        error = function(condition) condition
      )
      if (inherits(perturbed, "error")) {
        write_progress(state_count, "group", label, "invalid_start")
        failed <- list(
          model = NULL,
          log = tibble(
            state_count = state_count, family = "group", label = label,
            parent_label = paste0("group_structured_", route, "_start"),
            seed = perturbation_seed, max_iterations = config$max_iterations,
            random_start = FALSE,
            initialization_method = "provided_parameter_start",
            initial_parameter_sha256 = NA_character_,
            final_parameter_sha256 = NA_character_, converged = FALSE,
            log_likelihood = NA_real_,
            elapsed_seconds = 0, convergence_message = "", warnings = "",
            error = conditionMessage(perturbed), model_path = NA_character_,
            provenance_path = NA_character_
          )
        )
        perturbation_candidates[[perturbation_index]] <- failed
        perturbation_logs[[perturbation_index]] <- failed$log
      } else {
        fitted <- fit_candidate(
          perturbed, state_count, "group", label, perturbation_seed,
          config$max_iterations, save_model = TRUE,
          state_directory = state_directory,
          parent_label = paste0("group_structured_", route, "_start")
        )
        perturbation_candidates[[perturbation_index]] <- fitted
        perturbation_logs[[perturbation_index]] <- fitted$log
      }
    }
  }

  group_candidates <- c(
    group_random$candidates,
    list(group_direct, group_warm, group_best_simple_warm),
    perturbation_candidates
  )
  best_group <- select_best_candidate(
    group_candidates, state_count, "group", state_directory
  )
  all_logs <- bind_rows(
    simple_random$logs,
    simple_direct$log,
    group_random$logs,
    group_direct$log,
    group_warm$log,
    group_best_simple_warm$log,
    perturbation_logs
  )
  write_csv(all_logs, file.path(state_directory, "fit_log.csv"))

  state_profile <- map_dfr(seq_len(state_count), function(state) {
    ordered <- best_group$ordered_model
    parameters <- getpars(ordered@response[[state]][[1]])
    support <- seq(0, 1, length.out = 61)
    probabilities <- as.numeric(dens(vtdgaus(
      support, yield = rep(60, length(support)), pstart = parameters
    )))
    tibble(
      state_count = state_count,
      state = state,
      expected_return = sum(support * probabilities / sum(probabilities)),
      occupancy = mean(forwardbackward(ordered, return.all = TRUE)$gamma[, state]),
      sigma = unname(parameters[[2]])
    )
  })
  write_csv(state_profile, file.path(state_directory, "selected_state_profile.csv"))

  list(
    state_count = state_count,
    simple = best_simple,
    group = best_group,
    state_profile = state_profile
  )
}

results <- map(config$state_counts, run_state_search)
names(results) <- as.character(config$state_counts)

model_comparison <- map_dfr(results, function(result) {
  map_dfr(c("simple", "group"), function(family) {
    selected <- result[[family]]
    model <- selected$model
    likelihood <- as.numeric(logLik(model))
    parameters <- as.integer(attr(logLik(model), "df"))
    diagnostics <- candidate_diagnostics(model)
    tibble(
      state_count = result$state_count,
      family = family,
      selected_label = selected$label,
      log_likelihood = likelihood,
      parameters = parameters,
      AIC = -2 * likelihood + 2 * parameters,
      BIC_trials = -2 * likelihood + parameters * log(nrow(trustee_data)),
      BIC_sequences = -2 * likelihood + parameters * log(nrow(initial_data)),
      valid_candidates = sum(selected$comparison$valid),
      converged_candidates = sum(selected$comparison$converged),
      near_best_candidates = sum(
        selected$comparison$valid &
          selected$comparison$log_likelihood >= likelihood - 1
      ),
      minimum_occupancy = diagnostics$minimum_occupancy,
      minimum_sigma = diagnostics$minimum_sigma
    )
  })
})
write_csv(model_comparison, file.path(config$output_dir, "model_comparison.csv"))

nested_lrt <- map_dfr(results, function(result) {
  group_model <- result$group$model
  simple_model <- result$simple$model
  raw_statistic <- 2 * (
    as.numeric(logLik(group_model)) - as.numeric(logLik(simple_model))
  )
  if (raw_statistic < -1e-6) {
    stop(
      "Nested-model likelihood ordering failed for K=", result$state_count,
      ": group model is worse than the selected no-diagnosis model"
    )
  }
  statistic <- max(0, raw_statistic)
  degrees_of_freedom <-
    attr(logLik(group_model), "df") - attr(logLik(simple_model), "df")
  tibble(
    state_count = result$state_count,
    group_log_likelihood = as.numeric(logLik(group_model)),
    simple_log_likelihood = as.numeric(logLik(simple_model)),
    statistic = statistic,
    degrees_of_freedom = degrees_of_freedom,
    p_value = pchisq(statistic, degrees_of_freedom, lower.tail = FALSE)
  )
})
write_csv(nested_lrt, file.path(config$output_dir, "nested_lrt.csv"))
five_state_lrt <- nested_lrt %>% filter(state_count == 5L)
if (nrow(five_state_lrt)) {
  write_csv(five_state_lrt, file.path(config$output_dir, "five_state_lrt.csv"))
} else {
  five_state_lrt <- NULL
}

group_model_comparison <- model_comparison %>% filter(family == "group")
selection_summary <- tibble(
  criterion = c("AIC", "BIC_trials", "BIC_sequences"),
  selected_state_count = c(
    group_model_comparison$state_count[[which.min(group_model_comparison$AIC)]],
    group_model_comparison$state_count[[which.min(group_model_comparison$BIC_trials)]],
    group_model_comparison$state_count[[which.min(group_model_comparison$BIC_sequences)]]
  )
)
write_csv(selection_summary, file.path(config$output_dir, "selection_summary.csv"))

selected_model_provenance <- map_dfr(results, function(result) {
  map_dfr(c("simple", "group"), function(family) {
    selected <- result[[family]]
    selected_row <- selected$comparison %>% filter(selected)
    state_directory <- file.path(
      config$output_dir, "states", paste0("states", result$state_count)
    )
    raw_path <- file.path(state_directory, paste0("best_", family, "_model.rds"))
    ordered_path <- file.path(
      state_directory, paste0("best_", family, "_model_ordered.rds")
    )
    best_provenance_path <- file.path(
      state_directory, paste0("best_", family, "_model.provenance.rds")
    )
    tibble(
      state_count = result$state_count,
      family = family,
      selected_label = selected$label,
      parent_label = selected_row$parent_label,
      source_candidate_path = selected_row$model_path,
      source_candidate_provenance_path = selected_row$provenance_path,
      raw_model_path = raw_path,
      ordered_model_path = ordered_path,
      best_model_provenance_path = best_provenance_path,
      raw_file_sha256 = sha256_file(raw_path),
      ordered_file_sha256 = sha256_file(ordered_path),
      raw_parameter_sha256 = sha256_object(getpars(selected$model)),
      ordered_parameter_sha256 = sha256_object(getpars(selected$ordered_model)),
      log_likelihood = as.numeric(logLik(selected$model)),
      parameters = as.integer(attr(logLik(selected$model), "df")),
      valid = selected_row$valid,
      minimum_occupancy = selected_row$minimum_occupancy,
      minimum_sigma = selected_row$minimum_sigma
    )
  })
})
write_csv(
  selected_model_provenance,
  file.path(config$output_dir, "selected_model_provenance.csv")
)

saveRDS(
  list(
    config = config,
    model_comparison = model_comparison,
    selection_summary = selection_summary,
    nested_lrt = nested_lrt,
    five_state_lrt = five_state_lrt,
    selected_models = map(results, function(result) {
      list(
        simple = result$simple$model,
        group = result$group$model,
        ordered_group = result$group$ordered_model,
        simple_label = result$simple$label,
        group_label = result$group$label
      )
    })
  ),
  file.path(config$output_dir, "standardized_search_result.rds")
)

if (progress_completed != progress_total) {
  stop(
    "Fit-count mismatch: completed ", progress_completed,
    " of ", progress_total, " expected fits"
  )
}

source_hashes_at_completion <- c(
  script = sha256_file(script_path),
  helper = sha256_file(helper_path),
  data = sha256_file(data_path)
)
source_hashes_at_start <- c(
  script = manifest$value[manifest$field == "script_sha256"],
  helper = manifest$value[manifest$field == "helper_sha256"],
  data = manifest$value[manifest$field == "data_sha256"]
)
if (!identical(source_hashes_at_completion, source_hashes_at_start)) {
  stop("The script, helper, or governed data changed during the model search")
}

critical_artifacts <- list.files(
  config$output_dir, recursive = TRUE, full.names = TRUE
)
critical_artifacts <- critical_artifacts[
  !basename(critical_artifacts) %in% c(
    "progress.json", "run_manifest.csv", "artifact_hash_index.csv"
  )
]
artifact_hash_index <- tibble(
  relative_path = sub(
    paste0("^", config$output_dir, "/?"), "", critical_artifacts
  ),
  bytes = file.info(critical_artifacts)$size,
  sha256 = vapply(critical_artifacts, sha256_file, character(1))
) %>% arrange(relative_path)
write_csv(
  artifact_hash_index,
  file.path(config$output_dir, "artifact_hash_index.csv")
)

manifest <- bind_rows(
  manifest %>% filter(!field %in% c("status", "expected_fit_count")),
  tibble(
    field = c(
      "status", "expected_fit_count", "actual_fit_count", "completed_utc",
      "final_script_sha256", "final_helper_sha256", "final_data_sha256",
      "artifact_hash_index_sha256"
    ),
    value = c(
      "complete", progress_total, progress_completed,
      format(Sys.time(), tz = "UTC", usetz = TRUE),
      source_hashes_at_completion,
      sha256_file(file.path(config$output_dir, "artifact_hash_index.csv"))
    )
  )
)
atomic_write_csv(manifest, file.path(config$output_dir, "run_manifest.csv"))

writeLines(
  paste0('{"current":', progress_total, ',"total":', progress_total,
         ',"state_count":null,"family":"complete","label":"complete",',
         '"status":"complete","updated_utc":"',
         format(Sys.time(), tz = "UTC", usetz = TRUE), '"}'),
  progress_path
)
message("Standardized trustee model search complete: ", config$output_dir)
