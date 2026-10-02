# Reproduction orchestration. Numerical fits and tests are defined in the
# accompanying manuscript and function modules.
with_treatment_contrasts <- function(operation) {
  previous <- getOption("contrasts")
  on.exit(options(contrasts = previous), add = TRUE)
  options(contrasts = c("contr.treatment", "contr.poly"))
  operation()
}

read_analysis_chunks <- function(path = "analysis_main.Rmd") {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^```\\{r", lines)
  lapply(seq_along(starts), function(i) {
    start <- starts[[i]]
    end <- which(seq_along(lines) > start & lines == "```")[[1]]
    header <- lines[[start]]
    label <- sub("^```\\{r\\s*", "", header)
    label <- trimws(sub("[,}].*$", "", label))
    if (!nzchar(label)) label <- paste0("unnamed_", i)
    list(label = label, header = header,
         code = paste(lines[seq.int(start + 1L, end - 1L)], collapse = "\n"))
  })
}

run_batches <- function(seeds, operation) {
  cores <- as.integer(Sys.getenv("RTG_WORKERS", "4"))
  if (is.na(cores) || cores < 1L) stop("RTG_WORKERS must be positive")
  cores <- min(cores, length(seeds))
  results <- parallel::mclapply(seeds, operation, mc.cores = cores,
                               mc.set.seed = FALSE, mc.preschedule = FALSE)
  if (any(vapply(results, inherits, logical(1), "try-error"))) {
    stop("A bootstrap worker failed; no pooled result was accepted")
  }
  results
}

pool_transition_batches <- function(parts, seeds, per_batch, summarise_draws) {
  result <- parts[[1]]
  draws <- dplyr::bind_rows(lapply(seq_along(parts), function(i) {
    stopifnot(parts[[i]]$bootstrap_status$successful_draws == per_batch,
              identical(parts[[i]]$point, result$point))
    dplyr::mutate(parts[[i]]$bootstrap$draws,
                  iteration = iteration + (i - 1L) * per_batch)
  }))
  total <- per_batch * length(parts)
  stopifnot(dplyr::n_distinct(draws$iteration) == total,
            !anyDuplicated(draws[c("iteration", "diagnosis", "investment", "metric")]))
  result$bootstrap <- summarise_draws(draws,
    dplyr::select(result$point, -dplyr::any_of("state_count")))
  result$bootstrap_status <- dplyr::bind_rows(lapply(parts, `[[`, "bootstrap_status")) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), sum))
  result$bootstrap_seeds <- as.integer(seeds)
  result
}

release_pool_primary <- function(model, data, bootstrap_draws, seed) {
  stopifnot(bootstrap_draws %% 8L == 0L)
  seeds <- as.integer(seed + 0:7)
  per_batch <- as.integer(bootstrap_draws / 8L)
  parts <- run_batches(seeds, function(s) {
    trustee_run_corrected_five_state_analysis(model, data, per_batch, s)
  })
  pool_transition_batches(parts, seeds, per_batch, trustee_summarise_bootstrap)
}

release_pool_state_count <- function(model, data, bootstrap_draws, seed) {
  stopifnot(bootstrap_draws %% 8L == 0L)
  seeds <- as.integer(seed + 1000L * (0:7))
  per_batch <- as.integer(bootstrap_draws / 8L)
  parts <- run_batches(seeds, function(s) {
    trustee_run_state_count_sensitivity(model, data, per_batch, s)
  })
  pool_transition_batches(parts, seeds, per_batch, trustee_summarise_sensitivity_draws)
}

compute_demographic_cooperation <- function(model, B) {
  stopifnot(B %% 8L == 0L)
  prepared <- prepare_demographics(model)
  scenarios <- c("all_modes_age_gender_adjusted", "online_age_gender_adjusted",
                 "online_only", "online_complete_case")
  outputs <- list()
  for (i in seq_along(scenarios)) {
    name <- scenarios[[i]]
    config <- prepared$configs[[name]]
    seeds <- as.integer(21011L + 10000L * (i - 1L) + 1000L * (0:7))
    parts <- run_batches(seeds, function(s) sample_demographics(config, s, B / 8L))
    draws <- dplyr::bind_rows(lapply(seq_along(parts), function(j) {
      stopifnot(parts[[j]]$status$successful == B / 8L)
      dplyr::mutate(parts[[j]]$draws, iteration = iteration + (j - 1L) * B / 8L)
    }))
    stopifnot(dplyr::n_distinct(draws$iteration) == B,
              !anyDuplicated(draws[c("iteration", "diagnosis", "investment", "metric")]))
    outputs[[name]] <- list(sensitivity = name, point = config$point,
      bootstrap = trustee_summarise_bootstrap(draws, config$point),
      status = dplyr::bind_rows(lapply(parts, `[[`, "status")) |>
        dplyr::summarise(dplyr::across(dplyr::everything(), sum)) |>
        dplyr::mutate(sensitivity = name),
      seeds = seeds, gender_weights = config$gender_weights,
      rejected = lapply(parts, `[[`, "rejected"),
      roundoff = lapply(parts, `[[`, "roundoff"))
  }
  list(results = outputs, counts = prepared$counts)
}

validate_saved_search_inputs <- function(result, templates, response_values, order_function = NULL) {
  same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 0, check.attributes = FALSE))
  for (key in names(result$selected_models)) {
    k <- as.integer(key)
    families <- intersect(c("simple", "group", "ordered_simple", "ordered_group"),
                          names(result$selected_models[[key]]))
    stopifnot(all(c("simple", "group") %in% families))
    for (family in families) {
      model <- result$selected_models[[key]][[family]]
      template <- templates[[sub("^ordered_", "", family)]][[k]]
      stopifnot(same(model@ntimes, template@ntimes),
                same(model@prior@x, template@prior@x))
      for (state in seq_len(k)) {
        stopifnot(same(model@response[[state]][[1]]@y,
                       template@response[[state]][[1]]@y),
                  same(as.numeric(model@response[[state]][[1]]@y),
                       as.numeric(response_values)),
                  same(model@transition[[state]]@x, template@transition[[state]]@x))
        if ("yield" %in% methods::slotNames(model@response[[state]][[1]])) {
          stopifnot(same(model@response[[state]][[1]]@yield,
                        template@response[[state]][[1]]@yield))
        }
      }
    }
    if (!is.null(order_function)) {
      stopifnot(!is.null(result$selected_models[[key]]$ordered_group),
        same(depmixS4::getpars(result$selected_models[[key]]$ordered_group),
             depmixS4::getpars(order_function(result$selected_models[[key]]$group))))
      if ("ordered_simple" %in% families) {
        stopifnot(same(depmixS4::getpars(result$selected_models[[key]]$ordered_simple),
                       depmixS4::getpars(order_function(result$selected_models[[key]]$simple))))
      }
    }
  }
  invisible(TRUE)
}
