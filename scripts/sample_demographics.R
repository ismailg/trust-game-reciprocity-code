sample_demographics <- function(config, seed, B) {
 name <- config$sensitivity
set.seed(seed); started <- Sys.time()
storage <- vector('list', B); success <- 0L; attempts <- 0L; rejection_log <- list()
roundoff_log <- list()
while (success < B && attempts < max(B * 5L, B + 100L)) {
  attempts <- attempts + 1L
  blocks <- lapply(config$blocks, if (config$adjusted) simulate_augmented_block else trustee_simulate_transition_block)
  failed <- which(vapply(blocks, is.null, logical(1)))
  candidate <- if (length(failed)) simpleError(paste('Failed refit in states', paste(failed, collapse = ','))) else
    tryCatch(long_run_grid(blocks, config$adjusted, config$gender_weights), error = function(e) e)
  if (inherits(candidate, 'error')) {
    rejection_log[[length(rejection_log) + 1L]] <- tibble(attempt = attempts, detail = conditionMessage(candidate))
    next
  }
  stopifnot(nrow(candidate) == 9L, all(is.finite(candidate$value)),
            all(candidate$value >= -1e-12 & candidate$value <= 1 + 1e-12))
  outside <- candidate$value < 0 | candidate$value > 1
  if (any(outside)) {
    roundoff_log[[length(roundoff_log) + 1L]] <- candidate[outside, ] %>%
      mutate(attempt = attempts, value_before = value, value_after = pmin(1, pmax(0, value)))
    candidate$value <- pmin(1, pmax(0, candidate$value))
  }
  success <- success + 1L
  storage[[success]] <- mutate(candidate, iteration = success)
  if (success %% 250L == 0L || success == B) message(name, ': ', success, '/', B)
}
stopifnot(success == B)
draws <- bind_rows(storage)

 list(sensitivity = name, seed = seed, point = config$point, draws = draws,
 status = tibble(requested = B, successful = success, attempted = attempts, rejected = attempts - success),
 rejected = bind_rows(rejection_log), roundoff = bind_rows(roundoff_log))
}
