# Keep the original multinomial fit and simulation exactly for unadjusted checks.
# The augmented versions add only age/gender predictors to the same likelihood.
fit_augmented_counts <- function(counts, predictors, row_totals) {
  stopifnot(nrow(counts) == nrow(predictors), ncol(counts) == 5L)
  long <- tidyr::expand_grid(row_id = seq_len(nrow(counts)),
                             to_state = factor(1:5, levels = 1:5))
  long <- bind_cols(long, predictors[long$row_id, , drop = FALSE]) %>%
    mutate(weight = as.vector(t(counts))) %>%
    trustee_ensure_destination_levels(1:5) %>% filter(weight > 0)
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(nnet::multinom(
    to_state ~ next_investment + BPD + AD + next_investment:BPD +
      next_investment:AD + age_sc + female,
    data = long, weights = weight, trace = FALSE, Hess = FALSE,
    maxit = 1000, decay = 1e-6, MaxNWts = 20000), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
    }), error = function(e) NULL)
  if (is.null(fit) || fit$convergence != 0L || any(!is.finite(coef(fit))) || length(warnings)) return(NULL)
  probe <- tryCatch(predict(fit, newdata = predictors, type = 'probs'), error = function(e) NULL)
  if (is.null(probe) || !identical(colnames(probe), as.character(1:5)) ||
      any(!is.finite(probe)) || max(abs(rowSums(probe) - 1)) > 1e-10) return(NULL)
  list(fit = fit, row_predictors = predictors, row_totals = row_totals,
       state_count = 5L, augmented = TRUE)
}

simulate_augmented_block <- function(block) {
  probabilities <- as.matrix(predict(block$fit, newdata = block$row_predictors, type = 'probs'))
  stopifnot(identical(colnames(probabilities), as.character(1:5)))
  simulated <- t(vapply(seq_len(nrow(probabilities)), function(row) {
    size <- rbinom(1, 1, min(block$row_totals[row] / 2000L, 1))
    if (!size) numeric(5L) else as.numeric(rmultinom(1, 1, probabilities[row, ]))
  }, numeric(5L)))
  fit_augmented_counts(simulated, block$row_predictors, block$row_totals)
}

strict_stationary <- function(P) {
  stopifnot(identical(dim(P), c(5L, 5L)), all(is.finite(P)), min(P) >= 0,
            max(abs(rowSums(P) - 1)) < 1e-10)
  # Independent linear solve checks the original eigensystem calculation.
  A <- t(P) - diag(5L); A[5, ] <- 1
  pi <- as.numeric(solve(A, c(0, 0, 0, 0, 1)))
  original <- trustee_stationary_distribution(P)
  stopifnot(all(is.finite(pi)), min(pi) > -1e-10,
            abs(sum(pi) - 1) < 1e-10,
            max(abs(pi %*% P - pi)) < 1e-9,
            max(abs(pi - original)) < 1e-8,
            sum(Mod(eigen(t(P), only.values = TRUE)$values - 1) < 1e-8) == 1L)
  original
}

long_run_grid <- function(blocks, adjusted, gender_weights = NULL) {
  if (adjusted) stopifnot(length(gender_weights) == 2L, abs(sum(gender_weights) - 1) < 1e-12)
  grid <- expand_grid(diagnosis = trustee_analysis_diagnoses,
                      investment = trustee_analysis_investments,
                      female = if (adjusted) c(0, 1) else 0) %>%
    mutate(next_investment = investment, BPD = as.integer(diagnosis == 'BPD'),
           AD = as.integer(diagnosis == 'AD'), age_sc = 0)
  rows <- lapply(blocks, function(b) as.matrix(predict(b$fit, newdata = grid, type = 'probs')))
  stopifnot(all(vapply(rows, function(x) identical(colnames(x), as.character(1:5)) && nrow(x) == nrow(grid), logical(1))))
  values <- vapply(seq_len(nrow(grid)), function(i) {
    P <- do.call(rbind, lapply(rows, function(x) x[i, ]))
    sum(strict_stationary(P)[3:5])
  }, numeric(1))
  grid %>% mutate(value = values,
                  averaging_weight = if (adjusted) gender_weights[female + 1L] else 1) %>%
    group_by(diagnosis, investment) %>% summarise(value = sum(value * averaging_weight), .groups = 'drop') %>%
    mutate(metric = 'stationary_cooperation') %>% dplyr::select(diagnosis, investment, metric, value)
}
