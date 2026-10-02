prepare_demographics <- function(model) {
raw <- prepare_rtg_data('.')$trustee_dat %>% arrange(subject_ID, round)
same <- function(x, y) isTRUE(all.equal(x, y, check.attributes = FALSE, tolerance = 0))
design <- model.matrix(~ next_investment * (BPD_ctrst + AD_ctrst), data = raw,
                      contrasts.arg = list(BPD_ctrst = 'contr.sum', AD_ctrst = 'contr.sum'))
initial <- raw %>% filter(round == 1)
stopifnot(nrow(raw) == 8940L, same(model@ntimes, rep(10, nrow(initial))),
          same(model@prior@x, model.matrix(~ investment, data = initial)))
for (state in 1:5) stopifnot(same(as.numeric(model@response[[state]][[1]]@y), raw$return_0),
                            same(model@response[[state]][[1]]@yield, 3 * raw$investment),
                            same(model@transition[[state]]@x, design))
demo <- read.csv('Data/demographics.csv', encoding = 'UTF-8', check.names = FALSE) %>%
  transmute(subject_ID = as.character(ID), age = suppressWarnings(as.numeric(Age)),
            female = recode(as.character(Gender), `1` = 0, `2` = 1, .default = NA_real_)) %>%
  mutate(age = if_else(!is.na(age) & age > 0 & age < 100, age, NA_real_))
conflicts <- demo %>% group_by(subject_ID) %>%
  summarise(ages = n_distinct(age, na.rm = TRUE), sexes = n_distinct(female, na.rm = TRUE), .groups = 'drop')
stopifnot(all(conflicts$ages <= 1L), all(conflicts$sexes <= 1L))
demo <- distinct(demo)
data <- raw %>% mutate(subject_ID_char = as.character(subject_ID)) %>%
  left_join(demo, by = c('subject_ID_char' = 'subject_ID'), relationship = 'many-to-one')
stopifnot(nrow(data) == nrow(raw))
expected <- trustee_expected_transitions(model, raw)
indices <- unlist(lapply(expected, function(e) head(e$indices, -1L)), use.names = FALSE)
stopifnot(length(indices) == 894L * 9L,
          all(raw$round[indices] %in% 1:9), all(is.finite(raw$next_investment[indices])))
age_mean <- mean(data$age[indices], na.rm = TRUE)
age_sd <- sd(data$age[indices], na.rm = TRUE)
covariates <- data[indices, ] %>% transmute(
  next_investment, BPD = as.integer(group == 'BPD'), AD = as.integer(group == 'AD'),
  age_sc = (age - age_mean) / age_sd, female,
  subject_ID = subject_ID_char, group = as.character(group), type)
full_frames <- lapply(1:5, function(state) {
  counts <- do.call(rbind, lapply(expected, function(e) matrix(e$xi[, state, , drop = FALSE], ncol = 5L)))
  integerized <- trustee_integerize_expected_counts(counts)
  list(counts = integerized$counts,
       predictors = covariates[integerized$kept_rows, ],
       totals = rowSums(integerized$counts), state_count = 5L)
})
original_frames <- trustee_build_transition_frames(model, raw, expected)
for (i in 1:5) {
  stopifnot(identical(full_frames[[i]]$counts, original_frames[[i]]$counts),
            identical(as.data.frame(dplyr::select(full_frames[[i]]$predictors, next_investment, BPD, AD)),
                      as.data.frame(original_frames[[i]]$predictors)))
}
spec_names <- c('all_reference', 'online_only', 'online_complete_case',
                'online_age_gender_adjusted', 'all_modes_age_gender_adjusted')
expected_n <- c(894L, 691L, 592L, 592L, 788L)
configs <- list(); counts <- list()
for (i in seq_along(spec_names)) {
  name <- spec_names[[i]]
  adjusted <- grepl('adjusted$', name)
  selection <- rep(TRUE, nrow(covariates))
  if (startsWith(name, 'online_')) selection <- selection & covariates$type == 'WEB'
  if (adjusted || name == 'online_complete_case') selection <- selection & !is.na(covariates$age_sc) & !is.na(covariates$female)
  subjects <- covariates %>% filter(selection) %>% distinct(subject_ID, group, type, female)
  stopifnot(nrow(subjects) == expected_n[[i]], !anyDuplicated(subjects$subject_ID))
  keep_ids <- subjects$subject_ID
  frames <- lapply(full_frames, function(f) {
    keep <- f$predictors$subject_ID %in% keep_ids
    cols <- c('next_investment', 'BPD', 'AD', if (adjusted) c('age_sc', 'female'))
    list(counts = f$counts[keep, , drop = FALSE], predictors = dplyr::select(f$predictors[keep, ], all_of(cols)),
         totals = f$totals[keep], state_count = 5L)
  })
  for (f in frames) {
    X <- model.matrix(as.formula(paste('~ next_investment * (BPD + AD)', if (adjusted) '+ age_sc + female' else '')), f$predictors)
    stopifnot(qr(X)$rank == ncol(X), all(is.finite(X)))
  }
  blocks <- lapply(frames, function(f) if (adjusted) fit_augmented_counts(f$counts, f$predictors, f$totals) else trustee_fit_transition_block(f))
  if (any(vapply(blocks, is.null, logical(1)))) stop('Point fit failed: ', name)
  weights <- if (adjusted) as.numeric(table(factor(subjects$female, levels = 0:1)) / nrow(subjects)) else NULL
  point <- long_run_grid(blocks, adjusted, weights)
  configs[[name]] <- list(sensitivity = name, adjusted = adjusted, frames = frames, blocks = blocks,
                         gender_weights = weights, subjects = subjects, point = point)
  counts[[name]] <- subjects %>% count(group, name = 'participants') %>% mutate(sensitivity = name)
  message('Prepared ', name, ': ', nrow(subjects), ' participants; all five fits converged')
}

 list(configs = configs, counts = bind_rows(counts))
}
