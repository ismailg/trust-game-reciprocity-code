suppressPackageStartupMessages({
  library(tidyverse)
  library(depmixS4)
  library(MASS)
  library(numDeriv)
})

prepare_rtg_data <- function(base_dir = ".") {
  full_dat <- read.csv(file.path(base_dir, "Data", "full_RTG_data.csv")) %>%
    dplyr::filter(!is.na(group), group != "ASPD")

  problematic_subjects <- full_dat %>%
    count(subject_ID, round) %>%
    tidyr::pivot_wider(names_from = round, values_from = n, values_fill = 0) %>%
    dplyr::filter(dplyr::if_any(-subject_ID, ~ .x != 1)) %>%
    dplyr::pull(subject_ID)

  full_dat <- full_dat %>%
    dplyr::filter(!(subject_ID %in% problematic_subjects)) %>%
    arrange(subject_ID, round) %>%
    group_by(subject_ID) %>%
    mutate(game_id = cur_group_id()) %>%
    ungroup()

  clean_dat <- full_dat %>%
    arrange(subject_ID, round) %>%
    group_by(subject_ID) %>%
    mutate(
      return_0 = dplyr::if_else(is.na(return_ratio) & investment == 0, 0, return_ratio),
      fair_return = dplyr::if_else(investment >= 5, 2 * investment - 10, 0),
      better_than_fair = return > fair_return,
      worse_than_fair = return < fair_return,
      fair_gap = return - fair_return,
      lag_investment = lag(investment),
      lag_return = lag(return),
      lag_better_than_fair = lag(better_than_fair),
      lag_worse_than_fair = lag(worse_than_fair),
      investor_responsiveness_event = dplyr::case_when(
        role != "investor" ~ NA_real_,
        is.na(lag_better_than_fair) ~ NA_real_,
        !lag_better_than_fair ~ NA_real_,
        investment > lag_investment ~ 1,
        TRUE ~ 0
      ),
      investor_retaliation_event = dplyr::case_when(
        role != "investor" ~ NA_real_,
        is.na(lag_worse_than_fair) ~ NA_real_,
        !lag_worse_than_fair ~ NA_real_,
        investment < lag_investment ~ 1,
        TRUE ~ 0
      ),
      trustee_coaxing_event = dplyr::case_when(
        role != "trustee" ~ NA_real_,
        is.na(investment) | investment <= 0 ~ NA_real_,
        investment > 5 ~ NA_real_,
        is.na(return) ~ NA_real_,
        return >= investment ~ 1,
        TRUE ~ 0
      ),
      trustee_encouragement_event = dplyr::case_when(
        role != "trustee" ~ NA_real_,
        is.na(lag_investment) | is.na(lag_return) ~ NA_real_,
        investment >= lag_investment ~ NA_real_,
        return > lag_return ~ 1,
        TRUE ~ 0
      ),
      trustee_coop_pickup_event = dplyr::case_when(
        role != "trustee" ~ NA_real_,
        is.na(lag_investment) | is.na(lag_return) ~ NA_real_,
        investment <= lag_investment ~ NA_real_,
        return > lag_return ~ 1,
        TRUE ~ 0
      ),
      lag_fair_gap = lag(fair_gap),
      headroom_up = dplyr::if_else(role == "investor" & !is.na(lag_investment), as.numeric(lag_investment < 20), NA_real_),
      headroom_down = dplyr::if_else(role == "investor" & !is.na(lag_investment), as.numeric(lag_investment > 0), NA_real_)
    ) %>%
    ungroup() %>%
    group_by(game_id) %>%
    mutate(
      next_investment = lead(investment, default = 0),
      BPD_ctrst = factor(dplyr::if_else(group == "BPD", 1, 0)),
      AD_ctrst = factor(dplyr::if_else(group == "AD", 1, 0))
    ) %>%
    ungroup() %>%
    mutate(
      group = factor(group, levels = c("CC", "BPD", "AD")),
      subject_ID = factor(subject_ID)
    )

  inv_dat <- clean_dat %>% dplyr::filter(role == "investor")
  trustee_dat <- clean_dat %>% dplyr::filter(role == "trustee")

  list(
    full_dat = full_dat,
    clean_dat = clean_dat,
    inv_dat = inv_dat,
    trustee_dat = trustee_dat,
    priordat_inv = inv_dat %>% dplyr::filter(round == 1),
    priordat_trst = trustee_dat %>% dplyr::filter(round == 1)
  )
}

register_hmm_response_classes <- function() {
  if (!methods::isClass("discgaus")) {
    setClass("discgaus", contains = "response", slots = c(breaks = "numeric"))
  }
  if (!methods::isGeneric("discgaus")) {
    setGeneric("discgaus", function(y, pstart = NULL, fixed = NULL, ...) standardGeneric("discgaus"))
  }
  if (!methods::hasMethod("discgaus", signature(y = "ANY"))) {
    setMethod(
      "discgaus",
      signature(y = "ANY"),
      function(y, pstart = NULL, fixed = NULL, breaks = c(-Inf, seq(0, 19) + 0.5, Inf), ...) {
        y <- matrix(y, length(y))
        x <- matrix(1)
        parameters <- list()
        npar <- 2
        if (is.null(fixed)) fixed <- as.logical(rep(0, npar))
        if (!is.null(pstart)) {
          if (length(pstart) != npar) stop("length of 'pstart' must be ", npar)
          parameters$mu <- pstart[1]
          parameters$sigma <- pstart[2]
        } else {
          parameters <- list(mu = 10, sigma = 3)
        }
        new("discgaus", parameters = parameters, fixed = fixed, x = x, y = y, npar = npar, breaks = breaks)
      }
    )
  }
  if (!methods::hasMethod("show", "discgaus")) {
    setMethod("show", "discgaus", function(object) {
      cat("Gaussian with discrete support\n")
      cat("Parameters:\n")
      cat("mu:", object@parameters$mu, "\n")
      cat("sigma:", object@parameters$sigma, "\n")
    })
  }
  if (!methods::hasMethod("dens", "discgaus")) {
    setMethod("dens", "discgaus", function(object, log = FALSE) {
      p <- pnorm(object@breaks[-1], mean = object@parameters$mu, sd = object@parameters$sigma) -
        pnorm(object@breaks[-length(object@breaks)], mean = object@parameters$mu, sd = object@parameters$sigma)
      if (log) log(p[as.numeric(cut(object@y, breaks = object@breaks))]) else p[as.numeric(cut(object@y, breaks = object@breaks))]
    })
  }
  if (!methods::hasMethod("setpars", "discgaus")) {
    setMethod("setpars", "discgaus", function(object, values, which = "pars", ...) {
      npar <- npar(object)
      if (length(values) != npar) stop("length of 'values' must be", npar)
      nms <- names(object@parameters)
      switch(which,
        pars = {
          object@parameters$mu <- values[1]
          object@parameters$sigma <- values[2]
        },
        fixed = {
          object@fixed <- as.logical(values)
        }
      )
      names(object@parameters) <- nms
      object
    })
  }
  if (!methods::hasMethod("getpars", "discgaus")) {
    setMethod("getpars", "discgaus", function(object, which = "pars", ...) {
      switch(which,
        pars = unlist(object@parameters),
        fixed = object@fixed
      )
    })
  }
  if (!methods::hasMethod("fit", "discgaus")) {
    setMethod("fit", "discgaus", function(object, w) {
      if (missing(w)) w <- NULL
      negLL <- if (!is.null(w)) {
        function(pars) {
          object <- setpars(object, c(pars[1], exp(pars[2])))
          -sum(w * log(dens(object)))
        }
      } else {
        function(pars) {
          object <- setpars(object, c(pars[1], exp(pars[2])))
          -sum(log(dens(object)))
        }
      }
      pars <- optim(c(object@parameters$mu, log(object@parameters$sigma)), fn = negLL)$par
      setpars(object, c(pars[1], exp(pars[2])))
    })
  }

  if (!methods::isClass("truncdiscgaus")) {
    setClass("truncdiscgaus", contains = "discgaus", slots = c(min = "numeric", max = "numeric"))
  }
  # hasMethod() also returns TRUE for the inherited discgaus method. Require
  # an exact truncdiscgaus method so boundary mass is renormalized to 0--20.
  truncdiscgaus_dens_method <- methods::selectMethod(
    "dens", "truncdiscgaus"
  )
  if (!identical(
    as.character(truncdiscgaus_dens_method@defined[["object"]]),
    "truncdiscgaus"
  )) {
    setMethod("dens", "truncdiscgaus", function(object, log = FALSE) {
      breaks <- c(object@min, object@breaks[object@breaks > object@min & object@breaks < object@max], object@max)
      prec <- pnorm(object@max, mean = object@parameters$mu, sd = object@parameters$sigma) -
        pnorm(object@min, mean = object@parameters$mu, sd = object@parameters$sigma)
      if (prec < 1e-12) {
        p <- rep(1 / (length(breaks) - 1), length(object@y))
      } else {
        p <- pnorm(breaks[-1], mean = object@parameters$mu, sd = object@parameters$sigma) -
          pnorm(breaks[-length(breaks)], mean = object@parameters$mu, sd = object@parameters$sigma)
        p <- p / sum(p)
        p <- p[as.numeric(cut(object@y, breaks = object@breaks))]
      }
      if (log) log(p) else p
    })
  }
  stopifnot(identical(
    as.character(
      methods::selectMethod("dens", "truncdiscgaus")@defined[["object"]]
    ),
    "truncdiscgaus"
  ))
  if (!methods::isGeneric("truncdiscgaus")) {
    setGeneric("truncdiscgaus", function(y, pstart = NULL, fixed = NULL, ...) standardGeneric("truncdiscgaus"))
  }
  if (!methods::hasMethod("truncdiscgaus", signature(y = "ANY"))) {
    setMethod(
      "truncdiscgaus",
      signature(y = "ANY"),
      function(y, pstart = NULL, fixed = NULL, breaks = c(-Inf, seq(0, 19) + 0.5, Inf), min = -0.5, max = 20.5, ...) {
        y <- matrix(y, length(y))
        x <- matrix(1)
        parameters <- list()
        npar <- 2
        if (is.null(fixed)) fixed <- as.logical(rep(0, npar))
        if (!is.null(pstart)) {
          if (length(pstart) != npar) stop("length of 'pstart' must be ", npar)
          parameters$mu <- pstart[1]
          parameters$sigma <- pstart[2]
        }
        new("truncdiscgaus", parameters = parameters, fixed = fixed, x = x, y = y, npar = npar, breaks = breaks, min = min, max = max)
      }
    )
  }

  if (!methods::isClass("vtdgaus")) {
    setClass("vtdgaus", contains = "response", slots = c(yield = "numeric"))
  }
  if (!methods::isGeneric("vtdgaus")) {
    setGeneric("vtdgaus", function(y, pstart = NULL, fixed = NULL, ...) standardGeneric("vtdgaus"))
  }
  if (!methods::hasMethod("vtdgaus", signature(y = "ANY"))) {
    setMethod(
      "vtdgaus",
      signature(y = "ANY"),
      function(y, yield, pstart = NULL, fixed = NULL, ...) {
        y <- matrix(y, length(y))
        x <- matrix(1)
        parameters <- list()
        npar <- 2
        if (is.null(fixed)) fixed <- as.logical(rep(0, npar))
        if (!is.null(pstart)) {
          if (length(pstart) != npar) stop("length of 'pstart' must be ", npar)
          parameters$mu <- pstart[1]
          parameters$sigma <- pstart[2]
        } else {
          parameters <- list(mu = 0.5, sigma = 1)
        }
        new("vtdgaus", parameters = parameters, fixed = fixed, x = x, y = y, npar = npar, yield = yield)
      }
    )
  }
  if (!methods::hasMethod("show", "vtdgaus")) {
    setMethod("show", "vtdgaus", function(object) {
      cat("Gaussian with variable discrete support for percentage responses\n")
      cat("Parameters:\n")
      cat("mu:", object@parameters$mu, "\n")
      cat("sigma:", object@parameters$sigma, "\n")
    })
  }
  if (!methods::hasMethod("dens", "vtdgaus")) {
    setMethod("dens", "vtdgaus", function(object, log = FALSE) {
      prec <- pnorm(1 + 0.5 * (1 / 60), mean = object@parameters$mu, sd = object@parameters$sigma) -
        pnorm(0 - 0.5 * (1 / 60), mean = object@parameters$mu, sd = object@parameters$sigma)
      if (prec < 1e-12) {
        p <- 1 / (object@yield + 1)
      } else {
        p <- pnorm(object@y + 0.5 * (1 / object@yield), mean = object@parameters$mu, sd = object@parameters$sigma) -
          pnorm(object@y - 0.5 * (1 / object@yield), mean = object@parameters$mu, sd = object@parameters$sigma)
        norm <- pnorm(1 + 0.5 * (1 / object@yield), mean = object@parameters$mu, sd = object@parameters$sigma) -
          pnorm(0 - 0.5 * (1 / object@yield), mean = object@parameters$mu, sd = object@parameters$sigma)
        p <- p / norm
      }
      p[object@yield == 0] <- 1
      if (log) log(p) else p
    })
  }
  if (!methods::hasMethod("vcov", "vtdgaus")) {
    setMethod("vcov", signature(object = "vtdgaus"), function(object) {
      pars_now <- getpars(object)
      if (is.null(names(pars_now))) names(pars_now) <- c("mu", "sigma")
      u0 <- c(mu = pars_now[1], log_sigma = log(pars_now[2]))
      nll_u <- function(u) {
        mu <- u[1]
        sigma <- exp(u[2])
        obj <- setpars(object, c(mu, sigma))
        dens_vals <- dens(obj)
        -sum(log(pmax(dens_vals, 1e-12)))
      }
      H <- numDeriv::hessian(nll_u, u0)
      Vu <- tryCatch(solve(H), error = function(e) MASS::ginv(H))
      J <- matrix(
        c(1, 0, 0, exp(u0[2])),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("mu", "sigma"), c("mu", "log_sigma"))
      )
      Vp <- J %*% Vu %*% t(J)
      dimnames(Vp) <- list(c("mu", "sigma"), c("mu", "sigma"))
      Vp
    })
  }
  if (!methods::hasMethod("setpars", "vtdgaus")) {
    setMethod("setpars", "vtdgaus", function(object, values, which = "pars", ...) {
      npar <- npar(object)
      if (length(values) != npar) stop("length of 'values' must be", npar)
      nms <- names(object@parameters)
      switch(which,
        pars = {
          object@parameters$mu <- values[1]
          object@parameters$sigma <- values[2]
        },
        fixed = {
          object@fixed <- as.logical(values)
        }
      )
      names(object@parameters) <- nms
      object
    })
  }
  if (!methods::hasMethod("getpars", "vtdgaus")) {
    setMethod("getpars", "vtdgaus", function(object, which = "pars", ...) {
      switch(which,
        pars = unlist(object@parameters),
        fixed = object@fixed
      )
    })
  }
  if (!methods::hasMethod("fit", "vtdgaus")) {
    setMethod("fit", "vtdgaus", function(object, w) {
      if (missing(w)) w <- NULL
      negLL <- if (!is.null(w)) {
        function(pars) {
          object <- setpars(object, c(pars[1], exp(pars[2])))
          -sum(w * log(dens(object)))
        }
      } else {
        function(pars) {
          object <- setpars(object, c(pars[1], exp(pars[2])))
          -sum(log(dens(object)))
        }
      }
      pars <- optim(c(object@parameters$mu, log(object@parameters$sigma)), fn = negLL)$par
      setpars(object, c(pars[1], exp(pars[2])))
    })
  }
}

label_switch <- function(mod, labels) {
  # labels[i] is the new state label for original state i.
  if (!is(mod, "depmix") && !is(mod, "depmix.fitted")) {
    stop("label_switch() is only defined for depmix models.")
  }

  n_states <- mod@nstates
  if (
    length(labels) != n_states ||
      length(unique(labels)) != n_states ||
      !all(labels %in% seq_len(n_states))
  ) {
    stop("labels must be unique integers from 1 to ", n_states, ".")
  }

  inv_labels <- vapply(seq_len(n_states), function(x) which(labels == x), integer(1))
  tmp <- mod

  prior_pars <- getpars(mod@prior)
  prior_fixed <- getpars(mod@prior, which = "fixed")
  out_pars <- as.numeric(t(matrix(
    prior_pars,
    nrow = length(prior_pars) / n_states,
    byrow = TRUE
  )[, inv_labels]))
  out_fixed <- as.logical(t(matrix(
    prior_fixed,
    nrow = length(prior_fixed) / n_states,
    byrow = TRUE
  )[, inv_labels]))

  if (!tmp@prior@family$link == "identity") {
    tmp@prior@family$base <- labels[tmp@prior@family$base]
  }

  for (state in seq_len(n_states)) {
    trans_pars <- getpars(mod@transition[[inv_labels[state]]])
    trans_fixed <- getpars(mod@transition[[inv_labels[state]]], which = "fixed")
    out_pars <- c(out_pars, as.numeric(t(matrix(
      trans_pars,
      nrow = length(trans_pars) / n_states,
      byrow = TRUE
    )[, inv_labels])))
    out_fixed <- c(out_fixed, as.logical(t(matrix(
      trans_fixed,
      nrow = length(trans_fixed) / n_states,
      byrow = TRUE
    )[, inv_labels])))
    tmp@transition[[state]] <- mod@transition[[inv_labels[state]]]
    if (!tmp@transition[[state]]@family$link == "identity") {
      tmp@transition[[state]]@family$base <- labels[tmp@transition[[state]]@family$base]
    }
  }

  for (state in seq_len(n_states)) {
    out_pars <- c(out_pars, unlist(lapply(mod@response[[inv_labels[state]]], getpars)))
    out_fixed <- c(out_fixed, unlist(lapply(
      mod@response[[inv_labels[state]]],
      getpars,
      which = "fixed"
    )))
  }

  tmp <- setpars(tmp, out_fixed, which = "fixed")
  tmp <- setpars(tmp, out_pars)
  if (is(tmp, "depmix.fitted")) tmp@posterior <- viterbi(tmp)
  tmp
}

order_mod_truncdiscgaus <- function(mod) {
  ns <- nstates(mod)
  expected_values <- rep(0, ns)
  for (i in seq_len(ns)) {
    tpars <- getpars(mod@response[[i]][[1]])
    dmod <- truncdiscgaus(seq(0, 20), pstart = tpars, min = -0.5, max = 20.5)
    expected_values[i] <- sum(seq(0, 20) * dens(dmod))
  }
  label_switch(mod, rank(expected_values, ties.method = "first"))
}

order_mod_vtdgaus <- function(mod) {
  ns <- nstates(mod)
  expected_values <- rep(0, ns)
  for (i in seq_len(ns)) {
    tpars <- getpars(mod@response[[i]][[1]])
    dmod <- vtdgaus(seq(0, 1, length = 61), pstart = tpars, yield = rep(60, 61))
    expected_values[i] <- sum(seq(0, 1, length = 61) * dens(dmod))
  }
  label_switch(mod, rank(expected_values, ties.method = "first"))
}
