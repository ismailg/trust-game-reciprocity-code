expected_state_means <- function(model) {
    purrr::map_dbl(seq_len(nstates(model)), function(i) {
        pars <- getpars(model@response[[i]][[1]])
        response <- vtdgaus(seq(0, 1, length.out = 61), pstart = pars, yield = rep(60, 61))
        sum(seq(0, 1, length.out = 61) * dens(response))
    })
}
