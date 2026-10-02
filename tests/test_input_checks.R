# Synthetic two-state models: altered observations/design or ordered parameters
# must be rejected before fitted objects can be reused. No empirical data.
suppressPackageStartupMessages(library(depmixS4))
source("scripts/reproduction_helpers.R")
d <- data.frame(y=c(1, 2, 3, 4), x=c(0, 1, 0, 1))
m <- depmix(y~1, data=d, nstates=2, transition=~x, ntimes=c(2,2))
templates <- list(simple=list(NULL,m), group=list(NULL,m))
result <- list(selected_models=list(`2`=list(simple=m,group=m,ordered_group=m)))
validate_saved_search_inputs(result, templates, d$y, identity)
rejects <- function(x) inherits(try(x, silent=TRUE), "try-error")
bad <- result
bad$selected_models[["2"]]$group@response[[1]][[1]]@y[1,1] <- 99
stopifnot(rejects(validate_saved_search_inputs(bad, templates, d$y, identity)))
bad <- result
bad$selected_models[["2"]]$ordered_group@transition[[1]]@x[1,2] <- 99
stopifnot(rejects(validate_saved_search_inputs(bad, templates, d$y, identity)))
bad <- result
p <- getpars(m); p[length(p)] <- p[length(p)] + 1
bad$selected_models[["2"]]$ordered_group <- setpars(m,p)
stopifnot(rejects(validate_saved_search_inputs(bad, templates, d$y, identity)))
message("Synthetic input and ordered-model checks passed")
