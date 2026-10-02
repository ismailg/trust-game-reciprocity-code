# Independent numerical checks with invented values, not participant records.
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
source("scripts/reproduction_helpers.R")
chunks <- read_analysis_chunks()
wanted <- c("trustee_bootstrap_correction_family", "trustee_make_pairwise_contrasts", "trustee_summarise_bootstrap",
            "trustee_summarise_sensitivity_draws")
for (chunk in chunks) for (expr in parse(text=chunk$code)) {
 if (is.call(expr) && identical(expr[[1]],as.name("<-")) &&
     is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) eval(expr)
}
set.seed(17)
metrics <- c("retention","repair_from_exploitative","exit_to_exploitative","stationary_cooperation")
draws <- expand_grid(iteration=1:100, diagnosis=c("CC","BPD","AD"),
                     investment=c(5,10,15), metric=metrics) %>%
 mutate(value=runif(n(),-.15,.15)+recode(diagnosis,CC=.4,BPD=.3,AD=.8))
point <- draws %>% group_by(diagnosis,investment,metric) %>%
 summarise(value=mean(value),.groups="drop")
check_summary <- function(fun, draws, point, expected_families) {
 out <- fun(draws,point)
 counts <- table(out$contrasts$correction_family)
 stopifnot(setequal(names(counts),names(expected_families)),
           identical(as.integer(counts[names(expected_families)]),as.integer(expected_families)))
 wide <- draws %>% pivot_wider(names_from=diagnosis,values_from=value)
 manual <- bind_rows(lapply(list(c("BPD","CC"),c("AD","CC"),c("AD","BPD")),function(pair){
  wide %>% mutate(delta=.data[[pair[1]]]-.data[[pair[2]]]) %>%
   group_by(investment,metric) %>% summarise(
    n=n(),left=sum(delta<=0),right=sum(delta>=0),
    lo=quantile(delta,.025),hi=quantile(delta,.975),.groups="drop") %>%
   mutate(contrast=paste(pair,collapse="-"),p=pmin(1,2*(pmin(left,right)+1)/(n+1)))
 }))
 joined <- inner_join(out$contrasts,manual,by=c("investment","metric","contrast"))
 stopifnot(nrow(joined)==nrow(out$contrasts),max(abs(joined$p_two_sided-joined$p))<1e-14,
           max(abs(joined$lower-joined$lo))<1e-14,max(abs(joined$upper-joined$hi))<1e-14)
 for (family in unique(joined$correction_family)) {
  x<-filter(joined,correction_family==family);o<-order(x$p)
  independent<-pmin(1,cummax((nrow(x):1)*x$p[o]))
  stopifnot(max(abs(independent-x$p_holm_family[o]))<1e-14)
 }
}
check_summary(trustee_summarise_sensitivity_draws,draws,point,
              c(retention_and_repair=18,exploitative_exit=9,stationary_cooperation=9))
primary_draws <- draws %>% mutate(metric=recode(metric,
 repair_from_exploitative="entry_from_clearly_low_return",
 exit_to_exploitative="exit_to_clearly_low_return"))
primary_point <- point %>% mutate(metric=recode(metric,
 repair_from_exploitative="entry_from_clearly_low_return",
 exit_to_exploitative="exit_to_clearly_low_return"))
primary_draws <- bind_rows(primary_draws,
 filter(primary_draws,metric=="retention") %>% mutate(metric="entry_from_near_one_third"),
 filter(primary_draws,metric=="retention") %>% mutate(metric="exit_to_near_one_third"))
primary_point <- bind_rows(primary_point,
 filter(primary_point,metric=="retention") %>% mutate(metric="entry_from_near_one_third"),
 filter(primary_point,metric=="retention") %>% mutate(metric="exit_to_near_one_third"))
check_summary(trustee_summarise_bootstrap,primary_draws,primary_point,
              c(retention_and_entry=27,exit_destination=18,stationary_cooperation=9))
check_summary(trustee_summarise_bootstrap,
 filter(primary_draws,metric=="stationary_cooperation"),
 filter(primary_point,metric=="stationary_cooperation"),c(stationary_cooperation=9))
message("Independent bootstrap-tail, interval and Holm-rule checks passed")

# Exercise pooling through the reporting boundary using invented batch draws.
# The supplement selects successful-draw counts by the scenario label.
scenario_names <- c("all_modes_age_gender_adjusted", "online_age_gender_adjusted",
                    "online_only", "online_complete_case")
demo_point <- filter(primary_point, metric == "stationary_cooperation")
prepare_demographics <- function(model) list(configs = setNames(lapply(scenario_names,
  function(name) list(sensitivity = name, point = demo_point, gender_weights = NULL)),
  scenario_names), counts = tibble())
run_batches <- function(seeds, worker) lapply(seeds, worker)
sample_demographics <- function(config, seed, B) list(
  draws = filter(primary_draws, metric == "stationary_cooperation", iteration <= B),
  status = tibble(requested = B, successful = B, attempted = B, rejected = 0L),
  rejected = tibble(), roundoff = tibble())
pooled <- compute_demographic_cooperation(NULL, 80L)
for (name in scenario_names) {
  status <- pooled$results[[name]]$status
  stopifnot(nrow(filter(status, sensitivity == name)) == 1L,
            status$requested == 80L, status$successful == 80L,
            status$attempted == 80L, status$rejected == 0L,
            n_distinct(pooled$results[[name]]$bootstrap$draws$iteration) == 80L)
}
message("Pooled scenario labels and reporting counts passed")
