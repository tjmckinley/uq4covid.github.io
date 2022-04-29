## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(sitmo)
library(furrr)

## function to check counts (mainly useful for error checking)
## u: matrix of counts with columns:
##       t, S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
## cu: matrix of cumulative counts
## N: population size
checkCounts <- function(u, cu, N) {
    
    ## remove time column
    u <- u[, -1]
    cu <- cu[, -1]
    
    ## check states match population size in each age-class
    stopifnot(all((rowSums(u) - N) == 0))
    
    ## check states are positive
    stopifnot(all(u >= 0) & all(cu >= 0))
    
    ## check cumulative counts match up
    if(!all((cu[, 1] + cu[, 2] - N) == 0)) {
        browser()
    }
    stopifnot(all((cu[, 2] - rowSums(cu[, c(3, 5)])) >= 0))
    stopifnot(all((cu[, 3] - cu[, 4]) >= 0))
    stopifnot(all((cu[, 5] - cu[, 6]) >= 0))
    stopifnot(all((cu[, 6] - rowSums(cu[, c(7, 8, 10)])) >= 0))
    stopifnot(all((cu[, 8] - cu[, 9]) >= 0))
    stopifnot(all((cu[, 10] - rowSums(cu[, c(11, 12)])) >= 0))
    
    # ## check infective states and new incidence are valid
    ## KEEPING FOR POSTERITY BUT DON'T NEED DUE TO IMPORTS
    # Einc <- cu[, 2] - cuprev[, 2]
    # print(which(Einc > 0))
    # stopifnot(all(rowSums(uprev[Einc > 0, c(3, 5, 6, 8)]) > 0))
}

## source Rcpp PF code
sourceCpp("../APF1.cpp")

## source function to run PF and return log-likelihood
source("../APF1.R")

## read in simulated data and generate incidence curves
data <- readRDS("../outputs/disSims.rds")

## read in parameters, remove guff and reorder
pars <- readRDS("../wave1/disease.rds") %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta"), repeats)) %>%
    select(nu, nuA, !output) %>%
    as.data.frame()

## read in contact matrix
contact <- read_csv("../inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## read in initial conditions
u1 <- readRDS("../outputs/u1.rds")
u1_moves <- readRDS("../outputs/u1_moves.rds")

## set seed for reproducibility
set.seed(666)

## set up plot data
plot_data <- pivot_longer(filter(data, t <= 100), !t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(obs = grepl("obs", var)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    group_by(t, var, age, obs) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    mutate(age = gsub("obs", "", age))
 
for(k in 6) {
    ## run model with no model discrepancy
    runs_nomd <- APF1(pars[k, ], C = contact, data = data, u1_moves = u1_moves,
        u1 = u1, ndays = 100, npart = 10, MD = FALSE, saveAll = TRUE)
    
    ## plot particle estimates of states (unweighted)
    sims_nomd <- map(runs_nomd$particles[[1]], ~{
        map(., ~{
            x <- apply(., c(1, 2), sum)
            colnames(x) <- paste0("age", 1:ncol(x))
            as_tibble(x) %>%
            mutate(var = c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"))
        }) %>%
        bind_rows(.id = "particle")
    }) %>%
    bind_rows(.id = "t") %>%
    pivot_longer(!c(particle, t, var), names_to = "age", values_to = "n") %>%
    mutate(age = as.numeric(gsub("age", "", age))) %>%
    mutate(t = as.numeric(t) - 1) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var))
    
#    ## checks
#    plan(multisession, workers = 8)
#    map(runs_nomd$particles[[1]], ~{
#       out <- future_map(., ~{
#           out <- list()
#           for(i in 1:dim(.)[3]) {
#               temp <- t(.[, , i]) %>%
#                   as_tibble()
#               colnames(temp) <- c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
#               temp <- mutate(temp, age = 1:nrow(temp))
#               out[[i]] <- temp
#           }
#           out
#       })
#       out <- transpose(out) %>%
#           future_map(~bind_rows(., .id = "t"))
#       ## now generate incidence counts
#       outcum <- future_map(out, ~{
#           group_by(., age) %>%
#               mutate(across(c(DI, RI, DH, RH, RA), ~. - lag(., default = 0))) %>%
#               mutate(H = H - lag(H, default = 0) + DH + RH) %>%
#               mutate(I2 = I2 - lag(I2, default = 0) + RI) %>%
#               mutate(I1 = I1 - lag(I1, default = 0) + I2 + DI + H) %>%
#               mutate(P = P - lag(P, default = 0) + I1) %>%
#               mutate(A = A - lag(A, default = 0) + RA) %>%
#               mutate(E = E - lag(E, default = 0) + A + P) %>%
#               mutate(across(E:DH, cumsum)) %>%
#               ungroup()
#       })
#       out <- future_map(out, ~{
#           group_by(., age) %>%
#               nest() %>%
#               pluck("data")
#       })
#       outcum <- future_map(outcum, ~{
#           group_by(., age) %>%
#               nest() %>%
#               pluck("data")
#       })
#       ## check counts
#       future_map2(out, outcum, function(u, cu) {
#           map2(u, cu, function(u, cu) {
#               checkCounts(u, cu, sum(u[1, -1]))
#           })
#       })
#    })
#    plan(multisession, workers = 1)

    gc()
    
    ## repeat but adding some model discrepancy
    runs_md <- APF1(pars[k, ], C = contact, data = data, u1_moves = u1_moves,
        u1 = u1, ndays = 100, npart = 10, MD = TRUE, a_dis = 0.05, b_dis = 0.05, 
        a1 = 0.01, a2 = 0.2, b = 0.001, saveAll = TRUE)
    
    ## plot particle estimates of states (unweighted)
    sims_md <- map(runs_md$particles[[1]], ~{
        map(., ~{
            x <- apply(., c(1, 2), sum)
            colnames(x) <- paste0("age", 1:ncol(x))
            as_tibble(x) %>%
            mutate(var = c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"))
        }) %>%
        bind_rows(.id = "particle")
    }) %>%
    bind_rows(.id = "t") %>%
    pivot_longer(!c(particle, t, var), names_to = "age", values_to = "n") %>%
    mutate(age = as.numeric(gsub("age", "", age))) %>%
    mutate(t = as.numeric(t) - 1) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var))
    
#    ## checks
#    plan(multisession, workers = 8)
#    map(runs_md$particles[[1]], ~{
#       out <- future_map(., ~{
#           out <- list()
#           for(i in 1:dim(.)[3]) {
#               temp <- t(.[, , i]) %>%
#                   as_tibble()
#               colnames(temp) <- c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
#               temp <- mutate(temp, age = 1:nrow(temp))
#               out[[i]] <- temp
#           }
#           out
#       })
#       out <- transpose(out) %>%
#           future_map(~bind_rows(., .id = "t"))
#       ## now generate incidence counts
#       outcum <- future_map(out, ~{
#           group_by(., age) %>%
#               mutate(across(c(DI, RI, DH, RH, RA), ~. - lag(., default = 0))) %>%
#               mutate(H = H - lag(H, default = 0) + DH + RH) %>%
#               mutate(I2 = I2 - lag(I2, default = 0) + RI) %>%
#               mutate(I1 = I1 - lag(I1, default = 0) + I2 + DI + H) %>%
#               mutate(P = P - lag(P, default = 0) + I1) %>%
#               mutate(A = A - lag(A, default = 0) + RA) %>%
#               mutate(E = E - lag(E, default = 0) + A + P) %>%
#               mutate(across(E:DH, cumsum)) %>%
#               ungroup()
#       })
#       out <- future_map(out, ~{
#           group_by(., age) %>%
#               nest() %>%
#               pluck("data")
#       })
#       outcum <- future_map(outcum, ~{
#           group_by(., age) %>%
#               nest() %>%
#               pluck("data")
#       })
#       ## check counts
#       future_map2(out, outcum, function(u, cu) {
#           map2(u, cu, function(u, cu) {
#               checkCounts(u, cu, sum(u[1, -1]))
#           })
#       })
#    })
#    plan(multisession, workers = 1)
    
    ## plot both together nationally
    p <- mutate(sims_nomd, type = "No MD") %>%
        rbind(mutate(sims_md, type = "MD")) %>%
        group_by(t, var, type) %>%
        summarise(
            LCI = quantile(n, probs = 0.025),
            LQ = quantile(n, probs = 0.25),
            median = median(n),
            UQ = quantile(n, probs = 0.75),
            UCI = quantile(n, probs = 0.975),
            .groups = "drop"
        ) %>%
        mutate(var = gsub("one", "1", var)) %>%
        mutate(var = gsub("two", "2", var)) %>%
        ggplot(aes(x = t)) +
            geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = type), alpha = 0.5) +
            geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = type), alpha = 0.5) +
            geom_line(aes(y = median, colour = type)) +
            geom_line(aes(y = n), data = filter(plot_data, !obs), col = "red", linetype = "dashed") +
            geom_line(aes(y = n), data = filter(plot_data, obs), col = "blue", linetype = "dashed") +
            facet_grid(var ~ age, scales = "free") +
            xlab("Days") +
            ylab("Counts")
    ggsave(paste0("sims_combined_", k, "_APF1.pdf"), p, width = 10, height = 10)
}
