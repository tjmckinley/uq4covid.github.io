## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(sitmo)
library(patchwork)

## source Rcpp PF code
sourceCpp("../BPF.cpp")
sourceCpp("../APF1.cpp")
sourceCpp("../TPF.cpp")

## source function to run PF and return log-likelihood
source("../BPF.R")
source("../APF1.R")
source("../IAPF.R")

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

## set number of replicates
nreps <- 10

## extract parameters
pars1 <- slice(pars, 1:9) %>%
    mutate(id = 1:n())

## expand to replicates
pars <- pars1
if(nreps > 1) {
    for(i in 2:nreps) {
        pars <- rbind(pars, pars1)
    }
}
rm(pars1)
pars <- arrange(pars, id)

## set number of particles
npart <- 10

## run BPF
runs_bpf <- BPF(pars, C = contact, data = data, u1_moves = u1_moves,
    u1 = u1, ndays = 50, npart = 10, a_dis = 0.05, b_dis = 0.05,
    a1 = 0.01, a2 = 0.2, b = 0.001, saveAll = NA)
    
save.image("runs.RData")
    
## run APF1
runs_apf1 <- APF1(pars, C = contact, data = data, u1_moves = u1_moves,
    u1 = u1, ndays = 50, npart = 10, a_dis = 0.05, b_dis = 0.05,
    a1 = 0.01, a2 = 0.2, b = 0.001, saveAll = NA)
    
save.image("runs1.RData")

## run TPF
runs_tpf <- IAPF(pars, C = contact, data = data, u1_moves = u1_moves,
    u1 = u1, ndays = 50, npart = 10, a_dis = 0.05, b_dis = 0.05, kmax = 3,
    a1 = 0.01, a2 = 0.2, b = 0.001, saveAll = NA)
    
save.image("runs2.RData")
  
## collapse to data frame and plot
runs <- tibble(BPF = runs_bpf, APF1 = runs_apf1, TPF = runs_tpf) %>%
    cbind(select(pars, id)) %>%
    pivot_longer(!id, names_to = "Type")

## density plot of PF estimates
p <- ggplot(runs) +
    geom_density(aes(x = value, fill = Type), alpha = 0.5) +
    facet_wrap(~id, scales = "free") +
    ggtitle(paste0("Num. of replicates = ", nreps, " Num. particles = ", npart)) +
    ylab("Density")
    
ggsave(paste0("PF_dens_", nreps, "_ll.pdf"), p, width = 10, height = 10)
