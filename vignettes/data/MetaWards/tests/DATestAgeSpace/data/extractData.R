## load libraries
library(tidyverse)
library(lubridate)
library(abind)

## set dates
tstart <- dmy("15/01/2020")
tstop <- dmy("23/03/2020")

## load data, extract relevant dates, fill missing
## values - assume deaths are zero when not observed
## treat others as missing
death_lad <- readRDS("death_lad.rds") %>%
    filter(date >= tstart & date <= tstop) %>%
    complete(date = seq(tstart, tstop, by = 1)) %>%
    mutate(across(starts_with("D_"), ~ifelse(is.na(.), 0, .))) %>%
    mutate(date = as.numeric(date - tstart)) %>%
    rename(t = date)
    
death_region <- readRDS("death_region.rds") %>%
    filter(date >= tstart & date <= tstop) %>%
    complete(date = seq(tstart, tstop, by = 1)) %>%
    mutate(across(starts_with("D_"), ~ifelse(is.na(.), 0, .))) %>%
    mutate(date = as.numeric(date - tstart)) %>%
    rename(t = date)
    
nhsregion_cumadage <- readRDS("nhsregion_cumadage.rds") %>%
    filter(date >= tstart & date <= tstop) %>%
    complete(date = seq(tstart, tstop, by = 1)) %>%
    mutate(date = as.numeric(date - tstart)) %>%
    rename(t = date)
    
nhsregion_hosp <- readRDS("nhsregion_hosp.rds") %>%
    filter(date >= tstart & date <= tstop) %>%
    complete(date = seq(tstart, tstop, by = 1)) %>%
    mutate(date = as.numeric(date - tstart)) %>%
    rename(t = date)
    
## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## add age probabilities
ageProbs <- read_csv("../inputs/age_seeds.csv", col_names = FALSE)$X2

## read in commuter data
EW19 <- read_delim("../inputs/EW19.dat", delim = " ", col_names = FALSE)

## add seeding information
EW19 <- mutate(EW19, X4 = 0)
for(i in 1:10) {
    val <- 0
    while(val == 0) {
        seed <- sample(which(EW19[, 1] == 317), 1)
        if((EW19$X4[seed] + 1) <= EW19$X3[seed]) {
            EW19$X4[seed] <- EW19$X4[seed] + 1
            val <- 1
        }
    }
}
## expand to deal with age-classes
u1 <- apply(EW19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[2, ] <- smart_round(ageProbs * x[4])
        u[1, ] <- smart_round(ageProbs * x[3]) - u[2, ]
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)
u1_moves <- as.matrix(EW19[, 1:2])

## read in player data
PM19 <- read_delim("../inputs/PlayMatrix19.dat", delim = " ", col_names = FALSE)
PlaySize19 <- read_delim("../inputs/PlaySize19.dat", delim = " ", col_names = FALSE)
  
## expand to deal with age-classes
u2 <- apply(PlaySize19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[1, ] <- smart_round(ageProbs * x[2])
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)
    
## load lookup
lookup <- readRDS("lookup.rds")
lookup <- select(lookup, starts_with("FID"))
    
## save outputs
system("rm -r outputs")
dir.create("outputs")
saveRDS(death_lad, "outputs/cumDeath_lad.rds")
saveRDS(death_region, "outputs/cumDeath_age_region.rds")
saveRDS(nhsregion_cumadage, "outputs/cumHospAd_age_nhsregion.rds")
saveRDS(nhsregion_hosp, "outputs/hosp_nhsregion.rds")
saveRDS(u1, "outputs/u1.rds")
saveRDS(u1_moves, "outputs/u1_moves.rds")
saveRDS(u2, "outputs/u2.rds")
saveRDS(as.matrix(PM19), "outputs/u2_moves.rds")
saveRDS(lookup, "outputs/lookup.rds")
system("cp age_lookup.rds outputs")

