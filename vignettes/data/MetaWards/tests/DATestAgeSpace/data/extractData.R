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
    mutate(across(starts_with("deaths_"), ~ifelse(is.na(.), 0, .))) %>%
    mutate(date = as.numeric(date - tstart)) %>%
    rename(t = date)
    
death_region <- readRDS("death_region.rds") %>%
    filter(date >= tstart & date <= tstop) %>%
    complete(date = seq(tstart, tstop, by = 1)) %>%
    mutate(across(starts_with("deaths_"), ~ifelse(is.na(.), 0, .))) %>%
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
    rename(t = date)## read in commuter data
    
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
saveRDS(lookup, "outputs/lookup.rds")
system("cp age_lookup.rds outputs")

