## load libraries
library(tidyverse)
library(lubridate)
library(abind)

## set dates
tstart <- dmy("15/02/2020")
tstop <- dmy("23/03/2020")

## set output name
outnm <- "Feb"

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
system(paste0("rm -r outputs", outnm))
dir.create(paste0("outputs", outnm))
saveRDS(death_lad, paste0("outputs", outnm, "/cumDeath_lad.rds"))
saveRDS(death_region, paste0("outputs", outnm, "/cumDeath_age_region.rds"))
saveRDS(nhsregion_cumadage, paste0("outputs", outnm, "/cumHospAd_age_nhsregion.rds"))
saveRDS(nhsregion_hosp, paste0("outputs", outnm, "/hosp_nhsregion.rds"))
saveRDS(lookup, paste0("outputs", outnm, "/lookup.rds"))
system(paste0("cp age_lookup.rds outputs", outnm))

## write out top LADs
pivot_longer(death_lad, !t, names_to = "lad") %>%
    mutate(lad = as.numeric(gsub("deaths_", "", lad))) %>%
    group_by(lad) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    arrange(desc(value)) %>%
    slice(1:10) %>%
    select(lad) %>%
    write_delim(file = paste0("outputs", outnm, "/lads_outputs", outnm, ".txt"), col_names = FALSE)
    

