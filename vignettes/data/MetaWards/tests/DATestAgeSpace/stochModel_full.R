## load libraries
library(tidyverse)
library(truncnorm)

## set seed
seed <- 456
#set.seed(seed)

## create output directory
outputdir <- paste0("outputs", seed)
newoutputdir <- paste0("outputsFull", seed)
if(!dir.exists(outputdir)) {
    stop("'outputdir' doesn't exist")
}
if(dir.exists(newoutputdir)) {
    stop("Can't overwrite existing directory")
}
dir.create(newoutputdir)

## load in simulated data
disSims <- readRDS(paste0(outputdir, "/disSims.rds"))

## load OE terms
fixedInputs <- readLines("wave1/fixedInputs.txt")
a1 <- as.numeric(fixedInputs[6])
a2 <- as.numeric(fixedInputs[7])
b1 <- as.numeric(fixedInputs[8])
b2 <- as.numeric(fixedInputs[9])

## extract lookup
lookup <- readRDS(paste0(outputdir, "/lookup.rds")) %>%
    filter(!is.na(FID_death))

## extract relevant counts and sample observation error
disSims <- select(disSims, t, starts_with("D") | starts_with("H") | starts_with("RH")) %>%
    pivot_longer(!t) %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', name)) %>%
    mutate(lad = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(name = gsub('_.*', '', name)) %>%
    mutate(across(c(age, lad), ~as.numeric(.))) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    inner_join(select(lookup, FID, FID_death), by = c("lad" = "FID")) %>%
    select(!lad) %>%
    rename(lad = FID_death) %>%
    group_by(t, age, lad) %>%
    summarise(across(everything(), ~sum(.)), .groups = "drop") %>%
    mutate(deaths = DI + DH) %>%
    group_by(age, lad) %>%
    mutate(deaths_inc = deaths - lag(deaths, default = 0)) %>%
    mutate(DH = DH - lag(DH, default = 0)) %>%
    mutate(RH = RH - lag(RH, default = 0)) %>%
    mutate(H_inc = H - lag(H, default = 0) + DH + RH) %>%
    ungroup() %>%
    select(t, age, lad, deaths_inc, H_inc, H) %>%
    rename(hosp = H_inc, deaths = deaths_inc) %>%
    mutate(deaths = rtruncnorm(
        n(),
        a = 0,
        b = Inf,
        mean = (a1 - a2) + (b1 - b2 + 1) * deaths,
        sd = sqrt((a1 + a2) + (b1 + b2) * deaths)
    )) %>%
    mutate(hosp = rtruncnorm(
        n(),
        a = 0,
        b = Inf,
        mean = (a1 - a2) + (b1 - b2 + 1) * hosp,
        sd = sqrt((a1 + a2) + (b1 + b2) * hosp)
    )) %>%
    mutate(H = rtruncnorm(
        n(),
        a = 0,
        b = Inf,
        mean = (a1 - a2) + (b1 - b2 + 1) * H,
        sd = sqrt((a1 + a2) + (b1 + b2) * H)
    )) %>%
    mutate(across(c(deaths, hosp, H), ~round(.))) %>%
    arrange(t, age, lad)

## extract deaths data in the correct format
deaths <- mutate(disSims, name = paste0("deaths_", age, "_", lad)) %>%
    select(t, deaths, name) %>%
    pivot_wider(names_from = name, values_from = deaths) %>%
    mutate(across(starts_with("deaths"), ~cumsum(.)))
saveRDS(deaths, paste0(newoutputdir, "/cumDeath_age_lad.rds"))

## extract hospitalisation data in the correct format
hosp <- mutate(disSims, name = paste0("hosp_", age, "_", lad)) %>%
    select(t, hosp, name) %>%
    pivot_wider(names_from = name, values_from = hosp) %>%
    mutate(across(starts_with("hosp"), ~cumsum(.)))
saveRDS(hosp, paste0(newoutputdir, "/cumHosp_age_lad.rds"))

## extract hospital count data in correct format
H <- mutate(disSims, name = paste0("H_", age, "_", lad)) %>%
    select(t, H, name) %>%
    pivot_wider(names_from = name, values_from = H)
saveRDS(H, paste0(newoutputdir, "/H_age_lad.rds"))

## copy over other necessary files
files <- c(
    "death_lookup.rds",
    "disSims.rds",
    "Local_Authority_Districts_\\(December_2019\\)_Boundaries_UK_BUC.zip",
    "lookup.rds",
    "pars.rds",
    "u1*", "u2*",
    "*.pdf"
)
map(files, ~system(paste0("cp ", outputdir, "/", ., " ", newoutputdir, "/")))
system(paste0("cp ", outputdir, "/lads_", outputdir, ".txt ", newoutputdir, "/lads_", newoutputdir, ".txt"))

