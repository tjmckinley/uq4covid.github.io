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
disSims <- select(disSims, t, starts_with("D") | starts_with("H")) %>%
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
    group_by(age, lad) %>%
    mutate(DH = DH - lag(DH, default = 0)) %>%
    mutate(DI = DI - lag(DI, default = 0)) %>%
    ungroup() %>%
    mutate(across(c(DI, DH), ~rtruncnorm(
        n(),
        a = 0,
        b = Inf,
        mean = (a1 - a2) + (b1 - b2 + 1) * .,
        sd = sqrt((a1 + a2) + (b1 + b2) * .)
    ))) %>%
    mutate(H = rtruncnorm(
        n(),
        a = 0,
        b = Inf,
        mean = -(a1 - a2) + (b1 - b2 + 1) * H,
        sd = sqrt(3 * (a1 + a2) + (b1 + b2) * H)
    )) %>%
    mutate(across(c(DI, DH, H), ~round(.))) %>%
    arrange(t, age, lad)

## extract deaths data in the correct format
DI <- mutate(disSims, name = paste0("DI_", age, "_", lad)) %>%
    select(t, DI, name) %>%
    pivot_wider(names_from = name, values_from = DI) %>%
    mutate(across(starts_with("DI"), ~cumsum(.)))
saveRDS(DI, paste0(newoutputdir, "/cumDI_age_lad.rds"))

## extract deaths data in the correct format
DH <- mutate(disSims, name = paste0("DH_", age, "_", lad)) %>%
    select(t, DH, name) %>%
    pivot_wider(names_from = name, values_from = DH) %>%
    mutate(across(starts_with("DH"), ~cumsum(.)))
saveRDS(DH, paste0(newoutputdir, "/cumDH_age_lad.rds"))

## extract hospitalisation data in the correct format
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

