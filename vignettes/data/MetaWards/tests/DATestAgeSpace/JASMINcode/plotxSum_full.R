## load libraries
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 3)
        wave <- args[1]
        hash <- as.numeric(args[2])
        outputs <- args[3]
    } else {
        stop("No arguments")
    }
} else {
    ## set wave number and hash
    wave <- "1"
    hash <- 1
    outputs <- "outputs"
}

## read in LADs from file
## (lad = NA gives national plots, else give vector of lads)
if(file.exists(paste0(outputs, "/lads_", outputs, ".txt"))) {
    lads <- as.numeric(readLines(paste0(outputs, "/lads_", outputs, ".txt")))
} else {
    lads <- NA
}

## set writeExt (eventually pass in as argument)
writeExt <- TRUE

## extract file names
folder <- paste0("wave", wave, "/saveOut_wave", wave, "_", hash)

## lookup table for classes
class_lookup <- data.frame(var = c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH", "DIobs", "DHobs")) %>%
    mutate(class = 0:(n() - 1))

if(!writeExt) {
    stop("Not yet implemented for new model and writeExt = FALSE")
} else {
    ## extract file names
    files <- list.files(folder)
    files <- files[grep("p_[0-9]*.csv.bz2", files)]
    files <- paste0(folder, "/", files)
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_[0-9]*.csv.bz2", files1)]
    if(length(files1) > 0) {
        files1 <- paste0(folder, "_cont/", files1)
        files <- c(files, files1)
    }

    ## load in runs and group at the national level
    sims_md <- map(files, function(y) {
            read.csv(y, header = TRUE) %>%
            group_by(time, class) %>%
            summarise(across(!lad, sum), .groups = "drop")
        })
    names(sims_md) <- gsub("^.*_([0-9]*).csv.bz2", "\\1", files)    
    sims_md <- bind_rows(sims_md, .id = "particle") %>%    
        inner_join(class_lookup, by = "class") %>%
        select(!class) %>%
        rename(t = time) %>%
        pivot_longer(!c(particle, t, var), names_to = "age", values_to = "n") %>%
        mutate(age = as.numeric(gsub("age", "", age)))
        
    ## save output
    saveRDS(sims_md, paste0("wave", wave, "/plotSum_", hash, "_natFull.rds"))
            
    ###################################################
    #######     LAD-age-level observations      #######
    ###################################################

    ## extract file names
    files <- list.files(folder)
    files <- files[grep("p_age_lads_[0-9]*.csv.bz2", files)]
    files <- paste0(folder, "/", files)
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_age_lads_[0-9]*.csv.bz2", files1)]
    if(length(files1) > 0) {
        files1 <- paste0(folder, "_cont/", files1)
        files <- c(files, files1)
    }

    ## load in runs and group at the national level
    sims_md <- map(files, read.csv, header = TRUE)
    names(sims_md) <- gsub("^.*_([0-9]*).csv.bz2", "\\1", files)
    sims_md <- bind_rows(sims_md, .id = "particle") %>%
        pivot_longer(!c(time, lad, particle)) %>%
        mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', name)) %>%
        mutate(name = gsub('^([^_]*)_.*', '\\1', name)) %>%
        mutate(age = as.numeric(gsub("age", "", age))) %>%
        pivot_wider(names_from = name, values_from = value) %>%
        arrange(particle, time, age, lad) %>%
        group_by(particle, age, lad) %>%
        mutate(deaths = cumsum(deaths), hosp = cumsum(hosp)) %>%
        ungroup() %>%
        rename(t = time)
        
    ## produce lad-level plots if required
    if(!is.na(lads[1])) {
        ## save output
        saveRDS(sims_md, paste0("wave", wave, "/plotSum_", hash, "_age_lads.rds"))
    }
    
    ## now aggregate to national level 
    sims_md <- group_by(sims_md, particle, t, age) %>%
        summarise(deaths = sum(deaths), hosp = sum(hosp), .groups = "drop")
        
    ## save output
    saveRDS(sims_md, paste0("wave", wave, "/plotSum_", hash, "_natAgeDeathsHosp.rds"))
}

print("Finished")

