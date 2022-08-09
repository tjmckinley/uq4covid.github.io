## load libraries
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 2)
        wave <- as.numeric(args[1])
        hash <- as.numeric(args[2])
    } else {
        stop("No arguments")
    }
} else {
    ## set wave number and hash
    wave <- 1
    hash <- 1
}

## read in LADs from file
## (lad = NA gives national plots, else give vector of lads)
if(file.exists("JASMINcode/lads.txt")) {
    lads <- as.numeric(readLines("JASMINcode/lads.txt"))
} else {
    lads <- NA
}

## set writeExt (eventually pass in as argument)
writeExt <- TRUE

## extract file names
folder <- paste0("wave", wave, "/saveOut_", hash)

## lookup table for classes
class_lookup <- data.frame(var = c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH", "DIobs", "DHobs")) %>%
    mutate(class = 0:(n() - 1))

if(!writeExt) {
    stop("Not yet implemented for new model and writeExt = FALSE")
} else {
    ## extract file names
    files <- list.files(folder)
    files <- files[grep("p_[0-9]*.csv", files)]
    files <- paste0(folder, "/", files)
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_[0-9]*.csv", files1)]
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
    names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)    
    sims_md <- bind_rows(sims_md, .id = "particle") %>%    
        inner_join(class_lookup, by = "class") %>%
        select(!class) %>%
        rename(t = time) %>%
        pivot_longer(!c(particle, t, var), names_to = "age", values_to = "n") %>%
        mutate(age = as.numeric(gsub("age", "", age)))
        
    ## save output
    saveRDS(sims_md, paste0("wave", wave, "/plotSum_", hash, "_natFull.rds"))
            
    ###############################################
    #######     LAD-level observations      #######
    ###############################################

    ## extract file names
    files <- list.files(folder)
    files <- files[grep("p_lads_[0-9]*.csv", files)]
    files <- paste0(folder, "/", files)
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_lads_[0-9]*.csv", files1)]
    if(length(files1) > 0) {
        files1 <- paste0(folder, "_cont/", files1)
        files <- c(files, files1)
    }

    ## load in runs and group at the national level
    sims_md <- map(files, read.csv, header = TRUE)
    names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)
    sims_md <- bind_rows(sims_md, .id = "particle") %>%
        arrange(particle, time, lad) %>%
        group_by(particle, lad) %>%
        mutate(deaths = cumsum(deaths)) %>%
        group_by(particle, time) %>%
        summarise(deaths = sum(deaths), .groups = "drop") %>%
        rename(t = time, n = deaths)
        
    ## save output
    saveRDS(sims_md, paste0("wave", wave, "/plotSum_", hash, "_natDeaths.rds"))
              
    ###############################################
    #######  age/region-level observations  #######
    ###############################################

    ## load lookup for labels
    region_lookup <- readRDS("data/region_lookup.rds")
        
    ## extract file names
    files <- list.files(folder)
    files <- files[grep("p_age_region_[0-9]*.csv", files)]
    files <- paste0(folder, "/", files)
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_age_region_[0-9]*.csv", files1)]
    if(length(files1) > 0) {
        files1 <- paste0(folder, "_cont/", files1)
        files <- c(files, files1)
    }

    ## load in runs and group at the national level
    sims_md <- map(files, read.csv, header = TRUE)
    names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)
    sims_md <- bind_rows(sims_md, .id = "particle") %>%
        pivot_longer(!c(particle, time, region), names_to = "age", values_to = "n") %>%
        mutate(age = as.numeric(gsub("age", "", age))) %>%
        rename(t = time) %>%
        arrange(particle, region, age, t) %>%
        group_by(particle, region, age) %>%
        mutate(n = cumsum(n)) %>%
        ungroup() %>%
        inner_join(region_lookup, by = c("region" = "FID")) %>%
        select(particle, t, n, age, RGN19NM)
        
    ## save output
    saveRDS(sims_md, paste0("wave", wave, "/plotSum_", hash, "_ageRegionDeaths.rds"))
        
    ###############################################
    #######     NHS region observations     #######
    ###############################################

    ## load lookup for labels
    nhsregion_lookup <- readRDS("data/nhsregion_lookup.rds")

    ## extract file names
    files <- list.files(folder)
    files <- files[grep("p_nhsregion_[0-9]*.csv", files)]
    files <- paste0(folder, "/", files)
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_nhsregion_[0-9]*.csv", files1)]
    if(length(files1) > 0) {
        files1 <- paste0(folder, "_cont/", files1)
        files <- c(files, files1)
    }

    ## load in runs and group at the national level
    sims_md <- map(files, read.csv, header = TRUE)
    names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)

    ## load in runs and group at the national level
    sims_md <- bind_rows(sims_md, .id = "particle") %>%
        rename(t = time, n = hosp) %>%
        inner_join(nhsregion_lookup, by = c("nhsregion" = "FID")) %>%
        select(particle, t, n, areaName) %>%
        arrange(particle, areaName, t)
        
    ## save output
    saveRDS(sims_md, paste0("wave", wave, "/plotSum_", hash, "_nhsregionHosp.rds"))
        
    ###############################################
    #####  NHS age/region-level observations  #####
    ###############################################

    ## extract file names
    files <- list.files(folder)
    files <- files[grep("p_age_nhsregion_[0-9]*.csv", files)]
    files <- paste0(folder, "/", files)
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_age_nhsregion_[0-9]*.csv", files1)]
    if(length(files1) > 0) {
        files1 <- paste0(folder, "_cont/", files1)
        files <- c(files, files1)
    }

    ## load in runs and group at the national level
    sims_md <- map(files, read.csv, header = TRUE)
    names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)

    ## load in runs and group at the national level
    sims_md <- bind_rows(sims_md, .id = "particle") %>%
        pivot_longer(!c(particle, time, nhsregion), names_to = "age", values_to = "n") %>%
        mutate(age = as.numeric(gsub("age", "", age))) %>%
        rename(t = time) %>%
        arrange(particle, nhsregion, age, t) %>%
        group_by(particle, nhsregion, age) %>%
        mutate(n = cumsum(n)) %>%
        ungroup() %>%
        inner_join(nhsregion_lookup, by = c("nhsregion" = "FID")) %>%
        select(particle, t, n, age, areaName)
        
    ## save output
    saveRDS(sims_md, paste0("wave", wave, "/plotSum_", hash, "_ageNhsregionHosp.rds"))
    
    ## produce lad-level plots if required
    if(!is.na(lads[1])) {
      
        ###############################################
        #####          LAD-level plots            #####
        ###############################################

        ## extract file names
        files <- list.files(folder)
        files <- files[grep("p_[0-9]*.csv", files)]
        files <- paste0(folder, "/", files)
        ## extract file names
        files1 <- list.files(paste0(folder, "_cont"))
        files1 <- files1[grep("p_[0-9]*.csv", files1)]
        if(length(files1) > 0) {
            files1 <- paste0(folder, "_cont/", files1)
            files <- c(files, files1)
        }

        ## load in runs and group at the national level
        sims_md <- map(files, function(y, lads) {
                read.csv(y, header = TRUE) %>%
                    filter(lad %in% lads)
            }, lads = lads)
        names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)
        sims_md <- bind_rows(sims_md, .id = "particle") %>% 
            inner_join(class_lookup, by = "class") %>%
            select(!class) %>%
            rename(t = time) %>%
            pivot_longer(!c(particle, t, var, lad), names_to = "age", values_to = "n") %>%
            mutate(age = as.numeric(gsub("age", "", age))) %>%
            mutate(var = gsub("one", "1", var)) %>%
            mutate(var = gsub("two", "2", var)) %>%
            mutate(lad = as.character(lad))
            
        ## save output
        saveRDS(sims_md, paste0("wave", wave, "/plotSum_", hash, "_lads.rds"))
    }
}

print("Finished")

