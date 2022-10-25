## load libraries
library(tidyverse)

######################################
#######       Load data       ########
######################################

## load LAD data
lad <- read_csv("ltla_2022-07-06.csv")

## load regional data
region <- read_csv("region_2022-07-06.csv")

## load movement data
lookup <- read_csv("../inputs/LAD19_Lookup.csv")

## load lad to region lookup
region_lookup <- read_csv("Local_Authority_District_to_Region_(April_2019)_Lookup_in_England.csv")

## load NHS region data
nhsregion_hosp <- read_csv("nhsRegion_hospCases_2022-07-13.csv")
nhsregion_cumadage <- read_csv("nhsRegion_cumAdByAge_2022-07-13.csv")

######################################
#######        LAD data       ########
######################################

## check for LADs in LTLA table that are not in lookup
nomaplookup <- anti_join(lad, lookup, by = c("areaCode" = "LAD19CD"))
nomaplookup <- unique(nomaplookup$areaCode)

## check if these can be removed (i.e. if they are not in England or
## Wales then the model doesn't account for them)
table(substr(nomaplookup, 1, 1))

## remove lads in Scotland and NI from lad data
lad <- semi_join(lad, lookup, by = c("areaCode" = "LAD19CD"))

## check for LADs in lookup that are not in LTLA table 
nomaplookup <- anti_join(lookup, lad, by = c("LAD19CD" = "areaCode"))
nomaplookup <- unique(nomaplookup$LAD19CD)

## check English LADs
filter(lookup, LAD19CD %in% nomaplookup[grep("^E", nomaplookup)])

## Isles of Scilly and City of London not present in the lad data, but it 
## looks like these can be mapped to: Cornwall and Isles of Scilly
## and Hackney and City of London
lad[grep("Scilly", lad$areaName), 2][1, ]
lad[grep("City of London", lad$areaName), 2][1, ]

## amend code in lookup table to match ltla code for the
## two missing regions
lookup$LAD19CD[lookup$LAD19NM == "Isles of Scilly"] <- lad$areaCode[grep("Scilly", lad$areaName)][1]
lookup$LAD19CD[lookup$LAD19NM == "City of London"] <- lad$areaCode[grep("City of London", lad$areaName)][1]

## check for LADs in lookup that are not in LTLA table 
nomaplookup <- anti_join(lookup, lad, by = c("LAD19CD" = "areaCode"))
nomaplookup <- unique(nomaplookup$LAD19CD)

## there are some Welsh LADs that are not present in death data
## this could be because they have no recorded deaths

## see if there are any non-missing Welsh LADs in the death data
temp <- unique(lookup$LAD19CD)
temp <- temp[grep("^W", temp)]
identical(sort(unique(nomaplookup)), sort(temp))

## hence there are no Welsh LADs in the death data,
## but we know there were deaths in Wales, so most likely
## explanation is that they aren't recorded in the data
## and thus we can only calibrate to English LADs.
## Welsh LADs can remain in the model but must be removed
## when calibrating

## tidy up data set
lad <- select(lad, !areaType) %>%
    rename(cumDeaths = cumDeaths28DaysByDeathDate) %>%
    arrange(areaCode, date)
    
## fill in gaps to get complete counts
lad <- complete(lad, date, areaCode) %>%
    mutate(cumDeaths = ifelse(date == min(date) & is.na(cumDeaths), 0, cumDeaths)) %>%
    arrange(areaCode, date) %>%
    group_by(areaCode) %>%
    fill(cumDeaths) %>% 
    ungroup()
    
## set up unique identifiers
death_lookup <- distinct(lad, areaCode, areaName) %>%
    filter(!is.na(areaName)) %>%
    mutate(FID = as.numeric(factor(areaCode)))
    
## save death data lookup
saveRDS(death_lookup, "death_lookup.rds")
    
## amend lad data to use new ID
lad <- select(lad, !areaName) %>%
    inner_join(death_lookup, by = "areaCode") %>%
    select(FID, date, cumDeaths) %>%
    arrange(FID, date)

## add line to lookup file denoting which LADs we can calibrate to
lookup <- left_join(lookup, select(death_lookup, !areaName), by = c("LAD19CD" = "areaCode")) %>%
    rename(FID_death = FID.y, FID = FID.x)

## expand to correct format for model
lad <- mutate(lad, FID = paste0("deaths_", FID)) %>%
    pivot_wider(names_from = FID, values_from = cumDeaths)

## save lad-level data
saveRDS(lad, "death_lad.rds")

######################################
#######     regional data     ########
######################################

## check region lookup against regional data
anti_join(region_lookup, region, by = c("RGN19CD" = "areaCode"))
anti_join(region, region_lookup, by = c("areaCode" = "RGN19CD"))

## check region_lookup against lookup
anti_join(region_lookup, lookup, by = "LAD19CD")

## amend region lookup to match lookup for
## Isles of Scilly and City of London
region_lookup$LAD19CD[region_lookup$LAD19NM == "Isles of Scilly"] <- lookup$LAD19CD[lookup$LAD19NM == "Isles of Scilly"]
region_lookup$LAD19CD[region_lookup$LAD19NM == "City of London"] <- lookup$LAD19CD[lookup$LAD19NM == "City of London"]

## check region_lookup against lookup
anti_join(region_lookup, lookup, by = "LAD19CD")

## so all LADs in region lookup exist in lookup

## create region lookup for later
temp <- distinct(region_lookup, RGN19CD, RGN19NM) %>%
    arrange(RGN19CD) %>%
    mutate(FID = 1:n()) %>%
    select(FID, RGN19CD, RGN19NM)
    
## save region lookup
saveRDS(temp, "region_lookup.rds")

## add region to lookup
lookup <- distinct(region_lookup, LAD19CD, RGN19CD) %>%
    left_join(temp, by = "RGN19CD") %>%
    select(LAD19CD, FID_region = FID) %>%
    right_join(lookup, by = "LAD19CD") %>%
    select(FID, LAD19CD, LAD19NM, FID_death, FID_region) %>%
    arrange(FID)
region_lookup <- temp

## now tidy up regional death data
region <- select(region, !c(areaName, areaType, starts_with("rolling"))) %>%
    inner_join(region_lookup, by = c("areaCode" = "RGN19CD")) %>%
    select(FID, date, age, deaths) %>%
    complete(date, FID, age, fill = list(deaths = 0)) %>%
    arrange(date) %>%
    group_by(FID, age) %>%
    mutate(cumDeaths = cumsum(deaths)) %>%
    ungroup() %>%
    arrange(age, FID, date) %>%
    select(!deaths)
    
## convert to similar age-classes as model
region <- filter(region, age != "00_59") %>%
    filter(age != "60+") %>%
    mutate(age = ifelse(age == "00_04", "<5", age)) %>%
    mutate(age = ifelse(age == "05_09", "05-19", age)) %>%
    mutate(age = ifelse(age == "10_14", "05-19", age)) %>%
    mutate(age = ifelse(age == "15_19", "05-19", age)) %>%
    mutate(age = ifelse(age == "20_24", "20-29", age)) %>%
    mutate(age = ifelse(age == "25_29", "20-29", age)) %>%
    mutate(age = ifelse(age == "30_34", "30-39", age)) %>%
    mutate(age = ifelse(age == "35_39", "30-39", age)) %>%
    mutate(age = ifelse(age == "40_44", "40-49", age)) %>%
    mutate(age = ifelse(age == "45_49", "40-49", age)) %>%
    mutate(age = ifelse(age == "50_54", "50-59", age)) %>%
    mutate(age = ifelse(age == "55_59", "50-59", age)) %>%
    mutate(age = ifelse(age == "60_64", "60-69", age)) %>%
    mutate(age = ifelse(age == "65_69", "60-69", age)) %>%
    mutate(age = ifelse(age == "70_74", "70+", age)) %>%
    mutate(age = ifelse(age == "75_79", "70+", age)) %>%
    mutate(age = ifelse(age == "80_84", "70+", age)) %>%
    mutate(age = ifelse(age == "85_89", "70+", age)) %>%
    mutate(age = ifelse(age == "90+", "70+", age)) %>%
    group_by(FID, age, date) %>%
    summarise(cumDeaths = sum(cumDeaths), .groups = "drop") %>%
    arrange(age, FID, date)
    
## expand to correct format for model
region <- mutate(region, age = as.numeric(factor(age))) %>%
    mutate(age = paste0("deaths_", age)) %>%
    unite(FID, age, FID, sep = "_") %>%
    pivot_wider(names_from = FID, values_from = cumDeaths)

## save lad-level data
saveRDS(region, "death_region.rds")

######################################
#######   check death data    ########
######################################

tempDeaths <- pivot_longer(lad, !date, names_to = "lad", values_to = "deaths") %>%
    mutate(lad = as.numeric(gsub("deaths_", "", lad))) %>%
    inner_join(lookup, by = c("lad" = "FID_death")) %>%
    inner_join(region_lookup, by = c("FID_region" = "FID")) %>%
    group_by(RGN19NM, date) %>%
    summarise(deaths = sum(deaths), .groups = "drop")

tempDeathsRegion <- pivot_longer(region, !date, names_to = "region", values_to = "deaths") %>%
    mutate(region = gsub("deaths_", "", region)) %>%
    separate(region, c("age", "region"), sep = "_") %>%
    group_by(region, date) %>%
    summarise(deaths = sum(deaths), .groups = "drop") %>%
    mutate(region = as.numeric(region)) %>%
    inner_join(region_lookup, by = c("region" = "FID"))

p <- ggplot(tempDeaths, aes(x = date, y = deaths, colour = RGN19NM)) +
    geom_line() +
    geom_line(data = tempDeathsRegion, linetype = "dashed")
ggsave("aggregatedDeathsPlot.pdf", p)

######################################
#######    NHS region data    ########
######################################

## check admissions against age against total admissions
p <- group_by(nhsregion_cumadage, areaCode, areaName, date) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    inner_join(nhsregion_hosp, by = c("areaCode", "areaName", "date")) %>%
    select(areaName, date, value, cumAdmissions) %>%
    pivot_longer(!c(areaName, date), names_to = "type") %>%
    ggplot(aes(x = date, y = value, colour = areaName, linetype = type)) +
        geom_line()
ggsave("aggregatedHospAdmissionsPlot.pdf", p)

## currently treat age-aggregated counts as if they matched, and just
## add some more observation error (can do something better later on)

## generate NHS region lookup
nhsregion_lookup <- distinct(nhsregion_hosp, areaCode, areaName) %>%
    mutate(FID = 1:n())
saveRDS(nhsregion_lookup, "nhsregion_lookup.rds")

## update lookup
lookup <- left_join(lookup, 
    mutate(region_lookup, RGN19NM = ifelse(RGN19NM == "Yorkshire and The Humber", "North East and Yorkshire", RGN19NM)) %>%
        mutate(RGN19NM = ifelse(RGN19NM == "North East", "North East and Yorkshire", RGN19NM)) %>%
        mutate(RGN19NM = ifelse(RGN19NM == "East Midlands", "Midlands", RGN19NM)) %>%
        mutate(RGN19NM = ifelse(RGN19NM == "West Midlands", "Midlands", RGN19NM)) %>%
        inner_join(nhsregion_lookup, by = c("RGN19NM" = "areaName")) %>%
        select(FID.x, FID.y),
    by = c("FID_region" = "FID.x")) %>%
    rename(FID_nhsregion = FID.y)

## wrangle nhs regional data into correct format
nhsregion_hosp <- inner_join(nhsregion_hosp, nhsregion_lookup, by = "areaName") %>%
    select(FID, date, hospitalCases) %>%
    mutate(FID = paste0("hosp_", FID)) %>%
    pivot_wider(names_from = FID, values_from = hospitalCases) %>%
    arrange(date)
saveRDS(nhsregion_hosp, "nhsregion_hosp.rds")

## generate NHS age lookup
age_lookup <- data.frame(age = sort(unique(nhsregion_cumadage$age))) %>%
    mutate(age = ifelse(age == "6_to_17", "06_to_17", age)) %>%
    mutate(age = ifelse(age == "65_to_84", "65+", age)) %>%
    mutate(age = ifelse(age == "85+", "65+", age)) %>%
    distinct(age) %>%
    arrange(age) %>%
    mutate(FID = 1:n())
saveRDS(age_lookup, "nhsage_lookup.rds")

## wrangle nhs regional data into correct format
nhsregion_cumadage <- mutate(nhsregion_cumadage, age = ifelse(age == "6_to_17", "06_to_17", age)) %>%
    mutate(age = ifelse(age == "65_to_84", "65+", age)) %>%
    mutate(age = ifelse(age == "85+", "65+", age)) %>%
    group_by(areaName, date, age) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    inner_join(nhsregion_lookup, by = "areaName") %>%
    inner_join(age_lookup, by = "age") %>%
    select(FID.y, FID.x, date, value) %>%
    arrange(FID.y, FID.x, date) %>%
    mutate(FID.y = paste0("hospInc_", FID.y)) %>%
    unite(FID, FID.y, FID.x, sep = "_") %>%
    pivot_wider(names_from = FID, values_from = value) %>%
    arrange(date)
saveRDS(nhsregion_cumadage, "nhsregion_cumadage.rds")

######################################
#######      save lookup      ########
######################################

saveRDS(lookup, "lookup.rds")

## set up approximate age mapping
age_lookup <- data.frame(FID_death = 1:8, FID_nhsregion = c(1, 2, 3, 3, 3, 3, 4, 4))
saveRDS(age_lookup, "age_lookup.rds")

