## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(sitmo)
library(sf)
library(gganimate)
library(viridis)
library(patchwork)

## set seed
set.seed(666)

## create output directory
dir.create("outputs")

## source simulation function
sourceCpp("PF.cpp")
source("PF.R")

## source Skellam function for observation error
source("trSkellam.R")

## read in parameters, remove guff and reorder
pars <- readRDS("wave1/disease.rds") %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta"), repeats)) %>%
    select(nu, nuA, !output, output)

## read in contact matrix
contact <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## extract parameters for simulation   
pars <- select(slice(pars, 6), !output)

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
ageProbs <- read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2

## read in commuter data
EW19 <- read_delim("inputs/EW19.dat", delim = " ", col_names = FALSE)

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

## write inputs out
saveRDS(u1, "outputs/u1.rds")
saveRDS(u1_moves, "outputs/u1_moves.rds")

## set up stage names
stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH", "DIobs", "DHobs"), ~paste0(., "_", 1:8)) %>%
    map(~map(., ~paste0(., "_", 1:max(EW19[, 1])))) %>%
    reduce(c) %>%
    reduce(c)

## simulate discrete-time model
disSims <- PF(pars, C = contact, data = data, u1_moves = u1_moves,
    u1 = u1, ndays = 100, npart = 24, PF = FALSE)
        
## collapse to data frame
disSims <- map(1:length(disSims$particles[[1]]), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            ## collapse to vector in order: stage, age, LAD
            aperm(x[[i]], 3:1) %>%
            as.vector() %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims$particles[[1]]) %>%
    {do.call("rbind", .)}
colnames(disSims) <- c(stageNms, "rep", "t")
disSims <- as_tibble(disSims)

## extract simulation closest to median
medRep <- select(disSims, !contains("obs")) %>%
    pivot_longer(!c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        median = median(n),
        .groups = "drop"
    ) %>%
    inner_join(
        pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n"),
        by = c("t", "var")
    ) %>%
    mutate(diff = (median - n)^2) %>%
    group_by(rep) %>%
    summarise(diff = sum(diff), .groups = "drop") %>%
    arrange(diff) %>%
    slice(1) %>%
    inner_join(disSims, by = "rep") %>%
    select(!c(diff, rep))

## plot replicates at national level
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    mutate(var = gsub("_[0-9]*$", "", var)) %>%
    group_by(rep, t, var) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(t, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(age = gsub("[^0-9]", "", var)) %>%
    mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var))
    
p1 <- list()     
p1[[1]] <- filter(p, !(var %in% c("DHobs", "DIobs"))) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = median)) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(select(medRep, !contains("obs")), !t, names_to = "var", values_to = "n") %>%
                mutate(var = gsub("_[0-9]*$", "", var)) %>%
                group_by(t, var) %>%
                summarise(n = sum(n), .groups = "drop") %>%
                mutate(age = gsub("[^0-9]", "", var)) %>%
                mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
                mutate(var = gsub("one", "1", var)) %>%
                mutate(var = gsub("two", "2", var)),
            col = "red", linetype = "dashed"
        ) +
        facet_grid(var ~ age, scales = "free") +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle("Truth")
p1[[2]] <- filter(p, var %in% c("DHobs", "DIobs")) %>%
    mutate(var = gsub("obs", "", var)) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = median)) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(select(medRep, t, contains("obs")), !t, names_to = "var", values_to = "n") %>%
                mutate(var = gsub("obs", "", var)) %>%
                mutate(var = gsub("_[0-9]*$", "", var)) %>%
                group_by(t, var) %>%
                summarise(n = sum(n), .groups = "drop") %>%
                mutate(age = gsub("[^0-9]", "", var)) %>%
                mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
                mutate(var = gsub("one", "1", var)) %>%
                mutate(var = gsub("two", "2", var)),
            col = "red", linetype = "dashed"
        ) +
        facet_grid(var ~ age, scales = "free") +
        xlab("Days") + 
        ylab("Counts") + 
        ggtitle("Observed")
p1 <- wrap_plots(p1, nrow = 2, heights = c(0.8, 0.2)) +
    plot_layout(guides = "collect")
ggsave("outputs/simsNational.pdf", p1, width = 10, height = 10)

## plot medRep at LAD level for LADs with largest epidemic load
p <- select(medRep, !contains("obs")) %>%
    pivot_longer(!t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    filter(var == "DI") %>%
    filter(t == max(t)) %>%
    group_by(LAD) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    arrange(desc(n)) %>%
    slice(1:5) %>%
    select(-n)
pobs <- inner_join(p,
        select(medRep, t, contains("obs")) %>%
            pivot_longer(!t, names_to = "var", values_to = "n") %>%
            mutate(var = gsub("obs", "", var)) %>%
            mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
            mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)),
        by = "LAD"
    ) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) 
p <- inner_join(p,
        select(medRep, !contains("obs")) %>%
            pivot_longer(!t, names_to = "var", values_to = "n") %>%
            mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
            mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)),
        by = "LAD"
    ) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var))
p1 <- list()
p1[[1]] <- ggplot(p, aes(x = t, y = n, colour = LAD)) +
    geom_line() +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Truth")
p1[[2]] <- ggplot(pobs, aes(x = t, y = n, colour = LAD)) +
    geom_line() +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed")
p1 <- wrap_plots(p1, nrow = 2, heights = c(0.8, 0.2)) +
    plot_layout(guides = "collect")
ggsave("outputs/simsTopLADs.pdf", p1, width = 10, height = 10)

## save outputs
saveRDS(medRep, "outputs/disSims.rds")
saveRDS(pars, "outputs/pars.rds")

## spatial animation of simulation

## read in shapefile
lad19 <- st_read("inputs/LAD19_shapefile/LAD19_shapefile.shp")

## extract cases over time in each age-group
p <- select(medRep, t, starts_with("DH_")) %>%
    select(!contains("obs")) %>%
    filter(t <= 100) %>%
    pivot_longer(!t, values_to = "counts", names_to = "var") %>%
    separate(var, c("var", "age", "lad"), sep = "_") %>%
    select(!var) %>%
    mutate(age = as.numeric(age), lad = as.numeric(lad)) %>%
    group_by(lad) %>%
    mutate(tot = sum(counts)) %>%
    ungroup() %>%
    filter(tot > 0) %>%
    select(!tot) %>%
    group_by(lad, age) %>%
    mutate(tot = cumsum(counts)) %>%
    group_by(lad, t) %>%
    mutate(tot = sum(tot)) %>%
    ungroup() %>%
    mutate(counts = ifelse(tot == 0, NA, counts)) %>%
    select(!tot)
max_p <- max(p$counts, na.rm = TRUE)
p <- inner_join(lad19, p, by = c("objectid" = "lad")) %>%
    ggplot() +
    geom_sf(aes(fill = counts), colour = NA) +
    facet_wrap(~age) +
    scale_fill_viridis_c(limits = c(0, max_p)) +
    theme_bw()

## add transitions
p <- p + transition_time(t) + ggtitle("Day = {frame_time}")
spatial_gif <- animate(p, nframes = 50, fps = 1, renderer = gifski_renderer())
anim_save("outputs/simsSpatialDH.gif", spatial_gif)

## extract cases over time in each age-group
p <- select(medRep, t, starts_with("E")) %>%
    filter(t <= 100) %>%
    pivot_longer(!t, values_to = "counts", names_to = "var") %>%
    separate(var, c("var", "age", "lad"), sep = "_") %>%
    select(!var) %>%
    mutate(age = as.numeric(age), lad = as.numeric(lad)) %>%
    group_by(lad) %>%
    mutate(tot = sum(counts)) %>%
    ungroup() %>%
    filter(tot > 0) %>%
    select(!tot) %>%
    group_by(lad, age) %>%
    mutate(tot = cumsum(counts)) %>%
    group_by(lad, t) %>%
    mutate(tot = sum(tot)) %>%
    ungroup() %>%
    mutate(counts = ifelse(tot == 0, NA, counts)) %>%
    select(!tot)
max_p <- max(p$counts, na.rm = TRUE)
p <- inner_join(lad19, p, by = c("objectid" = "lad")) %>%
    ggplot() +
    geom_sf(aes(fill = counts), colour = NA) +
    facet_wrap(~age) +
    scale_fill_viridis_c(limits = c(0, max_p)) +
    theme_bw()

## add transitions
p <- p + transition_time(t) + ggtitle("Day = {frame_time}")
spatial_gif <- animate(p, nframes = 50, fps = 1, renderer = gifski_renderer())
anim_save("outputs/simsSpatialE.gif", spatial_gif)
