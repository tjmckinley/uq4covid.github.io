## load libraries
library(tidyverse)
library(patchwork)

## set wave, number of days and lads
## (If is.na(lads) then plots aggregated counts
wave <- 1
ndays <- 100

## read in ensemble summaries
runs_full <- readRDS(paste0("wave", wave, "/sumEns.rds"))

## read in simulated data and generate incidence curves
data <- readRDS("outputs/disSims.rds")

## set up plot data
plot_data <- pivot_longer(filter(data, t <= ndays), !t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(obs = grepl("obs", var)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    mutate(var = gsub("obs", "", var))
    
## NATIONAL LEVEL PLOT
    
## aggregate data
p <- filter(plot_data, !obs) %>%
    group_by(t, var, age) %>%
    summarise(n = sum(n), .groups = "drop")

## aggregate observed data   
pobs <- filter(plot_data, obs) %>%
    group_by(t, var, age) %>%
    summarise(n = sum(n), .groups = "drop")
    
## split simulations into truth and observed
runs_obs <- filter(runs_full, lad == "agg" & var %in% paste0(unique(pobs$var), "obs")) %>%
    mutate(var = gsub("obs", "", var))
    
runs <- filter(runs_full, lad == "agg" & !var %in% paste0(unique(pobs$var), "obs"))

## produce particle trajectories plot
p1 <- list()
p1[[1]] <- ggplot(runs, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, colour = NA) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.3, colour = NA) +
    geom_line(aes(y = Median)) +
    geom_line(aes(y = n), data = p, linetype = "dashed") +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Truth")
p1[[2]] <- ggplot(runs_obs, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, colour = NA) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.3, colour = NA) +
    geom_line(aes(y = Median)) +
    geom_line(aes(y = n), data = pobs, linetype = "dashed") +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed")
p1 <- wrap_plots(p1, nrow = 2, heights = c(0.8, 0.2))
ggsave(paste0("wave", wave, "/sumEns.pdf"), p1, width = 10, height = 10)

## LAD LEVEL PLOT

lads <- unique(runs_full$lad)
lads <- lads[-grep("agg", lads)]
plot_data <- filter(plot_data, LAD %in% lads)
    
## aggregate data
p <- filter(plot_data, !obs)

## aggregate observed data   
pobs <- filter(plot_data, obs)
    
## split simulations into truth and observed
runs_obs <- filter(runs_full, lad != "agg" & var %in% paste0(unique(pobs$var), "obs")) %>%
    mutate(var = gsub("obs", "", var)) %>%
    mutate(LAD = lad)
    
runs <- filter(runs_full, lad != "agg" & !var %in% paste0(unique(pobs$var), "obs")) %>%
    mutate(LAD = lad)

## plot of all states against simulations
p1 <- list()
p1[[1]] <- ggplot(runs, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = LAD), alpha = 0.3, colour = NA) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = LAD), alpha = 0.3, colour = NA) +
    geom_line(aes(y = Median, colour = LAD)) +
    geom_line(aes(y = n, colour = LAD), data = p, linetype = "dashed") +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Truth")
p1[[2]] <- ggplot(runs_obs, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = LAD), alpha = 0.3, colour = NA) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = LAD), alpha = 0.3, colour = NA) +
    geom_line(aes(y = Median, colour = LAD)) +
    geom_line(aes(y = n, colour = LAD), data = pobs, linetype = "dashed") +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed")
p1 <- wrap_plots(p1, nrow = 2, heights = c(0.8, 0.2)) +
    plot_layout(guides = "collect")
ggsave(paste0("wave", wave, "/sumEnsLADs.pdf"), p1, width = 10, height = 10)


