## load libraries
X11()
dev.off()
library(MASS)
library(tidyverse)
library(hmer)  ## experimental branch needed
library(mclust)
library(fields)
library(dgpsi) ## needs development version for set_seed()
init_py()

## read in emulator objects
source("emulator_fns.R")

## set vector of previous waves
wave_name <- "Feb"
wave_nos <- 1:2
prev_waves <- paste0("../wave", wave_nos, wave_name)

## read in seeds
seeds <- readRDS(paste0(prev_waves[1], "/seeds.rds"))

## build custom emulator for previous waves
emulator <- list()
for(i in 1:(length(prev_waves) - 1)) {

    ## set seeds for reproducibility
    m_gp <- read(paste0(prev_waves[i], "/m_gp.pkl"))
    dgpsi::set_seed(seeds[i])
    m_gp <- set_imp(m_gp, B = 30)
    
    ## set up emulator
    emulator[[i]] <- list(ll = Proto_emulator$new(
            readRDS(paste0(prev_waves[i], "/ranges.rds")),
            "ll",
            pred_func,
            pred_var_func,
            implausibility,
            print_func = function(em) {
                print(summary(em))
            },
            em = m_gp,
            hospStays = hospStays,
            hospThresh = hospThresh,
            pathwaysMod = pathwaysMod,
            pathwaysThresh = pathwaysThresh,
            ages = ages,
            pathwaysLimitFn = pathwaysLimitFn,
            NGM = NGM, 
            N = N, 
            S0 = S0, 
            C = C,
            range_list = readRDS(paste0(prev_waves[i], "/ranges.rds")),
            epsilon = 0.00001,
            actives = readRDS(paste0(prev_waves[i], "/actives.rds")),
            target = readRDS(paste0(prev_waves[i], "/targets.rds")),
            pcutoff = readRDS(paste0(prev_waves[i], "/cutoff.rds")),
            voronoi = readRDS(paste0(prev_waves[i], "/voronoi.rds"))
        )
    )
}

## read in original training data
inputs <- list()
ll <- list()
for(i in 1:(length(prev_waves) - 1)) {
    inputs[[i]] <- readRDS(paste0(prev_waves[i], "/inputs.rds"))
    ll[[i]] <- readRDS(paste0(prev_waves[i], "/ll.rds"))
}
inputs <- reduce(inputs, rbind)
ll <- reduce(ll, c)
inputs <- select(inputs, !output)

## check against implausibility at previous waves
new_set <- nth_implausible(emulator, inputs, list(ll = list(val = 1, sd = 1)), n = 1)
prev_points <- inputs[new_set <= 1, ]
prev_ll <- ll[new_set <= 1]

## check best point is in NROY
temp <- map_dbl(prev_waves[1:(length(prev_waves) - 1)], ~{
    readRDS(paste0(., "/targets.rds"))$ll$val
})
if(length(temp) > 1) {
    stopifnot(all(diff(temp) >= 0))
    temp <- max(temp)
}   
stopifnot(max(-exp(log(-prev_ll))) == temp)

## read in current runs
inputs <- readRDS("inputs.rds")
ll <- readRDS("ll.rds")

## extract training and test data
training <- mutate(inputs, ll = ll) %>%
    slice(1:200) %>%
    select(!output) %>%
    mutate(ll = log(-ll)) %>%
    as.matrix()
test <- mutate(inputs, ll = ll) %>%
    slice(-c(1:200)) %>%
    select(!output) %>%
    mutate(ll = log(-ll)) %>%
    as.matrix()
    
## bind to previous runs
if(length(prev_ll) > 0) {
    training <- rbind(training, as.matrix(cbind(prev_points, ll = log(-prev_ll))))
}

## separate inputs and outputs
training_ll <- training[, ncol(training)]
training <- training[, -ncol(training)]
test_ll <- test[, ncol(test)]
test <- test[, -ncol(test)]

## scale inputs to (0, 1)
ranges <- readRDS(paste0(prev_waves[1], "/ranges.rds"))
for(i in 1:ncol(training)) {
    temp_ranges <- ranges[[colnames(training)[i]]]
    training[, i] <- (training[, i] - temp_ranges[1]) / diff(temp_ranges)
    test[, i] <- (test[, i] - temp_ranges[1]) / diff(temp_ranges)
}
stopifnot(all(apply(training, 2, function(x) all(x > 0 & x < 1))))
stopifnot(all(apply(test, 2, function(x) all(x > 0 & x < 1))))

## set target
targets <- list(ll = list(val = max(-exp(c(training_ll, test_ll))), sigma = 1))

## find active variables using ARD in a GP
m_gp <- gp(training, training_ll, name = "matern2.5", nugget_est = TRUE)
plot(m_gp)
m_gp$emulator_obj$kernel$length ## large values inactive
actives <- m_gp$emulator_obj$kernel$length < 1000 ## 1000 is arbitrary - don't use blindly

m_gp <- gp(training[, actives], training_ll, name = "matern2.5", nugget_est = TRUE) 
m_gp$emulator_obj$kernel$length
plot(m_gp)

## if you need to run more times use code below
#actives[which(actives)] <- m_gp$emulator_obj$kernel$length < 1000
#m_gp <- gp(training[, actives], training_ll, name = "matern2.5", nugget_est = TRUE) 
#m_gp$emulator_obj$kernel$length
#plot(m_gp)

## active variables
colnames(training)[actives]

## save emulator
saveRDS(targets, "targets.rds")
saveRDS(actives, "actives.rds")

## fit two-layer deep GP
m_gp <- dgp(
    training[, actives],
    training_ll, 
    name = "matern2.5",
    depth = 2,
    likelihood = "Hetero",
    N = 1500,
    B = 30
)
pdf("traces.pdf")
trace_plot(m_gp, layer = 1, node = 1)
trace_plot(m_gp, layer = 1, node = 2)
dev.off()

## save object
write(m_gp, "m_gp.pkl")

### continue training for longer if needed
#m_gp <- continue(m_gp, N = 1500)
#pdf("traces.pdf")
#trace_plot(m_gp, layer = 1, node = 1)
#trace_plot(m_gp, layer = 1, node = 2)
#dev.off()

### save object
#write(m_gp, "m_gp.pkl")

## validation plots
summary(m_gp)
m_gp <- set_imp(m_gp, B = 30)
m_gp <- validate(m_gp)
m_gp <- validate(
    m_gp, 
    test[, actives],
    test_ll
)
pdf("validation.pdf")
plot(m_gp)
plot(
    m_gp, 
    test[, actives],
    test_ll
)
dev.off()

## add test points into emulator object
m_gp <- update(
    m_gp, 
    X = rbind(m_gp$data$X, test[, actives]),
    Y = c(m_gp$data$Y[, 1], test_ll),
    B = 1
)

## save object
write(m_gp, "m_gp.pkl")

## set cutoff
saveRDS(0.95, "cutoff.rds")

## CHECK FOR DOUBT POINTS AND SET VORONOI REGIONS IF REQUIRED

## set seeds for reproducibility
dgpsi::set_seed(seeds[length(prev_waves)])
m_gp <- set_imp(m_gp, B = 30)

## extract failed points
imp_mean <- predict(m_gp, m_gp$data$X)$results$mean
imp_var <- predict(m_gp, m_gp$data$X)$results$var
XD <- pnorm(log(-log(0.00001) - targets$ll$val), mean = imp_mean, sd = sqrt(imp_var))
XF <- which((1 - XD) > readRDS("cutoff.rds"))

## calculate doubt points
voronoi <- NULL
if(length(XF) > 0) {
    XD <- XF[-exp(m_gp$data$Y[XF, 1]) - targets$ll$val > log(0.00001)]

    if(length(XD) > 0) {
        ## extract lengthscale
        lscale <- m_gp$specs$layer1$node1$lengthscales

        ## augment space of doubt points if necessary
        temp <- matern(
            m_gp$data$X[-XD, , drop = FALSE],
            m_gp$data$X,
            2.5,
            lscale
        )
        temp <- which(apply(temp, 1, function(x, XD) {
            inds <- which(x == max(x, na.rm = TRUE))
            any(inds %in% XD)
        }, XD = XD))
        XD <- sort(unique(c(XD, temp)))
        if(length(XD) > 0) {
            voronoi <- list(kernel = matern, XD = XD)
        }
    }
}   

## save Voronoi region information
saveRDS(voronoi, "voronoi.rds")

## save ranges information
ranges <- readRDS(paste0(prev_waves[1], "/ranges.rds"))
saveRDS(ranges, "ranges.rds")

## build custom emulator for current wave
i <- length(prev_waves)

## set seeds for reproducibility
dgpsi::set_seed(seeds[i])
m_gp <- set_imp(m_gp, B = 30)

## set up emulator
emulator[[i]] <- list(ll = Proto_emulator$new(
        readRDS(paste0(prev_waves[i], "/ranges.rds")),
        "ll",
        pred_func,
        pred_var_func,
        implausibility,
        print_func = function(em) {
            print(summary(em))
        },
        em = m_gp,
        hospStays = hospStays,
        hospThresh = hospThresh,
        pathwaysMod = pathwaysMod,
        pathwaysThresh = pathwaysThresh,
        ages = ages,
        pathwaysLimitFn = pathwaysLimitFn,
        NGM = NGM, 
        N = N, 
        S0 = S0, 
        C = C,
        range_list = readRDS(paste0(prev_waves[i], "/ranges.rds")),
        epsilon = 0.00001,
        actives = readRDS(paste0(prev_waves[i], "/actives.rds")),
        target = readRDS(paste0(prev_waves[i], "/targets.rds")),
        pcutoff = readRDS(paste0(prev_waves[i], "/cutoff.rds")),
        voronoi = readRDS(paste0(prev_waves[i], "/voronoi.rds"))
    )
)

## check best point is in NROY
temp <- list()
temp_ll <- list()
for(i in 1:length(prev_waves)) {
    temp[[i]] <- readRDS(paste0(prev_waves[i], "/inputs.rds"))
    temp_ll[[i]] <- readRDS(paste0(prev_waves[i], "/ll.rds"))
}
temp <- reduce(temp, rbind)
temp_ll <- reduce(temp_ll, c)
temp <- select(temp, !output)
stopifnot(max(-exp(log(-temp_ll))) == targets$ll$val)
new_set <- nth_implausible(emulator, temp, list(ll = list(val = 1, sd = 1)), n = 1)
temp_ll <- temp_ll[new_set <= 1]
temp <- temp[new_set <= 1, ]
stopifnot(max(-exp(log(-temp_ll))) == targets$ll$val)

## read in baseline samples from previous non-implausible space
targets <- list(ll = list(val = 1, sd = 1))
inputs <- readRDS(paste0(prev_waves[length(prev_waves) - 1], "/baseline_samples.rds"))

## calculate space removed
sprem <- space_removal(emulator[[length(emulator)]], targets, inputs, cutoff = 1)
sprem
saveRDS(sprem, "space_removed.rds")

## extract previous design and update to extract which of these are still non-implausible
## with new emulator / targets
new_set <- nth_implausible(emulator, inputs, targets, n = 1)
inputs <- inputs[new_set <= 1, ]

## generate baseline set of points for the next design (and subsequent waves)
baseline_points <- generate_new_design(emulator, 1000, targets, cutoff = 1, method = "slice", plausible_set = inputs)
saveRDS(baseline_points, "baseline_samples.rds")

## generate new points using maximin
new_set <- maximin(baseline_points, 200, readRDS("ranges.rds"))
new_points <- baseline_points[new_set, ]
baseline_points <- baseline_points[-new_set, ]
new_set <- maximin(baseline_points, 50, readRDS("ranges.rds"))
new_points <- rbind(new_points, baseline_points[new_set, ])
rownames(new_points) <- 1:nrow(new_points)
baseline_points <- baseline_points[-new_set, ]

## write new runs
write_csv(new_points, paste0("inputsWave", max(wave_nos) + 1, ".csv"))

## visualise new points
waves <- purrr::map(prev_waves, ~readRDS(paste0(., "/inputs.rds")))
waves[[length(waves) + 1]] <- new_points

## plot all waves design
md_time <- purrr::map(waves, ~dplyr::select(., starts_with("MD_time") | starts_with("MD_scale")))
ranges <- emulator[[1]]$ll$ranges
ranges <- ranges[names(ranges) %in% colnames(md_time[[1]])]
p <- wave_points(md_time, input_names = names(ranges), p_size = 1, zero_in = FALSE, 
        wave_numbers = 1:length(waves)) +
    labs(fill = "Wave", title = NULL)
ggsave("new_design_MD_time.pdf", p, width = 10, height = 10)

waves <- purrr::map(waves, ~dplyr::select(., !starts_with("MD_time") & !starts_with("MD_scale")))
ranges <- emulator[[1]]$ll$ranges
ranges <- ranges[names(ranges) %in% colnames(waves[[1]])]
p <- wave_points(waves, input_names = names(ranges), p_size = 1, zero_in = FALSE, 
        wave_numbers = 1:length(waves)) +
    labs(fill = "Wave", title = NULL)
ggsave("new_design.pdf", p, width = 10, height = 10)

## generate baseline set of points for the next wave
baseline_points <- generate_new_design(emulator, 1000, targets, cutoff = 1, method = "slice", plausible_set = baseline_points)
saveRDS(baseline_points, "baseline_samples.rds")

