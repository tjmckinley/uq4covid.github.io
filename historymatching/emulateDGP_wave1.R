## load libraries
X11()
dev.off()
library(MASS)
library(tidyverse)
library(hmer)
library(mclust)
library(fields)
library(dgpsi)
init_py()

## read in emulator objects
source("emulator_fns.R")

## set parameter ranges
ranges <- list(
    beta_scale = c(0, 1),
    MD_scale = c(0.1, 1),
    MD_time_EE = c(0, 37),
    MD_time_EM = c(0, 37),
    MD_time_L = c(0, 37),
    MD_time_NE = c(0, 37),
    MD_time_NW = c(0, 37),
    MD_time_SE = c(0, 37),
    MD_time_SW = c(0, 37),
    MD_time_WM = c(0, 37),
    MD_time_YH = c(0, 37),
    nuA = c(0, 1),
    p_move = c(0, 1),
    R0 = c(2, 4.5),
    TE = c(0.1, 2),
    TI1 = c(2.8, 4.5),
    TI2 = c(0.0001, 0.5),
    TP = c(1.2, 3),
    alphaTH = c(-1.4, 1.8),
    etaTH = c(0.001, 0.08),
    alphaEP = c(-5, 0),
    alphaI1D = c(-20, 0),
    alphaHD = c(-20, 0),
    alphaI1H = c(-5, 0),
    eta = c(0, 0.05),
    etaI_scale = c(0.5, 2),
    etaH_scale = c(0.5, 2)
)

## save ranges
saveRDS(ranges, "ranges.rds")

## test it out with log-likelihood estimation
inputs <- readRDS("inputs.rds")
ll <- readRDS("ll.rds")

## build BL emulators
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

## separate inputs and outputs
training_ll <- training[, ncol(training)]
training <- training[, -ncol(training)]
test_ll <- test[, ncol(test)]
test <- test[, -ncol(test)]

## scale inputs to (0, 1)
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
actives <- m_gp$emulator_obj$kernel$length < 500 ## 500 is arbitrary - don't use blindly

m_gp <- gp(training[, actives], training_ll, name = "matern2.5", nugget_est = TRUE) 
m_gp$emulator_obj$kernel$length
plot(m_gp)

## if you need to run more times use code below
#actives[which(actives)] <- m_gp$emulator_obj$kernel$length < 500 # 500 is arbitrary - don't use blindly
#m_gp <- gp(training[, actives], training_ll, name = "matern2.5", nugget_est = TRUE) 
#m_gp$emulator_obj$kernel$length
#plot(m_gp)

## active variables
colnames(training)[actives]

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

## add test points into emulator
m_gp <- update(
    m_gp, 
    X = rbind(training[, actives], test[, actives]),
    Y = c(training_ll, test_ll),
    B = 1
)
write(m_gp, "m_gp.pkl")

## generate random seeds for reproducibility of
## imputations in future waves
seeds <- round(runif(30, 0, 50000000))
saveRDS(seeds, "seeds.rds")

## save variables for emulator
saveRDS(targets, "targets.rds")
saveRDS(0.95, "cutoff.rds")
saveRDS(actives, "actives.rds")

## CHECK FOR DOUBT POINTS AND SET VORONOI REGIONS IF REQUIRED

## set seeds for reproducibility
dgpsi::set_seed(seeds[1])
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

## build custom emulator
emulator <- list(ll = Proto_emulator$new(
        readRDS("ranges.rds"),
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
        range_list = readRDS("ranges.rds"),
        epsilon = 0.00001,
        actives = readRDS("actives.rds"),
        target = readRDS("targets.rds"),
        pcutoff = readRDS("cutoff.rds"),
        voronoi = readRDS("voronoi.rds")
    )
)

## check best point is in NROY
temp <- readRDS("inputs.rds")
temp_ll <- readRDS("ll.rds")
temp <- select(temp, !output)
stopifnot(max(-exp(log(-temp_ll))) == targets$ll$val)
new_set <- nth_implausible(emulator, temp, list(ll = list(val = 1, sd = 1)), n = 1)
temp_ll <- temp_ll[new_set <= 1]
temp <- temp[new_set <= 1, ]
stopifnot(max(-exp(log(-temp_ll))) == targets$ll$val)

## generate large sample from input space

## source dataTools
source("../inputs/dataTools.R")
library(lhs)

parRanges <- data.frame(
    parameter = c("R0", "TE", "TP", "TI1", "TI2", "nuA", "beta_scale", "p_move", "MD_scale", "MD_time_EM", "MD_time_EE", "MD_time_L", "MD_time_NE", "MD_time_NW", "MD_time_SE", "MD_time_SW", "MD_time_WM", "MD_time_YH"),
    lower = c(2, 0.1, 1.2, 2.8, 0.0001, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    upper = c(4.5, 2, 3, 4.5, 0.5, 1, 1, 1, 1, 37, 37, 37, 37, 37, 37, 37, 37, 37),
    stringsAsFactors = FALSE
) 

## generate LHS design
ndesign <- 1000
design <- randomLHS(ndesign, nrow(parRanges))
colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## convert to input space
inputs <- convertDesignToInput(design, parRanges, "zero_one")

## generate space-filling design for other parameters

## generate design points for hospital stay lengths
hospStaysInput <- FMMmaximin(
    hospStays, 
    ndesign + 20, 
    matrix(c(-Inf, Inf, 0, Inf), ncol = 2, byrow = TRUE)
) %>%
    as_tibble() %>%
    rename(alphaTH = x1, etaTH = x2)
    
## check against prior density restrictions
hospStaysInput <- hospStaysInput[dens(as.matrix(hospStaysInput), hospStays$modelName, hospStays$parameters, logarithm = TRUE) > hospThresh, ]
if(nrow(hospStaysInput) < ndesign) stop("Can't generate enough valid hospital points")
hospStaysInput <- hospStaysInput[1:ndesign, ]

## generate design points for other transition probabilities

## produces design points subject to constraints
pathwaysInput <- FMMmaximin(
    pathwaysMod, 
    ndesign + 20,
    matrix(c(rep(c(-20, 0), times = 4), 0, 1), ncol = 2, byrow = TRUE),
    pathwaysLimitFn,
    eta_scale = c(0.5, 2),
    ages = c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5)
) %>%
    as_tibble() %>%
    rename(alphaEP = x1, alphaI1D = x2, alphaHD = x3, alphaI1H = x4, eta = x5, etaI_scale = x6, etaH_scale = x7)
    
## check against prior density restrictions
pathwaysInput <- pathwaysInput[dens(as.matrix(pathwaysInput)[, -c(6, 7)], pathwaysMod$modelName, pathwaysMod$parameters, logarithm = TRUE) > pathwaysThresh, ]
pathwaysInput <- pathwaysInput[pathwaysInput[, 6] > 0.5 & pathwaysInput[, 6] < 2, ]
pathwaysInput <- pathwaysInput[pathwaysInput[, 7] > 0.5 & pathwaysInput[, 7] < 2, ]
if(nrow(pathwaysInput) < ndesign) stop("Can't generate enough valid pathways points")
pathwaysInput <- pathwaysInput[1:ndesign, ]

## bind to design
inputs <- cbind(inputs, hospStaysInput, pathwaysInput)

## add unique hash identifier
## (at the moment don't use "a0" type ensembleID, because MetaWards
## parses to dates)
inputs$output <- ensembleIDGen(ensembleID = "Ens1", nrow(inputs))

## convert input to disease
disease <- convertInputToDisease(inputs, C, N, S0, ages)
inputs <- semi_join(inputs, disease, by = "output")
    
## bind to design
inputs <- select(inputs, !!colnames(training))
saveRDS(inputs, "initial_samples.rds")

## calculate space removed
targets <- list(ll = list(val = 1, sd = 1))
sprem <- space_removal(emulator, targets, inputs, cutoff = 1)
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
write_csv(new_points, "inputsWave2.csv")

## visualise new points
waves <- list(readRDS("inputs.rds"))
waves[[length(waves) + 1]] <- new_points

## plot all waves design
md_time <- purrr::map(waves, ~dplyr::select(., starts_with("MD_time") | starts_with("MD_scale")))
ranges <- emulator$ll$ranges
ranges <- ranges[names(ranges) %in% colnames(md_time[[1]])]
p <- wave_points(md_time, input_names = names(ranges), p_size = 1, zero_in = FALSE, 
        wave_numbers = 1:length(waves)) +
    labs(fill = "Wave", title = NULL)
ggsave("new_design_MD_time.pdf", p, width = 10, height = 10)

waves <- purrr::map(waves, ~dplyr::select(., !starts_with("MD_time") & !starts_with("MD_scale")))
ranges <- emulator$ll$ranges
ranges <- ranges[names(ranges) %in% colnames(waves[[1]])]
p <- wave_points(waves, input_names = names(ranges), p_size = 1, zero_in = FALSE, 
        wave_numbers = 1:length(waves)) +
    labs(fill = "Wave", title = NULL)
ggsave("new_design.pdf", p, width = 10, height = 10)

## generate baseline set of points for the next design (and subsequent waves)
baseline_points <- generate_new_design(emulator, 1000, targets, cutoff = 1, method = "slice", plausible_set = baseline_points)
saveRDS(baseline_points, "baseline_samples.rds")

