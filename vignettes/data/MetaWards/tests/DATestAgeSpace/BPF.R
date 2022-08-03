## log-sum-exp function
log_sum_exp <- function(x, mn = FALSE) {
    maxx <- max(x)
    y <- maxx + log(sum(exp(x - maxx)))
    if(mn) y <- y - log(length(x))
    y
}

## run model for each set of design points in pars
## pars: matrix / data frame of parameters in order: 
##       nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD
##       each parameter except nu/nuA have entries for
##       each age-group e.g. nu, nuA, pE1, pE2, pE3 etc.
## C1:    contact matrix for mixing between age-classes
## C2:    contact matrix for mixing between age-classes post lockdown
## lockdown_day: the day of first lockdown
## cumDeath_lad: matrix of form D_1, D_2, ..., D_L along
##               columns with cumulative deaths per lad
## cumDeath_age_region: matrix of form D_a1_1, D_a1_2, ..., D_AR along
##               columns with cumulative deaths per age and region
## hosp_nhsregion: matrix of form H_1, H_2, ..., H_R along
##               columns with hospital counts per nhsregion
## cumHospAd_age_nhsregion: matrix of form H_a1_1, H_a1_2, ..., H_AR along
##               columns with cumulative admissions per age and nhsregion
## lookup: a lookup table mapping lads to different spatial hierarchies
##               of form: LAD, death ID, region ID, NHS region ID
## age_lookup: a lookup table of form: death age, nhs age
## u1_moves: matrix of commuter data, with columns:
##            home LAD, work LAD
## u1: 3D array of commuter data, with dimensions:
##            nclasses x nages x nrow(u1_moves)
## u2_moves: matrix of commuter data, with columns:
##            home LAD, play LAD, prob of moving
## u2: 3D array of player data, with dimensions:
##            nclasses x nages x nlads
## ndays: the number of days to fit to
## npart: the number of particles
## niter: number of Metropolis-Hastings steps to counter particle impoverishment
## obsScale: scaling parameter for Skellam observation process (see code)
## a1, a2, b: parameters for Skellam observation process
## saveAll: a logical specifying whether to return all states (if FALSE then returns just observed states))
## writeExt: a logical denoting whether to save particles externally or not
## PF:      a logical denoting whether to run a particle filter, or just simulate from the model
## ncores:  the number of cores for OpenMP parallelisation (if NA then defaults to all available cores)
## parEnsemble: decides whether to parallelise across or within ensemble

BPF <- function(pars, C1, C2, lockdown_day, cumDeath_lad, cumDeath_age_region, hosp_nhsregion, cumHospAd_age_nhsregion,
        lookup, age_lookup, u1_moves, u1, u2_moves, u2, ndays, npart = 10, niter = 0, 
        a1 = 0.01, a2 = 0.2, b = 0.1, a_dis = 0.05, b_dis = 0.05, 
        sigma2_lad = 1, sigma2_age_region = 1, sigma2_nhsregion = 1, sigma2_age_nhsregion = 1,
        saveAll = NA, writeExt = FALSE, outputName = "saveOut", PF = TRUE, ncores = NA, parEnsemble = FALSE) {
               
    ## set default for saveAll if PF = FALSE
    if(!PF & is.na(saveAll)) saveAll <- TRUE
    
    ## convert indicators
    if(is.na(saveAll)) {
        saveAllint <- 0
    } else {
        saveAllint <- ifelse(saveAll, 2, 1)
    }
    PFint <- ifelse(PF, 1, 0)
    writeExtint <- ifelse(writeExt, 1, 0)
    
    ## check number of requested cores
    if(is.na(ncores)) {
        ncores <- parallel::detectCores()
    }
    
    ## set up dummy "data" if just simulation required
    if(!PF) {
        cumDeath_lad <- matrix(NA, 1, 1)
        cumDeath_age_region <- matrix(NA, 1, 1)
        hosp_nhsregion <- matrix(NA, 1, 1)
        cumHospAd_age_nhsregion <- matrix(NA, 1, 1)
    }
    
    ## check no missing lads
    if(
        !identical(sort(unique(u1_moves[, 1])), sort(unique(u1_moves[, 1]))) |
        !identical(sort(unique(u1_moves[, 1])), sort(unique(u1_moves[, 2]))) |
        !identical(sort(unique(u1_moves[, 1])), sort(unique(u2_moves[, 1]))) |
        !identical(sort(unique(u1_moves[, 1])), sort(unique(u2_moves[, 2])))
    ) {
        stop("Can't match LADS")
    }
    if(any((1:length(unique(u1_moves[, 1])) - sort(unique(u1_moves[, 1]))) != 0)) stop("Can't match all LADs")
    if(dim(u2)[3] != length(unique(u2_moves[, 1]))) stop("u2 does not have entries for all LADS")
    if(!identical(dim(u1)[1:2], dim(u2)[1:2])) stop("u1 and u2 have different numbers of ages / classes")
    
    ## check probabilities
    if(!all((tapply(u2_moves[, 3], u2_moves[, 1], sum) - 1) < 1e-15)) stop("u2 probs don't sum to one")
    
    ## create dummy cohorts for players and bind with workers
    temp <- array(0, c(dim(u1)[1], dim(u1)[2], nrow(u2_moves)))
    u1 <- abind(u1, temp, along = 3)
    u1_moves <- rbind(cbind(u1_moves, rep(NA, nrow(u1_moves))), u2_moves)
    inds <- order(u1_moves[, 1], -u1_moves[, 3])
    u1_moves <- u1_moves[inds, ]
    u1 <- u1[, , inds, drop = FALSE]
    
    ## convert play probs to conditional probs
    playprobs <- tapply(u1_moves[, 3], u1_moves[, 1], function(y) {
        x <- y[!is.na(y)]
        tp <- x[1]
        for(i in 2:length(x)) {
            x[i:length(x)] <- x[i:length(x)] / (1 - tp)
            tp <- x[i]
            ## deal with rounding errors
            if(sum(x[i:length(x)]) != 1) {
#                print(sum(x[i:length(x)]))
                x[i:length(x)] <- x[i:length(x)] / sum(x[i:length(x)])
            }
        }
        x[length(x)] <- 1
        y[!is.na(y)] <- x
        y
    })
    playprobs <- do.call("c", playprobs)
    
    ## generate number of cohorts
    ncohorts1 <- tapply(u1_moves[, 1], u1_moves[, 1], length)
    ncohorts1 <- c(0, cumsum(ncohorts1))
    names(ncohorts1) <- NULL
    
    ## generate number of play cohorts
    ncohorts2 <- tapply(u1_moves[, 3], u1_moves[, 1], function(x) sum(!is.na(x)))
    names(ncohorts2) <- NULL
    
    ## remove extraneous column
    u1_moves <- u1_moves[, -3]   
    
    ## set up output folder
    if(writeExt) {
        if(dir.exists(outputName)) system(paste0("rm -rf ", outputName))
        dir.create(outputName, recursive = TRUE)
    }
    
    ## set up auxiliary parameters guiding parallelisation
    if(parEnsemble) {
        ncoresEns <- ncores
        ncores <- 1
    } else {
        ncoresEns <- 1
    }
    
    ## run particle filter for each set of inputs
    runs <- mclapply(1:nrow(pars), function(k, pars, C1, C2, lockdown_day, u1_moves, ncohorts1, u1, u2, playprobs, ncohorts2, npart, niter, ndays, cumDeath_lad, cumDeath_age_region, hosp_nhsregion, cumHospAd_age_nhsregion, lookup, age_lookup, a1, a2, b, a_dis, b_dis, sigma2_lad, sigma2_age_region, sigma2_nhsregion, sigma2_age_nhsregion, saveAll, writeExt, outputName, PF, ncores) {
    
        ## set up lookups and numbers of regions
        lookup <- as.matrix(lookup)
        age_lookup <- as.matrix(age_lookup)
        ndeathlads <- max(lookup[, 2], na.rm = TRUE)
        nregions <- max(lookup[, 3], na.rm = TRUE)
        nnhsages <- max(age_lookup[, 2], na.rm = TRUE)
        nnhsregions <- max(lookup[, 4], na.rm = TRUE)
        lookup[is.na(lookup)] <- -1
        
        if(PF == 1) {
            ## reformat observations
            deathInc_lad <- mutate(cumDeath_lad, across(!t, ~. - lag(., default = 0))) %>%
                select(!t) %>%
                as.matrix()
            deathInc_age_region <- mutate(cumDeath_age_region, across(!t, ~. - lag(., default = 0))) %>%
                select(!t) %>%
                as.matrix()
            hosp_nhsregion <- as.matrix(hosp_nhsregion)[, -1]
            hospInc_age_nhsregion <- mutate(cumHospAd_age_nhsregion, across(!t, ~. - lag(., default = 0))) %>%
                select(!t) %>%
                as.matrix()
        } else {
            ## set dummies if required
            deathInc_lad <- matrix(0, 1, 1)
            deathInc_age_region <- matrix(0, 1, 1)
            hosp_nhsregion <- matrix(0, 1, 1)
            hospInc_age_nhsregion <- matrix(0, 1, 1)
        }
        
        ## set pars
        pars <- unlist(pars[k, ])
        
        ## extract number of stages, age classes and lads
        nclasses <- dim(u1)[1]
        nages <- dim(u1)[2]
        nlads <- max(u1_moves[, 1])
    
        ## do garbage collection (seems to solve allocation issue)
        gc()
        
        ## if just simulating
        if(PF == 0) {
            ## check return outputs
            if(saveAll == 0) stop("Must set 'saveAll' to something if not running a PF")
            
            ## run particle filter
            particles <- BPF_cpp(pars, C1, C2, lockdown_day, deathInc_lad, deathInc_age_region, hosp_nhsregion, 
                hospInc_age_nhsregion, lookup, age_lookup, nclasses, nages, nlads, ndeathlads, 
                nregions, nnhsages, nnhsregions, u1_moves, ncohorts1, u1, u2, playprobs, 
                ncohorts2, ndays, npart, niter, a1, a2, b, a_dis, b_dis, sigma2_lad, 
                sigma2_age_region, sigma2_nhsregion, sigma2_age_nhsregion, saveAll, writeExt, 
                outputName, PF, ncores)
            return(particles)
        }
        
        ## run particle filter
        particles <- BPF_cpp(pars, C1, C2, lockdown_day, deathInc_lad, deathInc_age_region, hosp_nhsregion, 
                hospInc_age_nhsregion, lookup, age_lookup, nclasses, nages, nlads, ndeathlads, 
                nregions, nnhsages, nnhsregions, u1_moves, ncohorts1, u1, u2, playprobs, 
                ncohorts2, ndays, npart, niter, a1, a2, b, a_dis, b_dis, sigma2_lad, 
                sigma2_age_region, sigma2_nhsregion, sigma2_age_nhsregion, saveAll, writeExt, 
                outputName, PF, ncores)
        
        if(saveAll != 0 & writeExt == 0) {
            return(list(ll = particles$ll, particles = particles$particles))
        } else {
            return(list(ll = particles$ll))
        }
    }, pars = pars, C1 = C1, C2 = C2, lockdown_day = lockdown_day, u1_moves = u1_moves, ncohorts1 = ncohorts1, u1 = u1, u2 = u2, playprobs = playprobs, ncohorts2 = ncohorts2, npart = npart, niter = niter, ndays = ndays, cumDeath_lad = cumDeath_lad, cumDeath_age_region = cumDeath_age_region, hosp_nhsregion = hosp_nhsregion, cumHospAd_age_nhsregion = cumHospAd_age_nhsregion, lookup = lookup, age_lookup = age_lookup, a1 = a1, a2 = a2, b = b, a_dis = a_dis, b_dis = b_dis, sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion, saveAll = saveAllint, writeExt = writeExtint, outputName = outputName, PF = PFint, ncores = ncores, mc.cores = ncoresEns)
    if(!is.na(saveAll)) {
        if(!writeExt) {
            if(PF) {
                ll <- map(runs, "ll")
                runs <- map(function(runs, ndays, npart) {
                        lads <- map(runs, "lads")
                        age_region <- map(runs, "age_region")
                        nhsregion <- map(runs, "nhsregion")
                        age_nhsregion <- map(runs, "age_nhsregion")
                        xlads <- list()
                        xage_region <- list()
                        xnhsregion <- list()
                        xage_nhsregion <- list()
                        for(i in 1:(ndays + 1)) {
                            xlads[[i]] <- lads[(i - 1) * npart + 1:npart]
                            xage_region[[i]] <- age_region[(i - 1) * npart + 1:npart]
                            xnhsregion[[i]] <- nhsregion[(i - 1) * npart + 1:npart]
                            xage_nhsregion[[i]] <- age_nhsregion[(i - 1) * npart + 1:npart]
                        }
                        list(lads = xlads, age_region = xage_region, nhsregion = xnhsregion, age_nhsregion = xage_nhsregion)
                    }, ndays = ndays, npart = npart)
                ll <- do.call("c", ll)
                return(list(ll = ll, particles = runs))
            } else {
                runs <- map(runs, "particles") %>%
                    map(function(runs, ndays, npart, saveAll) {
                        lads <- map(runs, "lads")
                        age_region <- map(runs, "age_region")
                        nhsregion <- map(runs, "nhsregion")
                        age_nhsregion <- map(runs, "age_nhsregion")
                        xlads <- list()
                        xage_region <- list()
                        xnhsregion <- list()
                        xage_nhsregion <- list()
                        if(saveAll) {
                            full <- map(runs, "full")
                            xfull <- list()
                        }
                        for(i in 1:(ndays + 1)) {
                            xlads[[i]] <- lads[(i - 1) * npart + 1:npart]
                            xage_region[[i]] <- age_region[(i - 1) * npart + 1:npart]
                            xnhsregion[[i]] <- nhsregion[(i - 1) * npart + 1:npart]
                            xage_nhsregion[[i]] <- age_nhsregion[(i - 1) * npart + 1:npart]
                            if(saveAll) {
                                xfull[[i]] <- full[(i - 1) * npart + 1:npart]
                            }
                        }
                        if(saveAll) {
                            return(list(full = xfull, lads = xlads, age_region = xage_region, nhsregion = xnhsregion, age_nhsregion = xage_nhsregion))
                        } else {
                            return(list(lads = xlads, age_region = xage_region, nhsregion = xnhsregion, age_nhsregion = xage_nhsregion))
                        }
                    }, ndays = ndays, npart = npart, saveAll = saveAll)
                return(list(particles = runs))
            }
        } else {
            if(PF) {
                ll <- map(runs, "ll")
                ll <- do.call("c", ll)
                return(list(ll = ll))
            } else {
                return(NULL)
            }
        }       
    } else {
        ll <- map(runs, "ll")
        ll <- do.call("c", ll)
        return(ll)
    }
}
