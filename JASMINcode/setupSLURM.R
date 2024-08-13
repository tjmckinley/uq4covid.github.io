## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) >= 3)
        wave <- args[1]
        runCode <- args[2]
        outputs <- args[3]
        if(length(args) == 4) {
            time <- args[4]
        } else {
            time <- "00:35:00"
        }
        if(length(args) == 5) {
            inds <- args[5]
        } else {
            inds <- NA
        }
    } else {
        stop("No arguments")
    }
} else {
    ## set name of directory to save outputs
    wave <- "1Feb"
    runCode <- "runDesign"
    outputs <- "Feb"
    time <- "00:35:00"
    inds <- NA
    #runCode <- "runPlotSum"
    #runCode <- "runPlotAgg"
}
if(!runCode %in% c("runDesign", "runForecasts", "runPlotSum", "runPlotAgg")) stop("Incorrect 'runCode'")

## append outputs
outputs <- paste0("outputs", outputs)

## read in input file
pars <- readRDS(paste0("../wave", wave, "/disease.rds"))

## add additional inds argument if specified
if(!is.na(inds)) {
    outputs <- paste(outputs, inds)
}

## read run code
code <- readLines("submit_job_template.sbatch")
code <- gsub("FILEDIR", wave, code)
code <- gsub("RUNCODE", runCode, code)
code <- gsub("OUTPUTS", outputs, code)
code <- gsub("TIME", time, code)

if(runCode != "runPlotAgg") {
    
    ## update run code
    code <- gsub("RANGES", paste0("1-", nrow(pars)), code)

    ## write csv to query
    write.table(data.frame(job = 1:nrow(pars)), paste0("job_lookup_wave", wave, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
    
} else {
    
    ## load representative run in
    runs <- readRDS(paste0("../wave", wave, "/plotSum_1_natFull.rds"))

    ## extract days / ages / classes
    days <- unique(runs$t)
    
    ## expand grid
    jobs <- data.frame(t = days)
    
    ## update run code
    code <- gsub("RANGES", paste0("1-", nrow(jobs)), code)

    ## write csv to query
    write.table(jobs, paste0("job_lookup_wave", wave, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
}   

## write run code
writeLines(code, paste0("submit_job_wave", wave, ".sbatch"))
    
print("All done.")

