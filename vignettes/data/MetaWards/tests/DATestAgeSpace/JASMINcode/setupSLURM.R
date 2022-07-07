## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) >= 2)
        wave <- as.numeric(args[1])
        runCode <- args[2]
        if(length(args) > 2) {
            time <- args[3]
        } else {
            time <- "00:35:00"
        }
    } else {
        stop("No arguments")
    }
} else {
    ## set name of directory to save outputs
    wave <- 1
    runCode <- "runDesign"
    time <- "00:35:00"
    #runCode <- "runPlotSum"
    #runCode <- "runPlotAgg"
}
if(!runCode %in% c("runDesign", "runPlotSum", "runPlotAgg")) stop("Incorrect 'runCode'")

## read in input file
pars <- readRDS(paste0("../wave", wave, "/disease.rds"))

## read run code
code <- readLines("submit_job_template.sbatch")
code <- gsub("FILEDIR", wave, code)
code <- gsub("RUNCODE", runCode, code)
code <- gsub("TIME", time, code)

if(runCode != "runPlotAgg") {
    
    ## update run code
    code <- gsub("RANGES", paste0("1-", nrow(pars)), code)

    ## write csv to query
    write.table(data.frame(job = 1:nrow(pars)), "job_lookup.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
} else {
    
    ## load representative run in
    runs <- readRDS(paste0("../wave", wave, "/plotSum_1.rds"))

    ## extract days / ages / classes
    days <- unique(runs$time)
    
    ## expand grid
    jobs <- data.frame(t = days)
    
    ## update run code
    code <- gsub("RANGES", paste0("1-", nrow(jobs)), code)

    ## write csv to query
    write.table(jobs, "job_lookup.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
}   

## write run code
writeLines(code, "submit_job.sbatch")
    
print("All done.")

