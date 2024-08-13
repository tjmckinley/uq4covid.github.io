#!/bin/bash

## check correct number of arguments
if [ $# -ne 3 ] && [ $# -ne 4 ]; then
    echo "Command requires 3 or 4 arguments."
    exit 1
fi

## load jaspy to access R
module load jasr

## read in jobs
readarray -t jobs < "job_lookup_wave$2.txt"

## extract path name
jobname=${jobs[$1-1]}
echo $jobname

## change to main directory
cd ..

## run R script
## generate command with correct arguments
if [ $# -eq 3 ]; then
    cmd="R CMD BATCH --no-restore --no-save --slave '--args $2 $3 ${jobname}' JASMINcode/plotxAgg_full.R plot$2Agg_full_T${jobname}.Rout"
else
    cmd="R CMD BATCH --no-restore --no-save --slave '--args $2 $3 ${jobname} $4' JASMINcode/plotxAgg_full.R plot$2Agg_full_T${jobname}.Rout"
fi

## evaluate command
eval $cmd

cmd="mv plot$2Agg_full_T${jobname}.Rout wave$2"
eval $cmd

exit 0

