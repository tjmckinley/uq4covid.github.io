#!/bin/bash

## load jaspy to access R
module load jasr

## read in jobs
readarray -t jobs < "job_lookup.txt"

## extract path name
jobname=${jobs[$1-1]}
echo $jobname

## run R script
cd ..
cmd="R CMD BATCH --no-restore --no-save --slave '--args $2 ${jobname} $3' JASMINcode/plotxSum.R plot$2Sum_${jobname}.Rout"
eval $cmd

cmd="mv plot$2Sum_${jobname}.Rout wave$2"
eval $cmd

