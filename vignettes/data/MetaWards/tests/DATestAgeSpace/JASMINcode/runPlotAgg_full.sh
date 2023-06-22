#!/bin/bash

## load jaspy to access R
module load jasr

## read in jobs
readarray -t jobs < "job_lookup_wave$2.txt"

## extract path name
jobname=${jobs[$1-1]}
echo $jobname

## run R script
cd ..
cmd="R CMD BATCH --no-restore --no-save --slave '--args $2 $3 ${jobname}' JASMINcode/plotxAgg_full.R plot$2Agg_full_T${jobname}.Rout"
eval $cmd

cmd="mv plot$2Agg_full_T${jobname}.Rout wave$2"
eval $cmd

