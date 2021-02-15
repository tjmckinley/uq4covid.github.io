#!/bin/bash

## set up path to data
mkdir -p raw_outputs
export METAWARDSDATA=$HOME/Documents/covid/MetaWardsData

## set number of threads and processors
nprocessors=1
nthreads=4

## remove old plot files
rm *.pdf

## set up user inputs
cp ../user_inputs_default.txt user_inputs.txt
printf '\n.nseeds = 10\n' >> user_inputs.txt
printf '\n.ward_seed_filename = \"ward_seeds.csv\"\n' >> user_inputs.txt
printf '\n.age_seed_filename = \"age_seeds.csv\"\n' >> user_inputs.txt
printf '\n\n.contact_matrix_filename = \"contact_matrix.csv\"\n' >> user_inputs.txt

## set up contact matrix
rm contact_matrix.csv
cp ../../inputs/contact_matrix.csv .

## set up seeds
rm ward_seeds.csv
touch ward_seeds.csv
printf '1,1' >> ward_seeds.csv
rm age_seeds.csv
touch age_seeds.csv
printf '1,1\n' >> age_seeds.csv
printf '2,0\n' >> age_seeds.csv
printf '3,0\n' >> age_seeds.csv
printf '4,0\n' >> age_seeds.csv
printf '5,0\n' >> age_seeds.csv
printf '6,0\n' >> age_seeds.csv
printf '7,0\n' >> age_seeds.csv
printf '8,0\n' >> age_seeds.csv
   
## SEPID model
cp ../diseaseTest.dat disease.dat

# .pE .pEA 
printf '0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0 0 0 0 0 0 0 0 ' >> disease.dat
## .pP 
printf '0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ' >> disease.dat
## .pA 
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pI1 .pI1H .pI1I2 
printf '0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ' >> disease.dat
## .pI2
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pH .pHR 
printf '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ' >> disease.dat
## .lock_1_restrict .lock_2_release 
printf '1 1 ' >> disease.dat
## beta[1] beta[2] beta[3] beta[6]
printf '0.300630350 0.300630350 0 0 ' >> disease.dat
## repeats output
printf '10 test' >> disease.dat

## run MetaWards
metawards --nproc $nprocessors --nthreads $nthreads -m single -P 100000\
    -d ncov_age.json\
    -D demographics_age.json --mixer ../../model_code/mix_pathways\
    --input disease.dat\
    --iterator ../../model_code/iterator\
    -u user_inputs.txt -o raw_outputs --force-overwrite-output \
    --extractor ../../model_code/ward_extractor\
    --start-date 2020/01/01 --theme simple --nsteps 250

## run R script
#R CMD BATCH --no-restore --no-save --slave test.R 
