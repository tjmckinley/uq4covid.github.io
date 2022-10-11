Here we can run various pieces of helper code on JASMIN.

To run a design or plotting code we can first run `setupSLURM.R`, passing in `wave` and
`runCode` arguments: `wave` should be numerical, matching where the ensemble design is 
held (the design must be in a folder called e.g. `wave1`). The `runCode` argument should 
be one of `runDesign`, `runPlotSum` or `runPlotAgg`.

**NOTE**: the default max wall time is 35 mins. You can add a third argument to the
`setupSLURM` calls below to change this if required e.g. `01:00:00` will make into an
hour etc.

## Check design

To check a design that has been passed from a previous wave of HM, and then to
convert this to the correct format for the model, make sure the `.csv` file
is stored in a file e.g. `wavex/FILENAME.csv`, where `x` is the current wave
and `FILENAME` is the name of the file (e.g. `wave2/inputsWave2.csv`). You need
to also give a `prevwave` argument to copy the `fixedInputs.txt` file across from.

Then call e.g. from the **main directory**:

```
R CMD BATCH --no-restore --no-save --slave '--args 2 1 inputsWave2.csv' checkWavexDesign.R
```

## Design

To run a design use e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1 runDesign outputs 18:00:00' setupSLURM.R
```

The first `args` argument is the wave number `x` (converted to `wavex` to specify the folder
to find the inputs and save the outputs into. The next argument specifies that we wish to
run a design (see other options below). Then we have the name of the file to find the data in
(in this case `outputs`), and finally we have an optional argument setting the maximum wall
time (defaults to `00:35:00` if left unspecified).

This sets up a file called `job_lookup_wavex.txt` containing IDs for running the code, and
a file called `submit_job_wavex.sbatch` (where `x` is replaced by the wave argument as above)
which can submitted to the SLURM scheduler e.g.

```
sbatch submit_job_wave1.sbatch
```

Jobs can be monitored using e.g.

```
squeue -u USERNAME
```

Once jobs are run, results can be checked using `concatenateRuns.R` e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1 outputs FALSE' concatenateRuns.R
```

This checks that runs have completed and tidies up the outputs if so. If successful,
then this returns an object `ll.rds` in the relevant folder with all the log-likelihood
estimates stored in it. If not, then it returns an error (see below).

The first argument is the `wave`, the second is `outputs`, and the third is `FALSE` 
if you simply wish to return an error if the script fails, or if set to `TRUE` then this 
also recreates `job_lookup_wavex.txt` and `submit_job_wavex.sbatch` with the failed runs 
so that they can be easily resubmitted to the scheduler. (You might also want to set the 
time as an input e.g. `'--args 1 outputs TRUE 18:00:00'` else it defaults to `00:35:00`.)

## Forecasts

To run a forecast use e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1 runForecasts' setupSLURM.R
```

This sets up a file called `job_lookup.txt` containing IDs for running the code, and
a file called `submit_job.sbatch` which can submitted to the SLURM scheduler.

```
sbatch submit_job.sbatch
```

Jobs can be monitored using e.g.

```
squeue -u USERNAME
```

Once jobs are run, results can be checked using `concatenateForecasts.R` e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1 FALSE' concatenateForecasts.R
```

This checks that runs have completed and tidies up the outputs if so. If not successful,
then this returns an error (see below).

The first argument is the `wave` and the second is `FALSE` if you simply wish to return 
an error if the script fails, or if set to `TRUE` then this also recreates `job_lookup.txt`
and `submit_job.sbatch` with the failed runs so that they can be easily resubmitted to
the scheduler. The optinal third argument is the wall-time to pass to `submit_job.sbatch`.

## Plotting summaries

Once a design has been run, we can also use JASMIN to generate summaries for plotting
an ensemble using e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1 runPlotSum 00:05:00' setupSLURM.R
```

which generates `job_lookup.txt` and `submit_job.sbatch` which can be submitted to the
scheduler in the usual way. Currently this code will generate files called e.g. 
`plotSum_1.rds` containing aggregated counts of different class / age / time combinations
for a given ensemble member.

**NOTE**: if you want a subset of LADs to be extracted, then include a file called
`lads.txt` where each row contains a single LAD ID to extract. This file needs to be
created manually. If it is not present then the subsequent code only returns
counts aggregated to the national level, if it is present then it produces counts
at the national level as well as for the subset of chosen LADs.

The file `checkEns.R` will check the design and update the scheduler files on failure
if required (see above description of `concatenateRuns.R`) e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1 FALSE' checkEns.R
```

## Plotting the design

Finally, we can generate particle summaries for all age / class / time combinations and
combine them together using e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1 runPlotAgg 00:03:00' setupSLURM.R
```

in the usual way. Once the corresponding `submit_job.sbatch` file has completed, then
`concatenateEns.R` can be used to check runs and concatenate together as seen fit e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1 FALSE' concatenateEns.R
```

Then the plot can be run by e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1 51 NA' plotEnsembleTrajectories.R
```

where the first argument is the wave, the second is the start of the forecasts (if forecasting,
else set as `NA`), and the third is the final time point to plot (if set as `NA` then defaults
to the longest point present in the simulations).

**Note**: the final `plotEnsembleTrajectories.R` call must use the `job_lookup.txt` file
from the `runPlotAgg`/`setupSLURM.R` call above. Hence you might need to re-run this
original call if e.g. any of the orignal runs had failed and needed to be re-run (in which
case `job_lookup.txt` would be different).


