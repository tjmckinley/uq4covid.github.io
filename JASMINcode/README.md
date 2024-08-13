Here we can run various pieces of helper code to run an ensemble and post-process the model 
runs on JASMIN.

## Check design

To check a design that has been generated from a previous wave of history matching, and then 
to convert this to the correct format for the model, make sure the `.csv` file
is stored in a file e.g. `waveX/FILENAME.csv`, where `X` is the current wave
and `FILENAME` is the name of the file (e.g. `wave2Feb/inputsWave2.csv`). You need
to also give a `prevwave` argument to copy the `fixedInputs.txt` file across from.

Then call from the **main directory** you can call e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 2Feb 1Feb inputsWave2.csv' checkWavexDesign.R
```

This will look for the file `inputsWave2.csv` in the folder `wave2Feb`, and use this to create
the necessary files for running the ensemble. It will also check that the design is valid. 
Then it will copy the `fixedInputs.txt` file across from the `wave1Feb` folder. This way we 
can run multiple ensembles simultaneously using different naming conventions.

**Note:** the wave 1 design is already generated in the correct form so does not need this 
step.

## Setup code

To run a design or plotting code we can first run `setupSLURM.R`, passing in three to four
arguments. The first argument is the name of the folder containing the current wave design
(this will be appended to the word "wave", so if you want to run the design in `wave1Feb`, 
then you should pass the argument `1Feb`). The second argument should be one of `runDesign`, 
`runPlotSum` or `runPlotAgg` (see sections below for more details). The third argument is
the name of the outputs folder where the data are stored (this will be appended to the word
"outputs", so if the data are stored in `outputsFeb`, then you will have to pass the argument
`Feb`).

**Note**: the default max wall time is 35 mins (`00:35:00`). You can add a third argument to the
`setupSLURM` calls below to change this if required e.g. `01:00:00` will make into an
hour etc.

See specific examples below for more details.

## Design

To run a design use e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb runDesign Feb 18:00:00' setupSLURM.R
```

This sets up a file called `job_lookup_waveX.txt` containing IDs for running the code, and
a file called `submit_job_waveX.sbatch` (where `X` is replaced by the wave argument as above)
which can submitted to the SLURM scheduler. For example, the code above produces a file called
`submit_job_wave1Feb.sbatch`, which can be submitted to SLURM using:

```
sbatch submit_job_wave1Feb.sbatch
```

Jobs can be monitored using e.g.

```
squeue -u USERNAME
```

where `USERNAME` is replaced by your JASMIN username. Once jobs are run, results can be checked using `concatenateRuns.R`. This takes three arguments: the first is the wave name, the second is the output
folder name, and the third is a logical determining whether to set up failed runs ready to
be re-run. For example, if the design is held in `wave1Feb`, and the data in `outputsFeb`, then running

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb Feb FALSE' concatenateRuns.R
```

checks that runs have completed and tidies up the outputs if so. If successful,
then this returns an object `ll.rds` in the relevant wave folder with all the log-likelihood
estimates stored in it. If not, then it returns an error (see below).

The third argument is `FALSE` if you simply wish to return an error if the script fails, 
or if set to `TRUE` then this also recreates `job_lookup_waveX.txt` and `submit_job_waveX.sbatch` 
with the failed runs so that they can be easily resubmitted to the scheduler. 
(You might also want to set the time as an input e.g. `'--args 1 outputs TRUE 18:00:00'` else 
it defaults to `00:35:00`.) Then you can resubmit the missing jobs as before e.g.

```
sbatch submit_job_wave1Feb.sbatch
```

noting that this will just run the missing points and not all runs.

Alternatively, you can remove the failed points from the design completely (to avoid re-running).
To do this, you can use the `removeRun.R` script, which takes two arguments, the wave name
and the input to remove. The failed runs are printed in the `concatenateRuns.Rout` file for ease
of reference.

**Note**: if you wish to remove multiple points, then please remove one at a time and in 
**reverse order** (this is because when a point is removed from a dataset, all the indices
of subsequent points are shifted down by one, so doing in reverse order ensures that the
indices printed in `concatenateRuns.Rout` are correct at the point of being removed). For
example, to remove points 20 and 180 from the `wave1Feb` design, you would remove point
180 first, and then point 20, e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb 180' removeRun.R
R CMD BATCH --no-restore --no-save --slave '--args 1Feb 20' removeRun.R
```

This amends all of the saved file names appropriately, and also updates the `wave1Feb/inputs.rds` 
and `wave1Feb/disease.rds` files to remove the failed points.

Once you have re-run, or removed failed runs, then you can run the `concatenateRuns.R` script
again to produce the aggregated `ll.rds` file. This also produces a file called `times.rds`
that contains the runtimes for each design point.

**Note:** at this point the `inputs.rds` and `ll.rds` files can be used to perform history matching.
The rest of the code below is for producing summary trajectory plots across the ensemble.

## Forecasts (optional, and can be run after history matching has been completed)

To run a forecast use e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb runForecasts Feb 18:00:00' setupSLURM.R
```

As before, this sets up a file called `job_lookup_waveX.txt` containing IDs for running the code, and
a file called `submit_job_waveX.sbatch` (where `X` is replaced by the wave argument as above)
which can submitted to the SLURM scheduler. For example, the code above produces a file called
`submit_job_wave1Feb.sbatch`, which can be submitted to SLURM using:

```
sbatch submit_job_wave1Feb.sbatch
```

After the runs have completed, then running e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb Feb FALSE' concatenateForecasts.R
```

checks that forecasts have completed and tidies up the outputs if so. Failed runs can
be set-up and re-run in an analogous way to the `runDesign` code above.

**Note**: if you wish to produce plots of the forecasts, and have previously generated
plots from the initial model runs, then you will have to repeat all the steps below.

## Plotting summaries

Once a design has been run, we can also use JASMIN to generate summaries for plotting
an ensemble using e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb runPlotSum Feb 00:20:00' setupSLURM.R
```

which generates `job_lookup_wave1Feb.txt` and `submit_job_wave1Feb.sbatch` which can be 
submitted to the scheduler in the usual way. Currently this code will generate files called e.g. 
`plotSum_1_XXXX.rds` containing aggregated counts of different class / age / time combinations
for a given ensemble member (where `XXXX` is e.g. `natDeaths`, `ageRegionDeaths` etc.).

**Note**: if you want a subset of LTLAs to be extracted, then include a file called e.g.
`lads_OUTPUTS.txt` in the corresponding `OUTPUTS` folder, where each row contains a single 
LTLA ID to extract. For example, we can include a file called `lads_outputsFeb.txt` in the
`outputsFeb/` folder in the examples above if required. If this file is not present then 
the subsequent code only returns counts aggregated to the national level, if it is present 
then it produces counts at the national level as well as for the subset of chosen LTLAs.

The file `checkEns.R` will check the design and update the scheduler files on failure
if required (see above description of `concatenateRuns.R`) e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb Feb FALSE' checkEns.R
```

Note that these take less time to run so it's usually fine to simply rerun failed runs.

## Plotting the design

Finally, we can generate particle summaries for all age / class / time combinations and
combine them together using e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb runPlotAgg Feb 00:20:00' setupSLURM.R
```

in the usual way. Once the corresponding `submit_job_WAVE.sbatch` file has completed, then
`concatenateEns.R` can be used to check runs and concatenate together as seen fit e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb Feb FALSE' concatenateEns.R
```

Failed runs can be identified and re-run as described above. The output is a set of 
files called `sumEns_XXXX.rds` containing the summary data for plotting.

Once completed, the plot can be generated by e.g.

```
R CMD BATCH --no-restore --no-save --slave '--args 1Feb Feb NA NA FALSE' plotEnsembleTrajectories.R
```

where the first argument is the wave, the second is the name of the output folder, the 
third is the start of the forecasts (if forecasting, else set as `NA`), the fourth is the
last time point to plot (if set as `NA` then defaults to the longest point present in the
simulations). The fifth argument is `TRUE` if the hidden states are to be plotted, or
`FALSE` otherwise (note that this is only valid for simulated data where the "true"
trajectories are known).

This produces a series of PDF files but also a series of R files: `plot.rds`, `plots_sp.rds` 
and `plots_lad.rds` containing the `ggplot` objects to create each of the PDF files, which
are useful if you wish to reformat any of these plots for publication without having
to re-run all the plot generation code.

**Note**: the final `plotEnsembleTrajectories.R` call must use the `job_lookup_WAVE.txt` file
from the `runPlotAgg`/`setupSLURM.R` call above. Hence you might need to re-run this
original call if e.g. any of the original runs had failed and needed to be re-run (in which
case `job_lookup_WAVE.txt` would be different).

## Full simulation study files

For the full simulation study all the files above have to be amended slightly due to different
formatting of the data and outputs, therefore a set of corresponding files exist:

```
setupSLURM_full.R
concatenateRuns_full.R
removeRun_full.R
concatenateEns_full.R
plotEnsembleTrajectories_full.R
```

The `plotEnsembleTrajectories_full.R` only has four arguments since it is only run on simulated
data and therefore automatically plots the simulated data against the ensemble summaries.



