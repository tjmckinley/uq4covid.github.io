This repository contains the code and data used in [McKinley *et al.* (2024)]():

"On real-time calibrated prediction for complex model-based decision support in pandemics: Part II"

Trevelyan J. McKinley, Daniel B. Williamson, Xiaoyu Xiong, James M. Salter,
Robert Challen, Leon Danon, Ben Youngman and Doug McNeall

Please note that the simulator model runs were performed on the [JASMIN](https://jasmin.ac.uk/)
HPC, and so the code used to produce the model runs will contain JASMIN-specific
components and so should not be expected to run on other systems without some
modifications. These aspects are highlighted and discussed in more detail in the
relevant sections below.

## File descriptions

A brief description of the files in the repository is given below, and expanded upon in later
sections.

* `data/`: folder containing raw data files and scripts for processing the data into
  usable form for the model code.
* `historymatching/`: folder containing prototype history matching code.
* `inputs/`: folder containing additional data and scripts necessary for running the 
  models.
* `JASMINcode/`: folder of code required to setup and run ensembles on the JASMIN HPC.
* `BPF.cpp`: Rcpp code to run the particle filter (called by functions in `BPF.R` below).
* `BPF.R`: R functions to run the particle filters. Arguments to the main `BPF` function
  are defined at the top of the file.
* `BPF_full.cpp`: as above but for full simulation study.
* `BPF_full.R`: as above but for full simulation study.
* `stochModel.R`: code to simulate data.
* `stochModel_full.R`: code to simulate data.
* `checkWavexDesign.R`: code to check the design points at each new wave and convert the 
   design into the correct format for use in the model.
* `tnormConsts.h`: C++ header file necessary for the fast truncated Gaussian sampling.
* `wave1Design.R`: code to generate a wave 1 design across the input space.
* `wavexForecasts.R`: code to produce forecasts at any given wave.
* `wavexRuns.R`: code to run the model at a given design point.
* `wavexRuns_full.R`: code to run the model at a given design point for the full 
  simulation study.
  
## Workflow

The basic workflow is described below, but please read all the corresponding sections 
for details on each of these steps.

1. Run the `dataProcess.R` and `extractData.R` scripts in the `data/` folder to clean
   up the raw data and extract it into a usable form for the model.
2. Generate an initial wave 1 design using `wave1Design.R`.
3. If conducting a simulation study, then use `stochModel.R` and `stochModel_full.R`
   to generate simulated data.
4. Run the ensemble. For all the examples in the paper we used the JASMIN HPC,
   and so the `JASMINcode/` folder provides all the necessary code for setting up
   and running the ensemble. This uses `wavexRuns.R` and other associated files
   to run the model for each design point. The raw model outputs are then checked and
   collated into the correct format for history matching using other files in the
   `JASMINcode` folder, along with some summary plots of the ensemble trajectories.
5. Perform history matching. Some example code for wave 1 and then subsequent waves is
   provided in the `historymatching/` folder
6. Check the new design and set up the files in the correct format for the next wave
   using `checkWavexDesign.R`.
7. Return to step 4 and repeat until sufficient waves have been run.

## Hospitalisation, death and movement data

The `data/` and `inputs/` folders contain the raw data files, and some scripts for tidying them
up and extracting relevant time periods over which to fit the models. We have put
the original links where the data were downloaded below, but note that some of these are not accessible anymore. Nevertheless, the resulting downloaded data are included.

The raw data consist of:

* Deaths within 28 days of a positive test by LTLA: `data/ltla_2022-07-06.csv`, downloaded from
[https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=cumDeaths28DaysByDeathDate&format=csv](https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=cumDeaths28DaysByDeathDate&format=csv) [link defunct].
* Deaths within 28 days of a positive test by age and region: `data/region_2022-07-06.csv`, downloaded from [https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv](https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv) [link defunct].
* Cumulative hospital admissions by age and NHS region: `data/nhsRegion_cumAdByAge_2022-07-13.csv`, downloaded from [https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=cumAdmissionsByAge&format=csv](https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=cumAdmissionsByAge&format=csv) [link defunct].
* Total hospital cases by NHS region: `data/nhsRegion_hospCases_2022-07-13.csv`, downloaded from [https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=hospitalCases&metric=cumAdmissions&format=csv](https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=hospitalCases&metric=cumAdmissions&format=csv) [link defunct].

There is also:

* a lookup table to map 2019 LTLAs in England to regions: `data/Local_Authority_District_to_Region_(April_2019)_Lookup_in_England.csv`, downloaded from [https://geoportal.statistics.gov.uk/datasets/ons::local-authority-district-to-region-april-2019-lookup-in-england/explore](https://geoportal.statistics.gov.uk/datasets/ons::local-authority-district-to-region-april-2019-lookup-in-england/explore) [link defunct];
* a shapefile matching the 2019 LTLAs: `data/Local_Authority_Districts_(December_2019)_Boundaries_UK_BUC.zip`, downloaded from [https://geoportal.statistics.gov.uk/datasets/local-authority-districts-december-2019-boundaries-uk-buc/explore](https://geoportal.statistics.gov.uk/datasets/local-authority-districts-december-2019-boundaries-uk-buc/explore).

The movement data are converted from the 2011 census data available in the [MetaWards](https://metawards.org/)
package (see also [Woods *et al.* (2022)](https://joss.theoj.org/papers/10.21105/joss.03914)). 
These data have been mapped from 2011 electoral wards to 2019 LTLAs as described in [McKinley *et al.* (2024)]() 
and are stored as follows:

* `inputs/EW19.dat`: worker movements, with columns (home LTLA, work LTLA, number of movements);
* `inputs/PlayMatrix19.dat`: play movements, with columns (home LTLA, work LTLA, probability of movement);
* `inputs/PlaySize.dat`: number of play movements, with columns (home LTLA, number of movements);
* `inputs/LAD19_Lookup.csv`: a lookup table mapping the LTLA number used in the simulation model,
  to the LTLA codes used in the raw data.
  
These data require some further processing to be usable in the model. As described in [McKinley *et al.* (2024)]()
the movement data are available for all LTLAs in England and Wales (339 areas), but the hospitalisation
and death data are only available in England (315 areas). Furthermore, there are a couple of mis-matched 
area names between the different data sources, which were merged manually (and documented in the code
below). As such, the `data/dataProcess.R` file extracts and cleans up the raw data, and creates lookup
tables that are used in the model and for post-processing of some of the outputs. This can be run in the
terminal by navigating to the `data` folder and running the script file as a batch call e.g.

```
cd data
R CMD BATCH --no-save --no-restore --slave dataProcess.R
cd ..
```

or by changing the working directory of R to `data` and then using e.g.

```
source("dataProcess.R")
```

This will produce the following outputs:

* `age_lookup.rds`: a lookup table mapping the 8 age-classes used in the model to the 4 NHS age-classes. 
  This is stored as a `tibble` object with columns: `FID_death` (model age-classes) and `FID_nhsregion` 
  (NHS age-classes);
* `death_lad.rds`: a `tibble` object with columns: `date`, `deaths_1`, `deaths_2`, ..., `deaths_315`, 
  storing the cumulative number of deaths in each of the 315 English LTLAs over time;
* `death_lookup.rds`: a lookup table mapping the raw LTLA area codes and names to the 315 LTLA ID numbers 
  used in the model. This is stored as a `tibble` object with columns: `areaCode`, `areaName` and `FID`.
* `death_region.rds`: a `tibble` object with columns: `date`, deaths_1_1`, `deaths_1_2`, ..., 
  `deaths_8_9`, storing the cumulative number of deaths in each of the 8 age-classes and 9 regions
  over time;
* `lookup.rds`: a lookup table mapping the 339 raw LTLA area codes and names to the 315 LTLA ID numbers,
  9 regions and 7 NHS regions used in the model, with `NA`s denoting any areas that cannot be matched. 
  This is stored as a `tibble` object with columns: `FID`, `LAD19CD`, `LAD19NM`, `FID_death` and
  `FID_region`;
* `nhsage_lookup.rds`: a `tibble` object with columns: `age` and `FID` mapping the NHS age-classes to
  their corresponding ID number;
* `nhsregion_cumadage.rds`: a `tibble` object with columns: `date`, `hospInc_1_1`, `hospInc_1_2`, ...,
  `hospInc_4_7`, storing the cumulative hospital incidence in each of the 4 NHS age-classes and
  7 NHS regions over time;
* `nhsregion_hosp.rds`: a `tibble` object with columns: `date`, `hosp_1`, `hosp_2`, ..., `hosp_7`, 
  storing the hospital cases in each of the 7 NHS regions over time;
* `nhsregion_lookup.rds`: a `tibble` object with columns: `areaCode`, `areaName` and `FID` mapping the 
  raw region codes and names of the 7 NHS regions to their corresponding ID number;
* `region_lookup.rds`: a `tibble` object with columns: `FID`, `RGN19CD` and `RGN19NM` mapping the 
  raw region codes and names of the 9 regions to their corresponding ID number.
  
Once these have been produced, the `data/extractData.R` file then extracts the relevant data between two
dates, and creates an output folder that contains all the necessary files to run the models. This script
takes three input arguments, a start date and end date (in the form `%d/%m/%Y`) and an identifying name
(which is appended onto the word "outputs" before being used as the folder name). **Note that the output
folder is stored in the main directory and not `data`**. For example, to extract all data between the 15th
February and 23rd March 2020, and store this in a folder called `outputsFeb`, we can run the following
in a terminal window:

```
cd data
R CMD BATCH --no-save --no-restore --slave '--args 15/02/2020 23/03/2020 Feb' extractData.R
cd ..
```

Running the script without any arguments defaults to the ones above, and you can also run the script
directly in R but you will have to change the inputs on L20--22 manually. The `outputsFeb` folder then
contains the following files:

* `age_lookup.rds` (as above);
* `cumDeath_age_region.rds` (as `death_region.rds` above but restricted to required time period);
* `cumDeath_lad.rds` (as `death_lad.rds` above but restricted to required time period);
* `cumHospAd_age_nhsregion.rds` (as `nhsregion_cumadage.rds` above but restricted to required time period);
* `death_lookup.rds` (as above);
* `hosp_nhsregion.rds` (as `nhsregion_hosp.rds` above but restricted to required time period);
* `lads_outputsFeb.txt`: a text file containing the IDs of the ten LTLAs with the largest number of deaths
  at the end of the required time period (historically used for some diagnostic plots but left in for 
  legacy reasons);
* `Local_Authority_Districts_(December_2019)_Boundaries_UK_BUC.zip` (as above);
* `lookup.rds` (as above).

**Before running any model be sure to unzip the shapefile in the outputs folder.**

## Other data files, scripts and objects necessary for sampling from the input spaces

The `inputs/` folder also contains additional files required to run the model. These are:

* `inputs/age_seeds.csv`: this contains the proportion of the UK population in each of the age-classes of the model, derived from the Office for National Statistics. (2020). 2011 Census: Aggregate Data. [data collection]. UK Data Service. SN: 7427, DOI: [http://doi.org/10.5257/census/aggregate-2011-2](http://doi.org/10.5257/census/aggregate-2011-2).
* `inputs/coMix_matrix.csv`: this is a contact matrix between different age-classes used *after* the first lockdown.
  The original data (`20200327_comix_social_contacts.xlsx`) are from [Jarvis *et al.* (2020)](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01597-8#availability-of-data-and-materials) and
  were downloaded from here 
  
  [https://cmmid.github.io/topics/covid19/comix-impact-of-physical-distance-measures-on-transmission-in-the-UK.html](https://cmmid.github.io/topics/covid19/comix-impact-of-physical-distance-measures-on-transmission-in-the-UK.html)
  
  and we used the `All_contacts_imputed` sheet with column `A` and row `1` removed.
* The `inputs/dataTools.R` file contains helper function required to run the model.
* `inputs/hospStays.rds`: this is a fitted finite mixture model (FMM) used to sample from the input
  space for the parameters guiding the length of hospital stays.
* `inputs/hospThresh.rds`: is a threshold on the p.d.f. of the FMM in `inputs/hospStays.rds`, which
  is used to define a boundary for the input space.
* `inputs/parRanges.rds`: contains the input ranges for the parameters with a hypercube input space.
* `inputs/pathways.rds`: this is a fitted FMM used to sample from the input space for the parameters 
  guiding the probabilities of transitioning down different epidemiological pathways.
* `inputs/pathThresh.rds`: is a threshold on the p.d.f. of the FMM in `inputs/pathways.rds`, which
  is used to define a boundary for the input space.
* `inputs/POLYMOD_matrix.csv`: this is a contact matrix between different age-classes used *before* 
  the first lockdown. This was extracted from the [`socialmixr`](https://cran.r-project.org/web/packages/socialmixr/vignettes/socialmixr.html) package, which uses data in [Mossong *et al.* (2008)](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050074).

## Wave 1 design

The `wave1Design.R` file contains code to generate a Wave 1 space-filling design
across the original input space. This uses a Latin Hypercube design for some 
parameters, and then a space-filling design for the parameters with non-rectangular input
spaces, as described in [Williamson *et al.* (2024)](). This script takes a single argument
which corresponds to a name for the outputs (which is appended to the word "wave" and the wave
number to store the model runs in a standardised way across the waves). For example, to produce
a wave 1 design in a folder called `wave1Feb`, you can run the following from a terminal window:

```
R CMD BATCH --no-save --no-restore --slave '--args Feb' wave1Design.R
```

Alternatively, the script can be run from within R in the usual way, noting that you will have to
change the default wave name manually. 

**Important**: this script also sets up file called `fixedInputs.txt`, which is stored within
e.g. the `wave1Feb` folder. This contains additional data-specific inputs necessary to run the models.
These are defined on lines 28--46 of `wave1Design.R`, and correspond to:

* `tstart`: the starting time for the simulations (where day 0 corresponds to the `start` time
  used in the `extractData.R` call above).
* `tstop`: the end time (in days) over which to fit the model (must be $\leq$ the `end` time
  used in the `extractData.R` call above).
* `lockdown_day`: the day at which lockdown is introduced (relative to `tstart`).
* `npart`: the number of particles.
* `niter`: the number of iterations used in the MCMC updates to alleviate particle degeneracy.
* `a1`: parameter defining observation error at the LTLA-level.
* `a2`: parameter defining observation error.
* `b1`: parameter defining observation error.
* `b2`: parameter defining observation error.
* `a_dis`: parameter defining the amount of model discrepancy in age-classes 1--7.
* `b_dis`: parameter defining the amount of model discrepancy in age-classes 1--7.
* `b_dis_8`: parameter defining the amount of model discrepancy in age-class 8.
* `sigma2_lad`: the variance of the aggregation error placed on the observations at the LTLA-level.
* `sigma2_age_region`: the variance of the aggregation error placed on the observations at the age/region-level.
* `sigma2_nhsregion`: the variance of the aggregation error placed on the observations at the NHS region-level.
* `sigma2_age_nhsregion`: the variance of the aggregation error placed on the observations at the NHS age/region-level.
* `saveAll`: a logical determining whether to save the particle trajectories.
* `snapshot`: a logical determining whether to save the state of the system at the end of the run (used to 
  initiate any forwards simulations).
* `writeExt`: a logical determining whether to write the particle trajectories to external files, rather 
  than keeping them in R's internal memory (necessary for large models).
  
These lines will have to be adjusted manually as required before running the script. You will also need 
to check the prior ranges on L57--62, the number of design and validation points on L68 and L70
respectively.

The contents of the `wave1Feb` folder are:

* `design.pdf`: a plot of the wave1 design.
* `disease.rds`: a `tibble` object containing the design points but converted to the necessary
  form for use in the model.
* `fixedInputs.txt`: the fixed inputs for use in the model as described above.
* `inputs.rds`: a `tibble` object containing the design points.

## Note on memory

For the rest of the commands below, if you are on Unix, then you *may* have to increase the 
stack limit before loading R e.g. something like:

```
ulimit -s unlimited
```

**We take no responsibility if you choose to do this, but it allowed us to run the code 
without any issues.**

To run in parallel using OpenMP, you will also have to have the correct compiler tools installed. 
However this hasn't been tested in anger, since for the full model we use the JASMIN HPC
and run each model on a single core and parallelise across nodes instead of using OpenMP.

## Simulation study

The code used to simulate a synthetic data set from the underlying model can be found
in the `stochModel.R` file. This also requires a set of input parameters to use for
the simulation. We did this by selecting a set from an initial set of design points generated
using the `wave1Design.R` file. 

The process is to run the `wave1Design.R` file to generate a Wave 1 LHS design, then run 
run `stochModel.R` to simulate data using one of the wave 1 design points. It is worth noting
that `stochModel.R` copies across the information contained in the corresponding `fixedInputs.txt`
file, and so if you want to simulate over longer time periods than you wish to fit to (perhaps
if you want to assess forecasts against the actual simulated data), then you will have to 
amend the `wave1Design.R` code as described in the previous section in order to change e.g. the 
`tstop` parameter.

The `stochModel.R` script takes three arguments, the first is the seed, and the second is the 
ID of the folder containing the parameter sets from which to choose from, and the third is the
specific parameter set to choose. For example, to generate a set of candidate parameter sets 
for a simulation over 54 days say, then we can amend the `wave1Design.R` file and set
`tstop <- 54` as described above, and then run the following in a terminal:

```
R CMD BATCH --no-save --no-restore --slave '--args Sim' wave1Design.R
```

We can then run

```
R CMD BATCH --no-save --no-restore --slave '--args 456 Sim 10' stochModel.R
```

which sets the seed to be `456`, chooses the set of design points in the `wave1Sim` folder
as candidate parameter sets, and then chooses the tenth of these (`10`) as the "true" 
parameters. This script will then create a new output folder called `outputs456`, containing 
the simulated data in the same format as the real data. This runs a small number of replicates 
(8 in total), and then picks the run that is the closest in L2 distance to the median across 
the 8 replicates as a representative model run. Since we have a small number of replicates we 
can keep all the trajectories in memory and so do not have to write them to external files.

The contents of the output folder (e.g. `outputs456`) is the same as for the real data, except for
additional files:

* `disSims.rds`: a `tibble` object storing the simulated hidden states of the system, with 
  the first column being the time point `t`, and the subsequent columns having the form 
  `CLASS_AGE_LTLA`, for example `DH_8_339` would be the number of people who have died in hospital
  from the 8th age-class in the 339th LTLA.
* `pars.rds`: a `tibble` object containing the parameters used to simulate the data.
* `simsNational.pdf`: a plot of the simulated data against ensemble summaries from the 8 replicates.
* `simsTopLADs.pdf`: a plot of the hidden states for the ten LTLAs with the highest number of deaths
  at the time the simulation stopped.

**Important**: if you want to fit to a shorter time period than that used for the simulation,
you can amend the `tstop` parameter in the wave 1 `fixedInputs.txt` file, and/or generate
a new set of wave 1 design points for use in the subsequent history matching.

### Full simulation study

For the simulation study where we fit to the same level of aggregation as the model, we do not
re-simulate the data, we simply take the previous simulated hidden states and create a set
of new data objects for use in the simulation study. To do this run the `stochModel_full.R`
file, which takes three arguments, the seed, the name of the original simulated output folder,
and the name of the new folder containing the expanded simulations. For example,

```
R CMD BATCH --no-save --no-restore --slave '--args 456 Sim FullSim' stochModel_full.R
```

will take the simulations stored in `outputs456` and the design points stored in `wave1Sim`,
and create a new output folder called `outputsFull456` and a new design folder called
`wave1FullSim` containing the necessary files. Note that since we now have observed data at
a different resolution (the age/LTLA level), we only have one source of aggregation error
now, operating at the age/region level, and so simplify the `fixedInputs.txt` file accordingly.

The `outputsFullSim` folder contains the following different files:

* `cumDH_age_lad.rds`: contains cumulative hospital deaths over time at the age / region level.
* `cumDI_age_lad.rds`: contains cumulative community deaths over time at the age / region level.
* `cumH_age_lad.rds`: contains cumulative hospital incidence over time at the age / region level.
* `H_age_lad.rds`: contains number of people in hospital over time at the age / region level.
* `simsComparison.pdf`: a plot comparing the aggregated counts at the levels observed in the
  real data for the original simulated data versus the new simulated data accounting for the
  slightly different aggregation errors. This is just a check that the counts are close to each 
  other, and that the change in aggregations hasn't changed the underlying dynamics significantly.
  
## Running the model

The `wavexRuns.R` file provides code for running a single design point. This takes three arguments,
the "wave" folder where the design points are stored, the "hash" number relating to which design point
to extract, and the "outputs" folder denoting where the data are stored. For example, to run the first
design point of the wave 1 design (stored in the `wave1Feb` folder), for the data stored in the
`outputsFeb` folder, we can run:

```
R CMD BATCH --no-save --no-restore --slave '--args 1Feb 1 Feb' wavexRuns.R
```

Note that this automatically saves the outputs to external files due to memory considerations. Each
set of trajectory files is saved in a folder called `saveOut_WAVE_X`, where `WAVE` is the wave name 
and `X` is the hash number. These folders are saved inside the corresponding `WAVE` folder, e.g. 
`wave1Feb/saveOut_wave1Feb_1` etc. The log-likelihood estimates are stored in files called 
`runs_md_X.rds` e.g. `wave1Feb/runs_md_1.rds` etc.

## Running and visualising an ensemble on JASMIN

The code in the `JASMINcode/` folder provides wrapper functions to setup the necessary files to run
the entire ensemble on JASMIN. This workflow is documented in a separate [README.md](JASMINcode/README.md).
This also provides code to return the estimated log-likelihoods to be used in the subsequent history
matching steps, as well as to produce summary plots of the ensemble trajectories. Note that occasionally
runs fail, and so we also provide code to check for runs that have not completed and re-run
or remove them accordingly (since the ensembles take a long time to run, it is often preferable to
remove failed runs rather than rerun them). This checking code then produces an updated input file that
matches the log-likelihood estimates with the failed points removed.

## History matching

The `historymatching/` folder contains protoype code for building and validating emulators, and for
generating new designs. 

**Note**: these files will have to be copied and amended for different waves, and the wave 1 code works
slightly differently to the subsequent waves. 

### Wave 1

For the Wave 1 emulators, copy the `historymatching/emulateDGP_wave1.R` and `historymatching/emulator_fns.R`
files into the folder containing the `inputs.rds` and `ll.rds` files containing the design and 
log-likelihood estimates for the ensemble respectively. Then navigate to the ensemble folder for the 
next steps e.g. in a terminal

```
cp historymatching/emulateDGP_wave1.R wave1Feb
cp historymatching/emulator_fns.R wave1Feb
cd wave1Feb
```

(Note that you only need these two files to run the history matching, so if you are running ensembles
on an external cluster, you only need to copy those files across. Please note also that the code matches
entries in `inputs.rds` to entries in `ll.rds` by order (so both files must have the same number of entries
and be in the same order). (The run code in the `JASMINcode/` folder does this automatically.)

Once these files have been copied across, you will need to amend some of the lines of the `emulateDGP_wave1.R` to ensure they match to your data/inputs. Firstly, check the prior ranges on L16--44,
then check the number of training and validation points on L55 and L61 respectively (these may have to
be amended if you removed any failed runs when running the ensemble). The code requires some manual
checking as it is running, and so it is best to run this code interactively in R.

To choose the active variables for the deep Gaussian Process (DGP) emulator, we first fit a standard
GP, and then use a threshold such that inputs with a length-scale parameter of less than the threshold
are considered "active". This threshold is set on L88 and should be set by the user and not just used 
blindly. We used values of 500 for most cases in our model runs. These choices can also be guided by
the validation plots.

Once we apply the threshold we then refit the GP and check that all length-scales are less than the 
threshold after refitting. If they are not, then we can apply the threshold again and carry on this 
process. If you need to do this then the code on L95--98 can be uncommented and repeated as required.

Once this has been done, and the active variables selected, then we can fit the DGP to the active 
variables using MCMC, which is done on L104--112. As with any MCMC method it is useful to assess
convergence, and so if the trace plots do not seem to have converged then the chain can be run for 
longer. If you need to do this then the code on L122--126 can be uncommented and repeated as required.

**Note**: for many of the runs we conducted, the mixing of the chains aren't as good as we would 
ideally like. However, the MCMC here is used mainly to find the area of high posterior mass, and the
predictions from the DGP actually only a point prediction corresponding to the median value from the 
last half of the chain. To improve mixing the chain could be run for longer and then thinned, but 
this would take a long time for large-scale problems, and so instead here we take a more pragmatic 
approach and run it until it looks like it has converged at least in the last half of the chain. 
We then produce validation plots, and if the validation looks OK we proceed with the rest of the 
process. We found in tests that running the chains long enough to get good mixing did not noticably 
improve the validation, as long as the median value was roughly correct.

The code then stores any "doubt points" on L175--208 (as described in the paper), and then builds 
a custom emulator object on L211--238. The `hmer` package allows us to build custom emulators, using
a function called `Proto_emulator()`, but since the whole history matching approach is different to the 
standard approach, we wrap all of the custom implausibility measures, target log-likelihoods, and
necessary files for sampling from the non-hypercube input spaces into this object. The `emulator_fns.R`
file contains custom functions to do all of these steps, which are used in the construction of the
`emulator` object.

Then we run a final check to ensure that the current best point is retained in the NROY space 
(L241--248). A new design is then generated by initially sampling a large number of points from the 
original input space (1,000 in this case; L253--323), and then extracing the subset of these that are
retained in the NROY space. The space removed at this wave is estimated on L326--329 (noting that the "target" log-likelihood is stored in the `emulator` object, and as such we have to set a dummy `targets`
object with `val = 1` in order for the code to work. This is because, as described above, we have many
context-specific things used here which are different to the usual HM approach (such as non-hypercube
input spaces, use of Voronoi regions to retain space around "doubt points", custom implausibility measures
etc.), and as such we have to employ a few tricks to allow `hmer` to handle all of these aspects correctly.
L333--334 then extract the subset of the baseline points that remain in the NROY space, which we use as
a basis for building a new set of 1,000 points in the NROY space by utilising the slice sampling method
implemented in `hmer` (L337--338). The Wave 2 training and validation points are then sub-sampled from 
this set of baseline points using a maximin design (L341--350). If you want to change the number of 
training and validation points, then you can do this on L341 and L344 respectively (noting that if you 
want more than 1,000 design points, then you will have to generate a new set of baseline points
accordingly).

Finally, the new design is plotted (L353--371), and then to aid subsequent waves we then augment the
set of baseline NROY points back up to 1,000 to replace the points removed for the new design (L374--373).
The `inputsWave2.csv` then contains the training and validation points for Wave 2, and so this can be
copied into a new folder with the same naming convention as used in Wave 1 (so if the Wave 1 runs are
stored in `wave1Feb`, then create a new folder called `wave2Feb` and copy the `inputsWave2.csv` file
into this folder before running the new wave).

### Wave 2 onwards

For the Wave 2+ emulators, copy the `historymatching/emulateDGP_wavex.R` and 
`historymatching/emulator_fns.R` files into the folder containing the `inputs.rds` and `ll.rds` files
containing the design and log-likelihood estimates for the ensemble respectively. Then navigate to the
ensemble folder for the next steps e.g. in a terminal

```
cp historymatching/emulateDGP_wavex.R wave2Feb
cp historymatching/emulator_fns.R wave2Feb
cd wave2Feb
```

This contains the same basic structure as the Wave 1 code described above, but with a few key differences.
As in Wave 1, the code requires some manual checking as it is running, and so it is best to run this code
interactively in R.

The first thing to note is that since a consistent naming convention is used for the folders containing 
the outputs of each wave, then the `emulateDGP_wavex.R` code can automatically link to the outputs of the
previous waves where required (and can copy across seeds, prior ranges etc.). As such, you will need to
amend L16--17 to set the naming convention (`wave_name`) and a vector of wave numbers evaluted so far
(**including the current wave**); these are defined in the `wave_nos` object. For example, if you are at
Wave 3, and the waves are stored in `wave1Feb`, `wave2Feb` and `wave3Feb`, then you would set 
`wave_name <- "Feb"` and `wave_nos <- 1:3`.

Then the code re-builds the emulators for the previous waves, which are required to sample new points and 
check the NROY space etc. Note that in standard `hmer` code, we could save the previous wave emulators and
simply load them back in again as a single object. The complexity here is that the emulators are built using
the `dgpsi` package, which uses Python, rather than R, under-the-hood. As such we can save the Python
objects as `.pkl` files, but we cannot include the Python objects embedded in an R object saved as an `.rds`
file. We overcome this by saving all of the relevant components at each wave, and then looping over the
previous waves to rebuild the emulators. This is done on L24--61.

Then, to ensure maximising the information in the training data, we load in all the previous wave design
points, and then extract any of these that remain in the previous wave NROY space. We can then augment the
training data with these additional points. This whole process is done on L64--127. As before, check the
number of training and validation points on L95 and L100 respectively (these may have to be amended if you
removed any failed runs when running the ensemble).

Then the code runs similarly to Wave 1, where we define active variables (L130--143), build the DGP emulator
using the active variables (L152--207), check and record "doubt points" (L219--252), set up the new
custom emulator object for use in `hmer` (L258--293) and check that the best point is retained in the
previous NROY space (L296--309).

Then we load in the set of baseline NROY points from the previous wave, and use these to estimate the
space removed at the current wave (L312--318). Then, as before we extract the subset of the previous
wave baseline points that remain in the updated NROY space, and then use these to generate 1,000 baseline
points in the updated NROY space (L322--327). These are then sub-sampled using a maximin design to generate
a new set of training and validation points (L330--339). Again, the number of training and validation
points can be amended on L330 and L333 respectively if required. The new design points are then plotted
(L342--360), before the updated set of baseline NROY points are augmented back up to 1,000 to replace the
points removed for the new design (L363--364).

The resulting `inputsWaveX.csv` file then contains the training and validation points for the new Wave `X`,
and so this can be copied into a new folder with the same naming convention as used in the previous waves
(so if we are at the Wave 3, and the runs are stored in `wave1Feb`, `wave2Feb` and `wave3Feb`, then create 
a new folder called `wave4Feb` and copy the `inputsWave4.csv` file into this folder before running the new
wave).

## Forecasts

The `wavexForecasts.R` file provides code for running forecasts for a single design point. For a given
design point, the states of the system at a given time point can be obtained from the saved particle
filters, and then from these we can simulate forwards in time from the underlying model in order to produce 
forecasts. To do this you will have to create a new file called `fixedInputs_cont.txt`, inside of each
`waveX` folder. This should be identical to the corresponding `fixedInputs.txt` file, except with different
`tstart` and `tstop` arguments, where the former denotes the final set of data points that the model
is currently fitted to, and the latter to when you would like the forecasting to stop. For example, to run
a set of 14-day forecasts for the wave 1 design (stored in the `wave1Feb` folder), firstly copy the
`wave1Feb/fixedInputs.txt` file to `wave1Feb/fixedInputs_cont.txt`, and then amend the first two rows of
`wave1Feb/fixedInputs_cont.txt` from 

```
0
37
```

to

```
37
51
```

Then, the `wavexForecasts.R` script takes three arguments, the "wave" folder where the design points are 
stored, the "hash" number relating to which design point to extract, and the "outputs" folder denoting where
the data are stored. For example, to run the first design point of the wave 1 design (stored in the `wave1Feb`
folder), with the data stored in the `outputsFeb` folder, we can run:

```
R CMD BATCH --no-save --no-restore --slave '--args 1Feb 1 Feb' wavexForecasts.R
```

As before, this automatically saves the outputs to external files due to memory considerations, and in fact
just adds the forecasts to the existing directories.

Similarly, the code in the `JASMINcode/` folder provides wrapper functions to setup the necessary files to run
the entire set of ensemble forecasts on JASMIN, as well as to visualise the forecasts. For full details, please
see the separate `JASMINcode/README.md` file.



