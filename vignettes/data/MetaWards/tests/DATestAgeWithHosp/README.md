This folder contains a simulation study to test the PF
approach to calibrate models and perform DA simultaneously.

Firstly, run `wave1Design.R` file to generate a 
Wave 1 LHS design, along with the non-standard samplers. 

Then run `stochModel.R` to simulate data 
using one of the MW design points. Outputs are placed in the 
`outputs` folder (which is created if it is not present).

You can run `wavexRuns.R` to run the Wave 1 design.

Then you have to build your own emulators etc., and ways of sampling
from the NROY space. The `wavexDesign.R` file gives some examples
of how to test any NROY samples against the prior constraints in order
to generate the subsequent wave designs (these are just to show the functions)
and guide you as to what you need to amend. Change the `wave <- 2` on L6
to set the correct wave number.

Once a design has been generated, running the `wavexRuns.R` file and changing
`wave <- 1` on L17 will point to the new wave design.

## Tests

Some tests are contained in the `tests` folder. The working directory
must be set to `tests` before running.

The file `testPF.R` gives some commands for reading
in the data and the current design and testing the PF. The code 
is shown for running the model and PF and generating
log-likelihood estimates for the design, as well as producing
some plots of the particle densities (unweighted, so only an 
approximation here) against the data for comparison, if
you wish to do this. The PF code is in the `PF.R` file, the 
truncated Skellam sampler is in the `trSkellam.R` file.

The file `testPF1.R` does the same but for the new particle
filter sampler, where the code is in the file `PF1.R`.

The file `testdiffPFs.R` contains test code to compare the empirical
log-likelihood distribution from the original PF, the new PF, and
the original PF with larger numbers of particles in order to
assess how the estimates change.
