#################################
# Author: Bulent Ozel
# e-mail: bulent.ozel@gmail.com
# Collaborators:
#				Mario Eboli
#				Andrea Teglio
#				Andrea Toto
#################################


#### System Specific Configurations :::::

# Set working directory:
WORKING_DIR = "/Users/toto/Documents/R/Sim_FinNet_BW/"
setwd(WORKING_DIR)
rm(list = ls())

# Load the model library:
LIBRARY = "/Users/toto/Documents/R/Sim_FinNet_BW/scripts/ModelLibrary.R"
source(LIBRARY)

# Load the experiment library:
SIMULATOR = "/Users/toto/Documents/R/Sim_FinNet_BW/scripts/Simulator.R"
source(SIMULATOR)

# Set reporting options:
DATA_MODE = 1
# Only .csv files of the summary results are produced.

SCAN = TRUE
# FALSE means it reports T1 and Tfin
# TRUE means it produces T1,T2,.... Tfin.

# A whole simulation setup:

# Configure output locations:
TABLES_DIR = "/Users/toto/Documents/R/Sim_FinNet_BW/tables/"

# This configuration can be used to produce
# additional outputs (Balance Sheets) for inspection.
INSPECT_MODE = 0
PLOTS_DIR = TABLES_DIR


#### Running a Plan B centralization experiment on the ring 
# Note that shock type is "banks" where all external assets of bank are hit
# incrementallt until moving to the next bank in the permutation.

N = 64
Total.External.Assets = 100
amplitude = Total.External.Assets / (2 * N)


connectivity.and.centralization(
  N,
  KMin = 1,
  KMax = 1,
  rCentMin = 63,
  rCentMax = 63,
  eps = 0.1,
  phi = 0.4,
  A0 = Total.External.Assets,
  E0 = 10,
  b.sheet.equity = "hetero",
  d.min = 1,
  nRuns = 1000,
  shock.type = "banks",
  unit.shock = amplitude,
  r.shock = 1.0,
  n.sources = N
)

