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
WORKING_DIR = "~/Documents/Sim_FinNet/"
setwd(WORKING_DIR)
rm(list = ls())

# Load the model library:
LIBRARY = "~/Documents/Sim_FinNet/scripts/ModelLibrary.R"
source(LIBRARY)

# Load the experiment library:
SIMULATOR = "~/Documents/Sim_FinNet/scripts/Simulator.R"
source(SIMULATOR)

# Set reporting options:
DATA_MODE = 1
# Only .csv files of the summary results are produced.

SCAN = TRUE
# FALSE means it reports T1 and Tfin
# TRUE means it produces T1,T2,.... Tfin.

# A whole simulation setup:


TABLES_DIR = "~/Documents/Sim_FinNet/tables/"
PLOTS_DIR = TABLES_DIR
# Configures the output locations.

INSPECT_MODE = 1
# This configuration can be used to produce
# additional outputs for an inspection:
# - Balance Sheet
# - Plot of the network structure
# - The final state of the contagion for the latest run for a given
# network and balance sheet initialization. This is the distribution of shock absorption in the network.
# Shock is eventually absorbed by share holders (equity, the r.Absorb.E column) and 
# households (deposits, r.Flow.Out column). Note that depending on the
# (1) network structure, (2) balance sheet initialization, (3) the shock application type and
# (4) the vector of the shock each run may have a different final state of the diffusion.
# Please note that the "Shock" column of the output denotes total internal and external shock received
# by the bank.


##### Explanations for Inputs :::::::::

## A::: Required inputs
# 1. Number of banks: N

# Connectivity::
# 2. Number of neighbors in the wheel: KMin
# 3. Number of neighbors in the wheel: KMax
# Note that KMin and KMax together can be used to trace a range of connectivity.

# Centralization::
# 4. The minumum centralization multiplicant: rCentMin
# 5. The minumum centralization multiplicant: rCentMax
# Note that these inputs together can be used to trace a range of centralization.

# Network Symmtery:
# 6. Whether the network is forced to be bi-directional: symmetric
# Note that chen rCentMax >= rCentMin > 0, symmetric is set to be 1.

# Balancesheet Configurations::
# 7. Ratio of equity to total asset: eps
# 8. Ratio of external debt to internal debt: phi
# 9. Exogenously set total external assets: A0

# 10. Exogenously set total capital: E0
# This is used only when b.sheet.equity is one of "homo" or "planA.spec". 
# It is a redundant input when b.sheet.equity is "hetero".

# 11. b.sheet.equity can be one of:
# - "homo", denoting that each bank has the same and equal equity, which is E0/N.
# - "hetero", denoting that equity of the bank depends on phi, eps and network position.
# - "planA.spec", holds only when K = 1, for K > 1 and it switches to "homo" configuration.

# When b.sheet.equity is "hetero" the set of equations for balance sheet initializations are as follows:
# A + C = H + D + E
# C = D
# D / H = phi
# E / ( A + C) = eps
# A = A0

# When b.sheet.equity is "homo" the set of equations for balance sheet initializations are as follows:
# A + C = H + D + E
# C = D
# A = A0
# E = E0
# D = d.min * n.links (number of links in the network)
# a_i = A / N
# e_i = E / N
# h_i = H / N
# c_i = d_i = out.degree_i * d.min, (in.degree_i = out.degree_i) 

# When b.sheet.equity is "planA.spec" the set of equations for balance sheet initializations are as follows:
# A + C = H + D + E
# C = D
# A = A0
# E = E0
# a_i = A / N
# e_i = E / N
# h_i = H / N
# d_ring_i = d_pendant_j = 2 * d.min 
# d_star = 2 * |pendant| * d.min + d.min * (N - |pendant|)
# This configuration is applied only when K = 1.

# 12. Minimum link weight: d.min
# This is used only when b.sheet.equity is one of "homo" or "planA.spec". 
# It is a redundant input when b.sheet.equity is "hetero".
# Note that when it is "hetero", d.min is calculated endogenously and 
# It is equivalent to D / n.links. 

# 13. Number of repetion for the same network config: nRuns
# This corresponds to number of shock permutation is used.

# 14. The shock application method: shock.type
# The shock.type is one of:
# - "assets", the set of external assets is permuted for a shock sequence.
#             Each bank owns a distinct subset of assets.
# - "banks", the set of banks in the system is permuted for a shock sequence.
#             Set of all assets of a bank which is picked is incrementally hit before moving to the
#             next bank in the sequence.
# - "full", the set of banks in the system is permuted for a shock sequence.
#             All assets of a bank which is picked all at once.


# 15. The amplitude of unit of a shock: unit.shock
# Note that this is used when shock.type is one of "assets" or "banks".
# When It is "full" this endogenously set to be r.shock of  a_i of the shock recepient bank.

# 16. The rate of external shocks received by a bank. r.shock
# Note that r.shock default value is 1.0 and is active only if shock.type is "full".
# If exogenoysly picked a partial shock to be applied this input can be used.

# 17. The number of banks to be hit: n.sources
# Note that the default value is N the number of banks in the system.
# If only a random subset of nodes to be hit is desired this input can be used.




# Example Configurations :::::

# Default: 
N = 64
Total.External.Assets = 100
amplitude = Total.External.Assets / (2 * N)


connectivity.and.centralization(
  N,
  KMin = 1,
  KMax = N,
  rCentMin = 0,
  rCentMax = N,
  symmetric = 1,
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

#### Running a connectivity experiment

N = 64
Total.External.Assets = 100
amplitude = Total.External.Assets / (2 * N)


connectivity.and.centralization(
  N,
  KMin = 1,
  KMax = N - 1,
  rCentMin = 0,
  rCentMax = 0,
  symmetric = 0,
  eps = 0.1,
  phi = 0.4,
  A0 = Total.External.Assets,
  E0 = 10,
  b.sheet.equity = "hetero",
  d.min = 1,
  nRuns = 1,
  shock.type = "full",
  unit.shock = amplitude,
  r.shock = 1.0,
  n.sources = N
)

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
  rCentMin = 0,
  rCentMax = N,
  symmetric = 1,
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

#### Running a Plan B centralization experiment on the ring 
# Note that shock type is "assets" whereassets are projects and 
# permutations are on the 'asset projects'.

N = 64
Total.External.Assets = 100
amplitude = Total.External.Assets / (2 * N)


connectivity.and.centralization(
  N,
  KMin = 1,
  KMax = 1,
  rCentMin = 0,
  rCentMax = N,
  symmetric = 1,
  eps = 0.1,
  phi = 0.4,
  A0 = Total.External.Assets,
  E0 = 10,
  b.sheet.equity = "hetero",
  d.min = 1,
  nRuns = 1000,
  shock.type = "assets",
  unit.shock = amplitude,
  r.shock = 1.0,
  n.sources = N
)

#### Running a Plan A centralization experiment on the ring
# where all assets of a selected bank is hit all at once.

N = 64
Total.External.Assets = 100
amplitude = Total.External.Assets / (2 * N)


connectivity.and.centralization(
  N,
  KMin = 1,
  KMax = 1,
  rCentMin = 0,
  rCentMax = N,
  symmetric = 1,
  eps = 0.1,
  phi = 0.4,
  A0 = Total.External.Assets,
  E0 = 10,
  b.sheet.equity = "planA.spec",
  d.min = 1,
  nRuns = 1000,
  shock.type = "full",
  unit.shock = amplitude,
  r.shock = 1.0,
  n.sources = N
)

#### Running a Plan A centralization experiment on the ring
# where all assets of a selected bank is hit incrementally.

N = 64
Total.External.Assets = 100
amplitude = Total.External.Assets / (2 * N)


connectivity.and.centralization(
  N,
  KMin = 1,
  KMax = 1,
  rCentMin = 0,
  rCentMax = N,
  symmetric = 1,
  eps = 0.1,
  phi = 0.4,
  A0 = Total.External.Assets,
  E0 = 10,
  b.sheet.equity = "planA.spec",
  d.min = 1,
  nRuns = 1000,
  shock.type = "banks",
  unit.shock = amplitude,
  r.shock = 1.0,
  n.sources = N
)


#### Running a Plan A centralization experiment on the ring
# where assets are permuted for the shock.

N = 64
Total.External.Assets = 100
amplitude = Total.External.Assets / (2 * N)


connectivity.and.centralization(
  N,
  KMin = 1,
  KMax = 1,
  rCentMin = 0,
  rCentMax = N,
  symmetric = 1,
  eps = 0.1,
  phi = 0.4,
  A0 = Total.External.Assets,
  E0 = 10,
  b.sheet.equity = "planA.spec",
  d.min = 1,
  nRuns = 1000,
  shock.type = "assets",
  unit.shock = amplitude,
  r.shock = 1.0,
  n.sources = N
)


#### Running a core-periphery experiment
# Note that while core is varied from a star to 4 core
# while the centralization is from 0 to maximum possible at
# each connectivity.

N = 64
Total.External.Assets = 100
amplitude = Total.External.Assets / (2 * N)


connectivity.and.centralization(
  N,
  KMin = 1,
  KMax = 4,
  rCentMin = 0,
  rCentMax = N,
  symmetric = 1,
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

