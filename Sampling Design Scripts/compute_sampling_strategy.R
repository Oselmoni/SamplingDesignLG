### This script guides to the choice of the sampling locations using one of the design proposed in the paper.


### Input file: the user must provide a .csv listing the landscape sites across which the sampling locations need to be selected
### The .csv must be organized as follows:
### 1st column: identifiers of landscape sites
### 2nd column: longitude of landscape sites
### 3rd column: latitutde of landscape sites
### 4th-Nth column: environmental variables of interest

### Dummy dataset

env = '../Data Landscape/ENV_landscape.csv' # replace the path to your csv here

### Load sampling functions from the R-pipeline

source('../R-pipeline/functions.R')

### Set the paramters of the Sampling Strategy you want to use.

size = 200 # total number of samples
locations = 20 # total number of sampling locations
design = 'hyb' # chose between 'geo' (geographic design), 'env' (environmental design), 'hyb' (hybrid design), 'ran' (random design)
k = 4 # number of environmental regions for hybrid design (to select this number, check the 'compute_k_for_hybrid.R' script)
corT = 0.5 # set a threshold for filtering correlated environmental variables in environmental and hybrid design


Sampling_Strategy = Sampling(Size = size, Sites = locations, Design = design, env = env, k = k, corT = corT)

Sampling_Strategy$CELLsampled # shows the id of the selected sampling locations 
