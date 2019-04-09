### Example on how ks for hybrid design can be calcualted out of environmental data

env = '../Data Landscape/ENV_landscape.csv' # set csv containing environmental data

library(NbClust) # load NbClust library

cells = read.csv(env, row.names = 1) # load csv

ie=names(redENV(cells[,3:ncol(cells)], 0.5)) # group environmental variables by clusters of correlated descriptors (R threshold = 0.5), choose one descriptor per cluster
enb = length(ie) # number of uncorrelated environmental descriptors

ENV = matrix(scale(cells[,ie]), ncol=enb) # keep only uncorrelated descriptors and scale them
rownames(ENV) = rownames(cells)

NbClust(ENV[,], min.nc=2, max.nc=10, method='complete') # check best number of cluster according to NbClust. In this examples, tests all possible values for k between 2 and 10. 
