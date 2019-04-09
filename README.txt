README

The supplementary material of the paper on sampling design in landscape genomics is organized in four folder.

If you want to create a sampling strategy for your data, use the scripts contained in the Sampling Design Scripts folder:

- Sampling Deisgn Scripts: contains scripts for creating a sampling strategy out of a user-defined input. 
> 'compute_sampling_strategy.R': the R script to compute sampling strategy. 
> 'compute_k_for_hybrid.R': the R script to compute k for the hybrid design using NbClust. 



If you want to reproduce the analysis of the paper, use the codes and data avbailable in the other folders.

- Data Landscape: contains the dataset used for the simulations of the paper. 
> 'ENV_landscape.csv' : table containing coordinates and enviromental values for every landscape site of Europe.
> 'PS_lansdscape.csv': table containing the population membership coefficient for every landscaep of europe. 
> 'Shapefile Grid Landscape': shapefile to visualize environmental data on a GIS. 
> 'PS_europe.png': geo-referenced digital image of population structure gradient. 


- R-pipeline: contains all the functions and codes for running the simulations. Note that the pipeline is constructed for a UNIX-derived file system. 
> 'console.R' : the central script of the pipeline. Runs the iterations.
> 'functions.R': contains the four main functions to run the pipeline. 
> 'auxiliary_functions.R': contains secondary function that are also requirred for the pipeline. 
> 'evaluate_output.R' : reads the output of the simulations and produces the plots contained in the publication. 
> 'Genotype_Py': contains the python script for the computation of genotypes called during the pipeline.
> 'pySAM_v3': contains the python script for the computation of landscape genomics analysis called during the pipeline.

- Simulation Results: contains the results of the simulation used for the paper. 
> 'output.txt' : is the result of the simulations ran for the paper. 
