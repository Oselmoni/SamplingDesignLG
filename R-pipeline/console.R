source('functions.R')

### Load Env data

env= '../Data Landscape/ENV_landscape.csv' # csv file with environmental data: col1 = id, col2 = longitude, col3=latitude, col3-n = environmental variables
ps= '../Data Landscape/PS_landscape.csv' # csv file with population structure coefficient: col1=id, col2=ps coefficient

hyb_k = c('5'=2, '10'=3, '20'=5, '40'=5, '50'=5) # ks for hybrid design, calculated with the compute_k_for_hybrid.R script

### Simulations parameters
Size_levels = c(50,100,200,400,800,1600) # sample sizes
Sites_levels = c(5,10,20,40,50) # sampling locations
Design_levels= c('env','hyb','geo','ran') # sampling designs
c_levels=c(0.2, 0.5) # population structures
rep=1:20 # replicates
nbcore = 2 # nb of cpus to use for computation
nbsnps = 1000 # nb of SNPs per iteration
pathTOpython = 'python2.7' # path for calling python form terminal. (nb: required libraries: pandas, stastmodels, numpy )


### Run the iterations

nbruns = length(Sites_levels)*length(Size_levels)*length(Design_levels)*length(c_levels)*length(rep) # total number of iterations
cou=1  

line=paste(c('itnb', 'S', 'rS', 'E', 'D', 'c', 'r', 'TPR', 'FDR', 'nbSNPS', 'nbADA', 's', 'gamma', 'envADA'), collapse=' ')
write(line,file="output.txt")
for (S in Size_levels) { # for every level of sample size...

  for (E in Sites_levels) { # for every number of sampling locations...
    
    for (D in Design_levels) { # for every design approach...
      
      # Create sampling stategy
      
      if (D!='ran') { # non-random designs strategies do not vary 
        if (D=='env') { sa = Sampling(Size = S, Sites = E, Design = D, env = env, corT = 0.5) }
        else if (D=='hyb') {  sa = Sampling(Size = S, Sites = E, Design = D, env = env, k=hyb_k[as.character(E)], corT=0.5) }
        else if (D=='geo') { sa = Sampling(Size = S, Sites = E, Design = D, env = env) } 
        }
        
      for (c in c_levels) { # for every scenario of PS...
          
        for (r in rep) { #... replicate 20 times...
          
          if (D=='ran') {sa = Sampling(Size = S, Sites = E, Design = D, env = env)} # random design varies at every iteration
          
            print(paste0('Run #',cou,' out of ',nbruns))
        
            # compute genotype matrix
            gt = Genotype(SamplingStrategy = sa , nbSNPS=nbsnps, ps= ps, env= env, c=c, nbcore = nbcore, pathTOpy = pathTOpython)
            
            # run SAMBADA
            SamBada(GenotypeOutput = gt, nbcores=nbcore, pathTOpython = pathTOpython)  
            
            # Analyse SAMBADA output
            stats = PostProcessing(gt)
            
            line=paste(c(cou, S, sa$RealSize, E, D, c, r, stats$TPR, stats$FDR, gt$nbSNPS, gt$nbADA, gt$s, gt$gamma, gt$envADA), collapse=' ')
            write(line,file="output.txt",append=TRUE)
            
            cou=cou+1
          
}}}}}

