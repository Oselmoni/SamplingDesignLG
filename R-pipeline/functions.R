library(qvalue)
source('auxiliary_functions.R')

### Function to compute sampling strategies


Sampling = function(Size, Sites, Design, env, k, corT) { # k = nb of environmental regions for hybrid design
  
  print('...Read Input...')

  cells=read.csv(env,row.names=1)
  
  par(mfrow=c(1,1))
  
  if (Design=='geo') {
    DM = dist(cells[,1:2])
    SCL = cutree(hclust(DM), Sites)
    mlong= by(cells[,1],SCL,mean)
    mlat=by(cells[,2],SCL,mean)
    
    CELLsampled=c()
    
    for (i in 1:Sites) {
      
      dis= sqrt((mlong[i]-cells[SCL==i,1])^2+(mlat[i]-cells[SCL==i,2])^2)
      CELLsampled= c(CELLsampled,rownames(cells)[which(SCL==i)[which.min(dis)]])
      
    }
    
    CELLsampled=as.character(CELLsampled)
    
    plot(cells[,1:2], pch='.', cex=1, col=SCL)
    points(cells[CELLsampled,1:2], pch=16, col=SCL[CELLsampled])
  } 
  if (Design=='env') {
    
    ie=names(redENV(cells[,3:ncol(cells)], corT))
    enb = length(ie)
    
    ENV = matrix(scale(cells[,ie]), ncol=enb)
    rownames(ENV) = rownames(cells)
    
    
    DM = dist(ENV)
    SCL = cutree(hclust(DM), Sites)
    
    CELLsampled=c()
    
    for (i in 1:Sites) {
      
      scells=matrix(ENV[SCL==i,], ncol=enb)
      rownames(scells) = names(which(SCL==i))
      
      if (nrow(scells)==1) { CELLsampled=c(CELLsampled, names(which(SCL==i))) 
      } else {
        
        mENV=apply(scells,2,mean)
        
        dis=as.matrix(dist(rbind(mENV,scells)))[1,2:nrow(scells)]
        
        CELLsampled= c(CELLsampled,rownames(scells)[which.min(dis)])
      }
    }
    
    plot(cells[,1:2], col=SCL, pch='.', cex=1)
    points(cells[CELLsampled,1:2], pch=16, cex=1, col=SCL[CELLsampled]) 
    
  } 
  if (Design=='ran') {
    
    CELLsampled = as.character(sample(rownames(cells), size = Sites, F))
    plot(cells[,1:2], pch='.', cex=1)
    points(cells[CELLsampled,1:2], pch=16, cex=1)
    
  }
  if (Design=='hyb'){

    ie=names(redENV(cells[,3:ncol(cells)], corT))
    enb = length(ie)
    
    ENV = matrix(scale(cells[,ie]), ncol=enb)
    rownames(ENV) = rownames(cells)
    
    DM = dist(ENV)
    ECL = cutree(hclust(DM), k)
    
    CELLsampled=c()
    
    sk = table(1+(1:Sites%%k)) # nb of sites per environmental cluster
    
    for (i in 1:length(unique(ECL))) {
      
      scells=cells[ECL==i,]
      
      if (nrow(scells)< sk[i]) {
        
        ind = ((1:(sk[i]))%%nrow(scells))+1 # if there are not enough cells in a cluster...
        CELLsampled = c(CELLsampled, rownames(scells)[ind]) # take (and re-take) those available
        
      } else {
        
        DM = dist(scells[,1:2])
        
        SCL = cutree(hclust(DM), sk[i])

        for (l in 1:length(unique(SCL))) {
          
          sscells=as.matrix(cells[names(which(SCL==l)),1:2])
          rownames(sscells) = names(which(SCL==l))
          
          if (nrow(sscells)==1) { CELLsampled=c(CELLsampled, names(which(SCL==l))) 
          } else {
            
            mlong= apply(sscells,2,mean)[1]
            mlat=apply(sscells,2,mean)[2]
            
            dis= sqrt((mlong-sscells[,1])^2+(mlat-sscells[,2])^2)
            CELLsampled= c(CELLsampled,rownames(sscells)[which.min(dis)])
            
          }
        }
      }
    } 
    plot(cells[,1:2], col=ECL, pch='.', cex=1)
    points(cells[CELLsampled,1:2], pch=16, col=as.numeric(rep(names(sk), times=sk )))
  }  
  

  IndPerSites = round(Size/Sites)
  
  return(list('CELLsampled'=CELLsampled, 'Sites'=Sites, 'IndPerSites'=IndPerSites, 'Size'=Size, 'RealSize'=Sites*IndPerSites,'Design'=Design))

}  

### Function to compute genotype matrix

Genotype = function(SamplingStrategy, nbSNPS=1000, c=0.2, rangeADA=5:20, rangeSEL=seq(0.1,0.3,0.01), rangeGAMMA=seq(0,0.5,0.01), ps, env, nbcore=1, pathTOpy='python2.7') {
  
    nbADA = sample(rangeADA,1)
    nbNEU = nbSNPS-nbADA
    
  
    print('Compute Neutral...')
    
    cells_ps = read.csv(ps, row.names=1)
    scells_ps = cells_ps[SamplingStrategy$CELLsampled,]
    
    pA_ps = LIN(scells_ps, cells_ps, c) # probability of allele knowning ps_membership
    
    write.table(pA_ps,'Genotype_Py/pA_ps.txt', col.names=F, row.names=F) # input for python script: pA_ps
    
    ### Create bash script to parallelize computation of neutral gt matrix

    commander = c('#!/usr/bin/env bash','set -e', 'rm -f Genotype_Py/GT*')
    if (nbcore == 1) {nbchunks = c("all"=nbNEU*SamplingStrategy$IndPerSites)} else {nbchunks = table(cut(1:(nbNEU*SamplingStrategy$IndPerSites), nbcore))}
    co=1
    for (i in nbchunks) {
      
      commander = c(commander, paste0(pathTOpy, ' Genotype_Py/compute_gt.py ',i,' ',co,' &'))
      
      co=co+1
    }
    commander = c(commander, 'wait', 'cat Genotype_Py/GT_* > Genotype_Py/GT.txt','echo Done!')
    
    write.table(commander, 'Genotype_Py/gt_commander.sh', col.names=F, row.names=F, quote=F)
    
    system('chmod 755 Genotype_Py/gt_commander.sh')
    system('./Genotype_Py/gt_commander.sh')
    
    GTneu = read.table('Genotype_Py/GT.txt', sep=',') # read output of python script
    dim(GTneu)
    
    GTneu = matrix(as.matrix(GTneu), ncol=(SamplingStrategy$IndPerSites*SamplingStrategy$Sites))
    
    
    #################################
    
    print('Compute Adaptive...')
  
    cells_env = read.csv(env, row.names=1)
    scells_env = cells_env[SamplingStrategy[[1]],]
    
    envada = sample(colnames(cells_env)[3:ncol(cells_env)], size = nbADA, T) # select adaptive forces
    s = sample(rangeSEL, size = nbADA, T) # sample selection strength
    gamma = sample(rangeGAMMA, size=nbADA, T) # gamma parameters: contribution of pA_ps to pA. 
  
    pA=c() # probability of allele knowning environmental variable (selection)
    for (i in 1:nbADA) {
      pA_env = LIN(scells_env[envada[i]], cells_env[envada[i]], s[i])
      pA = cbind(pA,   (pA_env*(1-gamma[i]))+(pA_ps*gamma[i]) )
    }
    
  
    GTada=c()
    for (i in 1:(SamplingStrategy[[3]])) { # compute GTs out of allele probabilities
      
      GTada=rbind(GTada, apply(pA, 1, function(x) { return(getGT(x))} ) )
      
    }
  
  
    GTada = matrix(GTada, ncol=(SamplingStrategy[[3]]*SamplingStrategy[[2]]))
    
    GT=t(rbind(GTada, GTneu)) # final gt table
    colnames(GT) = paste0('snps', 1:nbSNPS)
    rownames(GT) = paste0( rep(SamplingStrategy[[1]], each=SamplingStrategy[[3]]), '_', 1:SamplingStrategy[[3]])
    
    ENV=scells_env[rep(SamplingStrategy[[1]], each=SamplingStrategy[[3]]),]
    rownames(ENV) = rownames(GT)
    ENV=ENV[,-c(1:2)]
    
  
    trues = cbind(colnames(GT)[1:nbADA], envada)
    
    return(list('env'=ENV, 'gt'=GT, 'trues'=trues, 'nbSNPS'=nbSNPS, 'c'=c, 'nbADA'=nbADA, 's'=paste(as.character(s), collapse='_'), 'gamma'=paste(as.character(gamma), collapse='_'), 'envADA'=paste(as.character(envada), collapse='.')))
  
  }

### Function to prepare genotype matrix for analysis and launch SAMBADA (pySAM_v3)

SamBada = function(GenotypeOutput, nbcores=1, pathTOpython = 'python2.7') {
  
  print('...Process Genotype Matrix')
  
  ## Load Data
  
  fSNPS = GenotypeOutput$gt
  
  env = GenotypeOutput$env
  
  ### Minor Allele Frequency
  
  ## For each SNP, it computes the frequncy of the most rare allele across all individuals. 
  
  MAF <- apply(fSNPS, 2, function(x) {
    a=(sum(x=='1', na.rm=T)+sum(x=='0', na.rm=T)*2)/(sum(is.na(x)==F)*2)
    A=(sum(x=='1', na.rm=T)+sum(x=='2', na.rm=T)*2)/(sum(is.na(x)==F)*2)
    return(min(a,A))
  })
  fSNPS <- fSNPS[,MAF>0.05] # We apply a 5% threshold. 
  
  ### Major Genotype Frequency
  
  ## For each SNP, it computes the frequency of the most frequent genotype across all individuals. 
  
  MGF <- apply(fSNPS, 2, function(x) {
    aa=(sum(x=='0', na.rm=T))/(sum(is.na(x)==F))
    aA=(sum(x=='1', na.rm=T))/(sum(is.na(x)==F))
    AA=(sum(x=='2', na.rm=T))/(sum(is.na(x)==F))
    return(max(aa,aA,AA))
  })
  fSNPS <- fSNPS[,MGF<0.95] # We apply a 95% threshold 
  
  ##########################  Population Genetic Structure

  PCA= prcomp(fSNPS, center = T, scale=T)
  
  par(mfrow=c(1,1))
  #plot(PCA$sdev/sum(PCA$sdev), pch=16, main=' % the Explained Variation', ylab='% the Explained Variation', xlab='PC#')
  
  DIFFPC = diff(c(0,cumsum(PCA$sdev)/sum(PCA$sdev)))
  
  PSPC = NbPCs(DIFFPC,0.2) # PSPC is the number of PC suggested to be used to quantify the population structure. 
  
  if (PSPC > 0) { env$PS1 = PCA$x[,1]} # it there is a structure, use first PC in SamBada
  
  env= cbind('Name'=rownames(env), env)

  ######################### Write SamBada (PySAMv3) Inputs
  
  print('Run SamBada')
  
  write.table(env, 'pySAM_v3/env-data.txt', row.names = F, quote=F) # environmental input
  
  sGTm = t(fSNPS) # genetic input

  if (nbcores == 1) {chunks=rep(1, times=nrow(sGTm))} else {chunks= as.numeric(cut(1:nrow(sGTm), breaks=nbcores))} ## index by computational chunk

  bashscript =c('#!/usr/bin/env bash', 'set -e', 'rm -f pySAM_v3/AltLL* pySAM_v3/NullLL* pySAM_v3/BetaALT*')
  
  catlineN = 'cat '
  catlineA = 'cat '
  catlineB = 'cat '
  
  #### Create bash script to parallelize pysam_v3
  
  for (i in 1:nbcores) {
    
    write.table(sGTm[chunks==i,], paste0('pySAM_v3/sGTm_',i,'.txt'), col.names=F, row.names=F, quote=F)
    
    bashscript = c(bashscript, paste0(pathTOpython, ' pySAM_v3/pySAM_v3.py pySAM_v3/sGTm_',i,'.txt ', sum(chunks==i),' &'))
    
    catlineN = paste0(catlineN, 'pySAM_v3/NullLL_',i,'.txt ')
    catlineA = paste0(catlineA, 'pySAM_v3/AltLL_',i,'.txt ')
    catlineB = paste0(catlineB, 'pySAM_v3/BetaALT_',i,'.txt ')
  }
  
  bashscript = c(bashscript, 'wait', paste0(catlineN, '> pySAM_v3/NullLL.txt &'), 
                 paste0(catlineA, '> pySAM_v3/AltLL.txt &'), 
                 paste0(catlineB, '> pySAM_v3/BetaALT.txt &'), 'wait', 'echo Done!')
  
  write.table(bashscript, 'pySAM_v3/pySAM_commander.sh', col.names=F, quote=F, row.names=F)
  
  system('chmod 755 pySAM_v3/pySAM_commander.sh')
  t1=Sys.time()
  system('./pySAM_v3/pySAM_commander.sh')
  t2=Sys.time()
  print(paste0('SamBada Run lasted ',as.character(round(t2-t1, 2)),' seconds!'))
}  

### Function to read Sambada Output and Compute Statistics of Performance

PostProcessing = function(GenotypeOutput) {
  
  nullLL=read.table('pySAM_v3/NullLL.txt')
  altLL=read.table('pySAM_v3/AltLL.txt')
  pvaluesB=read.table('pySAM_v3/BetaALT.txt')
  ENV=read.table('pySAM_v3/env-data.txt', header=T)  
  
  PySAM = processPySam(nullLL = nullLL, altLL = altLL, pvaluesB = pvaluesB, ENV=ENV, GT = GenotypeOutput$gt, pg = 0.05, pb=0.05)
  
  positives = paste0(   substr(PySAM$sign$Marker,1,nchar(as.character(PySAM$sign$Marker))-2) , '__', PySAM$sign$Env_1)
  trues = paste0(GenotypeOutput$trues[,1],'__',GenotypeOutput$trues[,2])

  TPR = sum(positives%in%trues) / length(trues)

  FDR = 1- (sum(trues%in%positives) / length(positives))
  
  return(list('TPR'=TPR, 'FDR'=FDR))
  
}

