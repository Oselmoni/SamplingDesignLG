### Function to transform allelic frequencies to genotypes

getGT = function(pAs) {
  
  GT=c()
  for (ind in pAs) {  
    Alleles=sample(c(1,0),2,T,c(ind, 1-ind))
    Alleles
    GT = c(GT, sum(Alleles))
  }
  return(GT)
}

### Function to perform a linear transformation of a values (e.g. membership coefficient or environmental variable) to a probability between 0 and 1. 

LIN = function(x,range,coef) {
  # 'x' is the value to transform, 'range' the range from which this value is observed and 'coef' the paramter defining the transformation (c for membership coefficient, s for environmental variables)
  slo=(1-2*coef)/(max(range)-min(range))
  ori=1-coef-slo*max(range)
  
  return((as.matrix(x)*slo)+ori)
}

### Function to calculate number of PCs suggested to describe Pop Str

NbPCs = function(DIFFPC, cutoff=0.1) {
  
  co=1
  while (((DIFFPC[co]-DIFFPC[co+1])/DIFFPC[co])>cutoff) {
    co=co+1}
  
  print(paste0('Number of PC suggested for describing pop structure: ',co-1))
  return(co-1)
  
}

### Function to process outpuit of pySAM_v3

processPySam = function(nullLL=nulllLL, altLL=altLL, pvaluesB=pvaluesB, ENV=ENV, GT, pg=0.01, pb=0.01) {
  
  rownames(nullLL) = paste0(rep(colnames(GT), each=3), c('_0', '_1', '_2'))
  rownames(altLL) = paste0(rep(colnames(GT), each=3), c('_0', '_1', '_2'))
  rownames(pvaluesB) = paste0(rep(colnames(GT), each=3), c('_0', '_1', '_2'))
  
  
  Gscores=apply(altLL,2, function(x) { return( -2* (as.numeric(nullLL$V1)-as.numeric(x)) ) })
  rownames(Gscores)=rownames(nullLL)
  
  pvaluesG=apply(Gscores,2, function(x) { return(1-pchisq(x, df=1))})
  
  qvaluesG=apply(pvaluesG, 2, function(x) { return(qvalue(x)$qvalue)})
  
  qvaluesB=apply(pvaluesB, 2, function(x) { return(qvalue(x)$qvalue)})
  
  res=c()
  
  for (i in (1:ncol(qvaluesG))) {

      signmar=names(which(qvaluesG[,i]<pg&qvaluesB[,i]<pb))
    if (length(signmar)>0) {  res=rbind(res,cbind(  colnames(ENV)[i+1], signmar ,qvaluesG[signmar,i],qvaluesB[signmar,i]))}
    
  }
  
  rownames(res) = NULL
  
  best=c()
  
  genes=substr(res[,2], 1, nchar(res[,2])-2)
  
  for (un in unique(genes)) {
    
    subres=res[genes==un,]
    
    if (is.null(dim(subres))) {best=rbind(best, subres)
    }else { best=rbind(best,  subres[order(as.numeric(subres[,3]), as.numeric(subres[,4]))[1],])}
    
  }
  
  sign=data.frame(cbind("Marker"=best[,2], "Env_1"=best[,1]), row.names = c(substr(best[,2], 1, nchar(best[,2])-2)))
  
  return(list("sign"=sign, "best"=best, "res"=res, "topQG"=sort(qvaluesG), "topQW"=sort(qvaluesB), 'pvaluesG'=pvaluesG))
  
}

### Function to compute groups of correlated environmental variables

redENV = function(ienv, corcutoff=0.7) {
  
  listvar = colnames(ienv)
  
  co=1
  
  olist=list()
  
  while(length(listvar)>0) {
    
    C = abs(cor(ienv)[listvar[co],listvar[-co]]) 
    group=names(C[C>corcutoff]) 
    
    olist[listvar[co]]=paste(group, collapse = ', ')
    
    listvar=listvar[-c(1,which(listvar%in%group))]
    
  }
  
  return(olist)
  
}

