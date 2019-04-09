## load results of the run

output = read.table('../Simulation Results/output.txt', header=T)

### Plot overall results

par(mfrow=c(1,2))

boxplot(output$TPR~output$E)
boxplot(output$FDR~output$E)

boxplot(output$TPR~output$D)
boxplot(output$FDR~output$D)

boxplot(output$TPR~output$S)
boxplot(output$FDR~output$S)

### ANOVA analysis

# Compare TPR and FDR under different pop. str. scenarios

anova(aov(output$TPR~as.factor(output$c)))
anova(aov(output$FDR~as.factor(output$c)))

# Panmictic

i=output$c==0.5

anova(aov(output$TPR[i]~as.factor(output$S[i])*as.factor(output$E[i])*as.factor(output$D[i])))
anova(aov(output$FDR[i]~as.factor(output$S[i])*as.factor(output$E[i])*as.factor(output$D[i])))

# Structured

i=output$c==0.2

anova(aov(output$TPR[i]~as.factor(output$S[i])*as.factor(output$E[i])*as.factor(output$D[i])))
anova(aov(output$FDR[i]~as.factor(output$S[i])*as.factor(output$E[i])*as.factor(output$D[i])))

### Function to get Means and Sds of FDR and TPR for subset of data

get_means_sd = function(output, subsets, rate) {
  
  ss =paste0('output$',paste(subsets, collapse='&output$'))
  rates = output[eval(parse(text=ss)),rate]
  
  return(paste0(round(mean(rates), 2), '±', round(sd(rates),2)))
  
}

## example of use: get FDR of geographic design under structured population (c=0.2)

get_means_sd(output, subsets=c('c==0.2' , 'D=="geo"'), rate='FDR')


######### Plot TPR and FDR by couples of factors

### Panmictic
par(mfrow=c(1,2))

i=output$c==0.5

colt = c(rgb(1,0,0, alpha=0.3) , rgb(0,1,0, alpha=0.3), rgb(0,0,1, alpha=0.3), rgb(1,1,0, alpha=0.3), rgb(0.5,0.5,0.5, alpha=0.3) )

### To show overall effects of Sampling Size - Sample Effort

boxplot(output$TPR[i]~output$S[i], main='TPR', outline=F )
boxplot(output$TPR[i&output$S==50]~as.numeric(output$E[i&output$S==50]), add=TRUE, col=colt, at=seq(0.8,1.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$TPR[i&output$S==100]~as.numeric(output$E[i&output$S==100]), add=TRUE, col=colt, at=seq(1.8,2.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==200]~as.numeric(output$E[i&output$S==200]), add=TRUE, col=colt, at=seq(2.8,3.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==400]~as.numeric(output$E[i&output$S==400]), add=TRUE, col=colt, at=seq(3.8,4.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==800]~as.numeric(output$E[i&output$S==800]), add=TRUE, col=colt, at=seq(4.8,5.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==1600]~as.numeric(output$E[i&output$S==1600]), add=TRUE, col=colt, at=seq(5.8,6.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

boxplot(output$FDR[i]~output$S[i], main='FDR', outline=F )
boxplot(output$FDR[i&output$S==50]~as.numeric(output$E[i&output$S==50]), add=TRUE, col=colt, at=seq(0.8,1.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$FDR[i&output$S==100]~as.numeric(output$E[i&output$S==100]), add=TRUE, col=colt, at=seq(1.8,2.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==200]~as.numeric(output$E[i&output$S==200]), add=TRUE, col=colt, at=seq(2.8,3.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==400]~as.numeric(output$E[i&output$S==400]), add=TRUE, col=colt, at=seq(3.8,4.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==800]~as.numeric(output$E[i&output$S==800]), add=TRUE, col=colt, at=seq(4.8,5.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==1600]~as.numeric(output$E[i&output$S==1600]), add=TRUE, col=colt, at=seq(5.8,6.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

### To show overall effects of Sampling Size - Sample Design

colt = c(rgb(1,0,0, alpha=0.3) , rgb(0,0,1, alpha=0.3), rgb(0,1,0, alpha=0.3), rgb(1,1,0, alpha=0.3) )

boxplot(output$TPR[i]~output$S[i], main='TPR', outline=F )
boxplot(output$TPR[i&output$S==50]~(output$D[i&output$S==50]), add=TRUE, col=colt, at=seq(0.7,1.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$TPR[i&output$S==100]~(output$D[i&output$S==100]), add=TRUE, col=colt, at=seq(1.7,2.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==200]~(output$D[i&output$S==200]), add=TRUE, col=colt, at=seq(2.7,3.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==400]~(output$D[i&output$S==400]), add=TRUE, col=colt, at=seq(3.7,4.3,0.2), boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==800]~(output$D[i&output$S==800]), add=TRUE, col=colt, at=seq(4.7,5.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==1600]~(output$D[i&output$S==1600]), add=TRUE, col=colt, at=seq(5.7,6.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

boxplot(output$FDR[i]~output$S[i], main='FDR', outline=F )
boxplot(output$FDR[i&output$S==50]~(output$D[i&output$S==50]), add=TRUE, col=colt, at=seq(0.7,1.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$FDR[i&output$S==100]~(output$D[i&output$S==100]), add=TRUE, col=colt, at=seq(1.7,2.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==200]~(output$D[i&output$S==200]), add=TRUE, col=colt, at=seq(2.7,3.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==400]~(output$D[i&output$S==400]), add=TRUE, col=colt, at=seq(3.7,4.3,0.2), boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==800]~(output$D[i&output$S==800]), add=TRUE, col=colt, at=seq(4.7,5.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==1600]~(output$D[i&output$S==1600]), add=TRUE, col=colt, at=seq(5.7,6.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

### To show overall effects of Sampling Effort - Sample Design


colt = c(rgb(1,0,0, alpha=0.3) , rgb(0,0,1, alpha=0.3), rgb(0,1,0, alpha=0.3), rgb(1,1,0, alpha=0.3) )

boxplot(output$TPR[i]~output$E[i], main='TPR', outline=F)
boxplot(output$TPR[i&output$E==5]~(output$D[i&output$E==5]), add=TRUE, col=colt, at=seq(0.7,1.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$TPR[i&output$E==10]~(output$D[i&output$E==10]), add=TRUE, col=colt, at=seq(1.7,2.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$E==20]~(output$D[i&output$E==20]), add=TRUE, col=colt, at=seq(2.7,3.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$E==40]~(output$D[i&output$E==40]), add=TRUE, col=colt, at=seq(3.7,4.3,0.2), boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$E==50]~(output$D[i&output$E==50]), add=TRUE, col=colt, at=seq(4.7,5.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

boxplot(output$FDR[i]~output$E[i], main='FDR', outline=F )
boxplot(output$FDR[i&output$E==5]~(output$D[i&output$E==5]), add=TRUE, col=colt, at=seq(0.7,1.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$FDR[i&output$E==10]~(output$D[i&output$E==10]), add=TRUE, col=colt, at=seq(1.7,2.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$E==20]~(output$D[i&output$E==20]), add=TRUE, col=colt, at=seq(2.7,3.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$E==40]~(output$D[i&output$E==40]), add=TRUE, col=colt, at=seq(3.7,4.3,0.2), boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$E==50]~(output$D[i&output$E==50]), add=TRUE, col=colt, at=seq(4.7,5.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')


### To show differences between Sampling efforts

boxplot(output$TPR[i]~output$E[i], main='TPR', outline=F)
boxplot(output$FDR[i]~output$E[i], main='FDR', outline=F)


### Structured
par(mfrow=c(1,2))

i=output$c==0.2

colt = c(rgb(1,0,0, alpha=0.3) , rgb(0,1,0, alpha=0.3), rgb(0,0,1, alpha=0.3), rgb(1,1,0, alpha=0.3), rgb(0.5,0.5,0.5, alpha=0.3) )

### To show overall effects of Sampling Size - Sample Effort

boxplot(output$TPR[i]~output$S[i], main='TPR', outline=F )
boxplot(output$TPR[i&output$S==50]~as.numeric(output$E[i&output$S==50]), add=TRUE, col=colt, at=seq(0.8,1.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$TPR[i&output$S==100]~as.numeric(output$E[i&output$S==100]), add=TRUE, col=colt, at=seq(1.8,2.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==200]~as.numeric(output$E[i&output$S==200]), add=TRUE, col=colt, at=seq(2.8,3.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==400]~as.numeric(output$E[i&output$S==400]), add=TRUE, col=colt, at=seq(3.8,4.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==800]~as.numeric(output$E[i&output$S==800]), add=TRUE, col=colt, at=seq(4.8,5.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==1600]~as.numeric(output$E[i&output$S==1600]), add=TRUE, col=colt, at=seq(5.8,6.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

boxplot(output$FDR[i]~output$S[i], main='FDR', outline=F )
boxplot(output$FDR[i&output$S==50]~as.numeric(output$E[i&output$S==50]), add=TRUE, col=colt, at=seq(0.8,1.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$FDR[i&output$S==100]~as.numeric(output$E[i&output$S==100]), add=TRUE, col=colt, at=seq(1.8,2.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==200]~as.numeric(output$E[i&output$S==200]), add=TRUE, col=colt, at=seq(2.8,3.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==400]~as.numeric(output$E[i&output$S==400]), add=TRUE, col=colt, at=seq(3.8,4.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==800]~as.numeric(output$E[i&output$S==800]), add=TRUE, col=colt, at=seq(4.8,5.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==1600]~as.numeric(output$E[i&output$S==1600]), add=TRUE, col=colt, at=seq(5.8,6.2,0.1) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

### To show overall effects of Sampling Size - Sample Design

colt = c(rgb(1,0,0, alpha=0.3) , rgb(0,0,1, alpha=0.3), rgb(0,1,0, alpha=0.3), rgb(1,1,0, alpha=0.3) )

boxplot(output$TPR[i]~output$S[i], main='TPR', outline=F )
boxplot(output$TPR[i&output$S==50]~(output$D[i&output$S==50]), add=TRUE, col=colt, at=seq(0.7,1.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$TPR[i&output$S==100]~(output$D[i&output$S==100]), add=TRUE, col=colt, at=seq(1.7,2.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==200]~(output$D[i&output$S==200]), add=TRUE, col=colt, at=seq(2.7,3.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==400]~(output$D[i&output$S==400]), add=TRUE, col=colt, at=seq(3.7,4.3,0.2), boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==800]~(output$D[i&output$S==800]), add=TRUE, col=colt, at=seq(4.7,5.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$S==1600]~(output$D[i&output$S==1600]), add=TRUE, col=colt, at=seq(5.7,6.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

boxplot(output$FDR[i]~output$S[i], main='FDR', outline=F )
boxplot(output$FDR[i&output$S==50]~(output$D[i&output$S==50]), add=TRUE, col=colt, at=seq(0.7,1.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$FDR[i&output$S==100]~(output$D[i&output$S==100]), add=TRUE, col=colt, at=seq(1.7,2.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==200]~(output$D[i&output$S==200]), add=TRUE, col=colt, at=seq(2.7,3.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==400]~(output$D[i&output$S==400]), add=TRUE, col=colt, at=seq(3.7,4.3,0.2), boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==800]~(output$D[i&output$S==800]), add=TRUE, col=colt, at=seq(4.7,5.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$S==1600]~(output$D[i&output$S==1600]), add=TRUE, col=colt, at=seq(5.7,6.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

### To show overall effects of Sampling Effort - Sample Design


colt = c(rgb(1,0,0, alpha=0.3) , rgb(0,0,1, alpha=0.3), rgb(0,1,0, alpha=0.3), rgb(1,1,0, alpha=0.3) )

boxplot(output$TPR[i]~output$E[i], main='TPR', outline=F)
boxplot(output$TPR[i&output$E==5]~(output$D[i&output$E==5]), add=TRUE, col=colt, at=seq(0.7,1.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$TPR[i&output$E==10]~(output$D[i&output$E==10]), add=TRUE, col=colt, at=seq(1.7,2.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$E==20]~(output$D[i&output$E==20]), add=TRUE, col=colt, at=seq(2.7,3.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$E==40]~(output$D[i&output$E==40]), add=TRUE, col=colt, at=seq(3.7,4.3,0.2), boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$TPR[i&output$E==50]~(output$D[i&output$E==50]), add=TRUE, col=colt, at=seq(4.7,5.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')

boxplot(output$FDR[i]~output$E[i], main='FDR', outline=F )
boxplot(output$FDR[i&output$E==5]~(output$D[i&output$E==5]), add=TRUE, col=colt, at=seq(0.7,1.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0 ,xaxt='n')
boxplot(output$FDR[i&output$E==10]~(output$D[i&output$E==10]), add=TRUE, col=colt, at=seq(1.7,2.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$E==20]~(output$D[i&output$E==20]), add=TRUE, col=colt, at=seq(2.7,3.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$E==40]~(output$D[i&output$E==40]), add=TRUE, col=colt, at=seq(3.7,4.3,0.2), boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')
boxplot(output$FDR[i&output$E==50]~(output$D[i&output$E==50]), add=TRUE, col=colt, at=seq(4.7,5.3,0.2) , boxwex=0.05, names=F, outline=F, whisklty = 0, staplelty = 0,xaxt='n')


### To show differences between Sampling efforts

boxplot(output$TPR[i]~output$E[i], main='TPR', outline=F)
boxplot(output$FDR[i]~output$E[i], main='FDR', outline=F)
