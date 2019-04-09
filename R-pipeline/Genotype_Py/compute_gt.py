import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy
import sys

gtXind = int(sys.argv[1])
chunk = sys.argv[2]

PAps= open('Genotype_Py/pA_ps.txt', 'r').readlines()


def getGT(pAs):
    GT=[]
    for ind in pAs:

        p=float(ind.rstrip())

        Alleles = numpy.random.choice([1,0], size=2, replace=True, p=[p, 1-p])

        GT+=[str(sum(Alleles))]

    return(GT)


outGT=[]

for i in range(gtXind):

    outGT+=[','.join(getGT(PAps))+'\n']


o=open('Genotype_Py/GT_'+str(chunk)+'.txt', 'w')
o.writelines(outGT)
o.close()

