import os
os.environ["OMP_NUM_THREADS"] = "1"
import pandas as pd
import statsmodels.api as sm
import time
import warnings
import sys
warnings.filterwarnings("ignore")
import re

psGT=sys.argv[1]

print 'Opening environmental input'

chunknb = re.findall('sGTm_(.*).txt', psGT)[0]

### Check Mode: with or without PS


print 'Checking environmental variables'

ENV=pd.read_table('pySAM_v3/env-data.txt', sep=' ', header=0) # environmental input

try: PS=ENV.PS1 # check if population structure

except:

    # Mode without PS
    PS=pd.DataFrame()
    PS['PS1']= [1]*len(ENV.index)
    ENV=ENV.iloc[:,1:]

else:

    # Mode with PS
    lvar=(str(ENV.columns).split('u\''))
    nPS=0
    for i in lvar:
        if i[0:2]=='PS': nPS+=1

    PS =ENV.iloc[:,-nPS:]
    ENV=ENV.iloc[:,1:-nPS]



def toint(x): # function to convert string to list of numbers

    if (x.rstrip()== 'NA'):
        return float('NaN')
    else:
        return float(x.rstrip())

print 'Calculate Associations...'

### open writing stream to outputs

with open("pySAM_v3/NullLL_"+chunknb+".txt", "a") as nullLL:
  with open("pySAM_v3/AltLL_"+chunknb+".txt", "a") as altLL:
    with open("pySAM_v3/BetaALT_"+chunknb+".txt", "a") as betaALT:

## for each job (SNP) assigned to worker
     c=0
     with open(psGT) as inGTfile:

        for lineGT in inGTfile:

                c=c+1
            #    if (c%10==0):   print str(round(float(c)/float(int(sys.argv[2])),3))

                #################################################################################################################################################
                ## create a NULL MODEL  ####################################################################################################################
                #################################################################################################################################################

                GT = pd.Series([toint(x) for x in lineGT.split(' ')])

                marker = pd.to_numeric(GT)

                naMOL = pd.notnull(marker)

                naPS = pd.notnull(PS.PS1)


                for p in range(len(PS.columns)):
                    naPS = (naPS & PS.iloc[:, p])

                Marker = marker[naMOL & naPS]

                M0=((Marker==0)+0)
                M1=((Marker==1)+0)
                M2=((Marker==2)+0)

                ps = sm.add_constant(PS[naMOL & naPS])

                try:
                    logit_model = sm.Logit(M0, ps)
                    result = logit_model.fit(disp=0)

                except:
                    ll = 'NA'  # if error in model construction -> NA

                else:
                    ll = result.llf

                nullLL.write(str(ll)+'\n')

                try:
                    logit_model = sm.Logit(M1, ps)
                    result = logit_model.fit(disp=0)

                except:
                    ll = 'NA'  # if error in model construction -> NA

                else:
                    ll = result.llf

                nullLL.write(str(ll)+'\n')

                try:
                    logit_model = sm.Logit(M2, ps)
                    result = logit_model.fit(disp=0)

                except:
                    ll = 'NA'  # if error in model construction -> NA

                else:
                    ll = result.llf

                nullLL.write(str(ll)+'\n')

                ####################################################################################################################
                ### create ALT MODELS #######################################################################################
                ####################################################################################################################

                markLLs0=''
                b2s0=''
                markLLs1=''
                b2s1=''
                markLLs2=''
                b2s2=''

                for l in range(len(ENV.columns)):

                        E=ENV.iloc[:,l]

                        naE= pd.notnull(E)


                        Marker=marker[naMOL&naPS&naE]

                        M0=((Marker==0)+0)
                        M1=((Marker==1)+0)
                        M2=((Marker==2)+0)

                        ps=PS[naMOL&naPS&naE]
                        e=E[naMOL&naPS&naE]

                        exp = pd.concat([e, ps], axis=1, join_axes=[e.index])

                        exp=sm.add_constant(exp)

                        try:

                            logit_model = sm.Logit(M0, exp)
                            result=logit_model.fit(disp=0)


                        except:
                            markLLs0+='NA '
                            b2s0+='NA '

                        else:

                            ll=result.llf
                            markLLs0+=str(ll)+' '

                            b2=result.pvalues[ENV.columns[l]] # B of environmental variable
                            b2s0+=str(b2)+' '

                        try:

                            logit_model = sm.Logit(M1, exp)
                            result=logit_model.fit(disp=0)


                        except:
                            markLLs1+='NA '
                            b2s1+='NA '

                        else:

                            ll=result.llf
                            markLLs1+=str(ll)+' '

                            b2=result.pvalues[ENV.columns[l]] # B of environmental variable
                            b2s1+=str(b2)+' '


                        try:

                            logit_model = sm.Logit(M2, exp)
                            result=logit_model.fit(disp=0)


                        except:
                            markLLs2+='NA '
                            b2s2+='NA '

                        else:

                            ll=result.llf
                            markLLs2+=str(ll)+' '

                            b2=result.pvalues[ENV.columns[l]] # B of environmental variable
                            b2s2+=str(b2)+' '



                altLL.write(markLLs0[:-1]+'\n')
                betaALT.write(b2s0[:-1]+'\n')
                altLL.write(markLLs1[:-1]+'\n')
                betaALT.write(b2s1[:-1]+'\n')
                altLL.write(markLLs2[:-1]+'\n')
                betaALT.write(b2s2[:-1]+'\n')

print('Done')

