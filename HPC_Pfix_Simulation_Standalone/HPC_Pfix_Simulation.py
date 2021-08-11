# -*- coding: utf-8 -*-
"""
Created on Mon May 10 13:21:46 2021
@author: wyseq
"""

#import time
#from joblib import Parallel, delayed
#
#def countdown(n):
#    while n>0:
#        n -= 1
#    return n
#
#t = time.time()
#for _ in range(20):
#    print(countdown(10**7))
#    
#print(time.time() - t)  
## takes ~10.5 seconds on medium sized Macbook Pro
#
#t = time.time()
#results = Parallel(n_jobs=2)(delayed(countdown)(10**7) for _ in range(20))
#print(results)
#print(time.time() - t)
#
##def yourfunction(k):   
##    s=3.14*k*k
##    print "Area of a circle with a radius ", k, " is:", s
##element_run = Parallel(n_jobs=1)(delayed(yourfunction)(k) for k in range(1,10))
##n_jobs=-1: use all available cores

import time
from joblib import Parallel, delayed
import numpy as np
import fig_functions as myfun  
import math as math
import csv


parName = {'T':0,'b':1,'do':2,'sa':3,'UaMax':4,'Uad':5,'cr':6,'Ur':7,'Urd':8,'R':9,'part':10}
parValue = np.zeros([11])

with open('evoParameters.csv','r') as csvfile:
    parInput = csv.reader(csvfile)
    for (line,row) in enumerate(parInput):
        parValue[line] = float(row[0])

# calculate list of death terms, i_ext, and upper-bound on death term size
# for population to be viable
[di,i_ext,d_max] = myfun.get_death_terms_classes(parValue[parName['b']], \
                            parValue[parName['do']], parValue[parName['sa']])


# Adjusted all index terms to be at increments of 17/18 steps in the absolute fitness class ~10 samples
#T = 1e9         # carrying capacity
##do=100/98.0     # minimum death rate / death rate of optimal genotype
## di=np.arange(start=3, stop=1.0, step=-0.1)
# b = np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.82, 1.748, 1.676, 1.604, 1.532, 1.46, 1.388, 1.316]);
# di = np.flip(np.array([1.020408163,1.020824656 ,1.021249823 ,1.02168385 ,1.022126931 ,1.022579262 ,1.023041044 ,1.023512482 ,1.023993789 ,1.024485178 ,1.024986872 ,1.025499097 ,1.026022083 ,1.026556067 ,1.027101293 ,1.027658009 ,1.028226468 ,1.028806932 ,1.029399668 ,1.030004948 ,1.030623053 ,1.031254269 ,1.031898891 ,1.03255722 ,1.033229564 ,1.033916239 ,1.03461757 ,1.035333889 ,1.036065536 ,1.036812862 ,1.037576223 ,1.038355986 ,1.03915253 ,1.039966239 ,1.04079751 ,1.041646749 ,1.042514373 ,1.043400809 ,1.044306498 ,1.045231889 ,1.046177446 ,1.047143642 ,1.048130965 ,1.049139916 ,1.050171009 ,1.051224772 ,1.052301747 ,1.053402491 ,1.054527578 ,1.055677594 ,1.056853146 ,1.058054855 ,1.059283359 ,1.060539317 ,1.061823403 ,1.063136314 ,1.064478764 ,1.065851489 ,1.067255248 ,1.068690818 ,1.070159003 ,1.071660629 ,1.073196546 ,1.074767632 ,1.076374788 ,1.078018946 ,1.079701064 ,1.081422131 ,1.083183165 ,1.084985217 ,1.086829371 ,1.088716745 ,1.090648493 ,1.092625806 ,1.094649913 ,1.096722083 ,1.098843628 ,1.101015902 ,1.103240304 ,1.105518282 ,1.107851329 ,1.110240994 ,1.112688876 ,1.115196629 ,1.117765967 ,1.120398662 ,1.123096552 ,1.125861539 ,1.128695592 ,1.131600755 ,1.134579145 ,1.137632959 ,1.140764475 ,1.143976057 ,1.14727016 ,1.150649333 ,1.154116225 ,1.157673585 ,1.161324276 ,1.165071272 ,1.168917668 ,1.172866685 ,1.176921677 ,1.181086136 ,1.185363702 ,1.18975817 ,1.194273497 ,1.198913811 ,1.203683421 ,1.208586828 ,1.213628734 ,1.218814053 ,1.224147926 ,1.229635731 ,1.235283097 ,1.241095921 ,1.247080385 ,1.253242967 ,1.259590466 ,1.266130019 ,1.272869124 ,1.279815657 ,1.286977906 ,1.294364591 ,1.301984893 ,1.309848488 ,1.31796558 ,1.326346934 ,1.335003919 ,1.34394855 ,1.353193533 ,1.362752317 ,1.372639148 ,1.38286913 ,1.393458288 ,1.404423642 ,1.415783284 ,1.427556465 ,1.439763685 ,1.4524268 ,1.465569136 ,1.479215612 ,1.493392876 ,1.508129464 ,1.523455964 ,1.539405207 ,1.55601247 ,1.573315717 ,1.59135585 ,1.610177001 ,1.629826861 ,1.650357036 ,1.671823462 ,1.694286866 ,1.717813289 ,1.742474673 ,1.768349539 ,1.79552375 ,1.824091386 ,1.854155746 ,1.885830501 ,1.919241025 ,1.954525927 ,1.99183884 ,2.031350503 ,2.07325119 ,2.117753576 ,2.165096109 ,2.21554701 ,2.269409041 ,2.327025208 ,2.38878563 ,2.455135853 ,2.526586977 ,2.603728072 ,2.687241508 ,2.777922017 ,2.876700591 ,2.984674705]))
#sa = 1e-2       # selection coefficient of beneficial mutation in 
#sr = 0.175/3.0  # multiplicative increment to "c" is (1+sr)
#yi_option = 3   # select option to get equilibrium population density (=1,2,3)
#max_index = 400 #
#de = b+1        # extinction death rate
yi_option = 3   # option for determining equilibrium density
samp = 10      # 
part = int(parValue[parName['part']])
subsamp = len(di)/part - 1
for k in range(subsamp):
    print(di[k*part:(k*part + 2)])



t = time.time()
element_run_eqpop = Parallel(n_jobs=1)(delayed(myfun.get_eq_pop_density)(parValue[parName['b']],di[k*part],parValue[parName['sa']],yi_option) for k in range(subsamp)) #got rid of the minus 1
print(time.time() - t)  

parallel_eqpop = open("parallel_eqpop.txt","a")
for row in element_run_eqpop:
    parallel_eqpop.write(str(int(math.ceil(parValue[parName['T']]*row[0])))+'\n')
parallel_eqpop.close()

parallel_effsr = open("parallel_effsr.txt","a") 
for row in element_run_eqpop: # b,di,y,cr
    parallel_effsr.write(str(myfun.get_c_selection_coefficient_OLD(parValue[parName['b']],row[0],parValue[parName['cr']]))+'\n')
parallel_effsr.close()
#The idea would be to take this and apply to
d_Inc = 1
c_Inc = 0

t = time.time()
element_run_abs = Parallel(n_jobs=1)(delayed(myfun.modsimpop)(d_Inc,c_Inc,samp,parValue[parName['T']],parValue[parName['cr']],parValue[parName['b']],di[k*part:(k*part + 2)],parValue[parName['do']],((di[k*part]/di[k*part + 1])-1)/(di[k*part + 1]-1),d_max,yi_option) for k in range(subsamp)) #got rid of the minus 1
print(time.time() - t)  

parallel_abs = open("parallel_abs.txt","a")
for row in element_run_abs:
    parallel_abs.write(str(row)+'\n')
parallel_abs.close()

d_Inc = 0
c_Inc = 1
#d_pfixes[i] = myfun.modsimpop(d_Inc,c_Inc,samp,parValue[parName['T']],parValue[parName['cr']],parValue[parName['b']],di[i:(i+2)],parValue[parName['do']],parValue[parName['cr']],d_max,yi_option) #fix
t = time.time()
element_run_rel = Parallel(n_jobs=1)(delayed(myfun.modsimpop)(d_Inc,c_Inc,samp,parValue[parName['T']],parValue[parName['cr']],parValue[parName['b']],di[k*5:(k*5 + 2)],parValue[parName['do']],((di[k*5]/di[k*5 + 1])-1)/(di[k*5 + 1]-1),d_max,yi_option) for k in range(subsamp)) #got rid of the minus 1
print(time.time() - t)  

parallel_rel = open("parallel_rel.txt","a")
for row in element_run_rel:
    parallel_rel.write(str(row)+'\n')
parallel_rel.close()
