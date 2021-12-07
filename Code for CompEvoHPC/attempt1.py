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
import pfix_sim_functions as pfix_sim  # my functions in a seperate file 
import fig_functions as myfun  
import math as math
######Libre#######
######Datum#######

T = 1e9         # carrying capacity
b=2.0           # juvenile birth rate
do=100/98.0     # minimum death rate / death rate of optimal genotype
di=np.arange(start=3, stop=1.0, step=-0.1)
#b = np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.82, 1.748, 1.676, 1.604, 1.532, 1.46, 1.388, 1.316]);

sa = 1e-2       # selection coefficient of beneficial mutation in 
sr = 0.175      # multiplicative increment to "c" is (1+sr)
yi_option = 3   # select option to get equilibrium population density (=1,2,3)
max_index = 400 #
de = b+1        # extinction death rate
sa = 1e-2       # selection coefficient of beneficial mutation in 
sr = 0.175      # multiplicative increment to "c" is (1+sr)
yi_option = 3   # option for determining equilibrium density
samp = 1000       # 


t = time.time()
element_run_eqpop = Parallel(n_jobs=-1)(delayed(myfun.get_eq_pop_density)(b,di[k],sa,yi_option) for k in range(len(di)-1))
print(time.time() - t)  

parallel_eqpop = open("parallel_eqpop.txt","a")
for row in element_run_eqpop:
    parallel_eqpop.write(str(int(math.ceil(T*row[0])))+'\n')
parallel_eqpop.close()

parallel_effsr = open("parallel_effsr.txt","a")
for row in element_run_eqpop:
    parallel_effsr.write(str(myfun.get_c_selection_coefficient(b,row[0],sr))+'\n')
parallel_effsr.close()
#The idea would be to take this and apply to
d_Inc = 1
c_Inc = 0

t = time.time()
element_run_abs = Parallel(n_jobs=-1)(delayed(pfix_sim.modsimpop)(d_Inc,c_Inc,samp,T,sr,b,di[k:(k+2)],do,sa,de,yi_option) for k in range(len(di)-1))
print(time.time() - t)  

parallel_abs = open("parallel_abs.txt","a")
for row in element_run_abs:
    parallel_abs.write(str(row)+'\n')
parallel_abs.close()

d_Inc = 0
c_Inc = 1
#d_pfixes[i] = pfix_sim.modsimpop(d_Inc,c_Inc,samp,T,sr,b,di[i:(i+2)],do,sa,de,yi_option) #fix
t = time.time()
element_run_rel = Parallel(n_jobs=-1)(delayed(pfix_sim.modsimpop)(d_Inc,c_Inc,samp,T,sr,b,di[k:(k+2)],do,sa,de,yi_option) for k in range(len(di)-1))
print(time.time() - t)  

parallel_rel = open("parallel_rel.txt","a")
for row in element_run_rel:
    parallel_rel.write(str(row)+'\n')
parallel_rel.close()