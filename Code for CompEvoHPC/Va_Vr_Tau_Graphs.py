# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 08:33:25 2021

@author: wyseq
"""

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import fig_functions as myfun            # my functions in a seperate file
import csv
import math as math
# load parameters from evoParameters.csv file. The parameters must be provided in the following 
# order (numeric values). All input csv files must be available in working directory containing scripts.
#T   = 1e9              # carrying capacity
#b   = 2.0              # birth rate
#do = 100/98.0          # minimum death rate / death rate of optimal genotype
#sa = 1e-2              # selection coefficient of beneficial mutation in 
#UaMax = 1e-5           # max beneficial mutation rate in trait "d"
#Uad = 1e-5             # deleterious mutation rate in trait "d"
#cr  = 0.175            # increment to "c" is (1+cr)
#Ur = 1e-5              # beneficial mutation rate in trait "c"
#Urd = 1e-5             # deleterious mutation rate in trait "c"
#R1 = 1/130.0           # rate of environmental change
#R2 = 1/125.0           # rate of environmental change

parName = {'T':0,'b':1,'do':2,'sa':3,'UaMax':4,'Uad':5,'cr':6,'Ur':7,'Urd':8,'R':9,'part':10}
parValue = np.zeros([11])

with open('evoParameters.csv','r') as csvfile:
    parInput = csv.reader(csvfile)
    for (line,row) in enumerate(parInput):
        parValue[line] = float(row[0])
part = int(parValue[parName['part']])
T = int(parValue[parName['T']])

# calculate list of death terms, i_ext, and upper-bound on death term size
# for population to be viable
[di,i_ext,d_max] = myfun.get_death_terms_classes(parValue[parName['b']], \
                            parValue[parName['do']], parValue[parName['sa']])
di = np.array(list(reversed((di))))

Ua = myfun.get_mutation_rates_RunOutMutModel(parValue[parName['UaMax']],i_ext)
Ua = np.array(list(reversed((Ua))))

positions = (int(math.ceil(len(di)/part)-3))

pfix_a = np.zeros([positions])
with open('parallel_abs.txt','r') as csvfile:
    parInput = csv.reader(csvfile)
    for (line,row) in enumerate(parInput):
        pfix_a[line] = float(row[0])
        
pfix_r = np.zeros([positions])
with open('parallel_abs.txt','r') as csvfile:
    parInput = csv.reader(csvfile)
    for (line,row) in enumerate(parInput):
        pfix_r[line] = float(row[0])
        
pop_eq = np.zeros([positions])
with open('parallel_eqpop.txt','r') as csvfile:
    parInput = csv.reader(csvfile)
    for (line,row) in enumerate(parInput):
        pop_eq[line] = float(row[0])    

'''
parValue = np.zeros([11])
with open('evoParameters.csv','r') as csvfile:
    parInput = csv.reader(csvfile)
    for (line,row) in enumerate(parInput):
        parValue[line] = float(row[0])
'''
yi_option = 3   # option for determining equilibrium density     # 
subsamp = int(len(di)/part)- 1

for k in range(2,subsamp):
    print(di[k*part:(k*part + 2)])


def get_tau_ratio(pfix_a,Ua,pfix_r,Ur):
#calculate ratio of tau's
#inputs
#s,U - selection coeff & mutation rates a-abs fit, r-rel fit
#outputs
# ratio of taus, as defined in Nathan et al. 2021
    tau_a = (1/pfix_a)*np.log(pfix_a/Ua)
    tau_r = (1/pfix_r)*np.log(pfix_r/Ur)
    return tau_a/tau_r

va = np.zeros([positions])
vr = np.zeros([positions])
taus = np.zeros([positions])
taus_old = np.zeros([positions])
cr = parValue[parName['cr']]
sa = parValue[parName['sa']]
Ur = parValue[parName['Ur']]


for k in range(positions):
    va[k] = myfun.get_vDF_pfix_vers(pop_eq[k],sa,Ua[k*part + 2],pfix_a[k])
    eff_cr = myfun.get_c_selection_coefficient(parValue[parName['b']],di[k*part + 2], \
                                                    pop_eq[k]/T,parValue[parName['cr']])
    eff_cr_old = myfun.get_c_selection_coefficient_OLD(parValue[parName['b']],pop_eq[k]/T,cr)
    vr[k] = myfun.get_vDF_pfix_vers(pop_eq[k],eff_cr,Ur,pfix_r[k])
    taus[k] = get_tau_ratio(sa,Ua[k*part + 2],eff_cr,Ur)
    taus_old[k] = get_tau_ratio(sa,Ua[k*part + 2],eff_cr_old,Ur)
absFitClass = range(positions)
fig1,ax1 = plt.subplots(1,1,figsize=[7,7])
ax1.scatter(absFitClass,va,color="blue",linewidth=1.0,label=r'$Eff_Sr$')
ax1.scatter(absFitClass,vr,color="red",linewidth=1.0,label=r'$Eff_Sr_Rescaled$')
#ax1.scatter(absFitClass,taus,color="red",linewidth=1.0,label=r'$Eff_Sr_Rescaled$')
ax1.set_ylim(0,np.max(1.25*vr))

fig1,ax1 = plt.subplots(1,1,figsize=[7,7])
#ax1.scatter(absFitClass,va,color="blue",linewidth=1.0,label=r'$Eff_Sr$')
#ax1.scatter(absFitClass,vr,color="red",linewidth=1.0,label=r'$Eff_Sr_Rescaled$')
ax1.scatter(absFitClass,taus,color="red",linewidth=1.0,label=r'$Eff_Sr_Rescaled$')
#ax1.scatter(absFitClass,taus_old,color="blue",linewidth=1.0,label=r'$Eff_Sr_Rescaled$')
ax1.set_ylim(0,np.max(1.25*taus))