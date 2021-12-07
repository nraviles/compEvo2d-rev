# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 08:36:01 2021

@author: wyseq
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 11:56:52 2020
@author: Kevin Gomez, Nathan Aviles
Masel Lab
see Bertram, Gomez, Masel 2016 for details of Markov chain approximation
see Bertram & Masel 2019 for details of lottery model
"""

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import fig_functions as myfun            # my functions in a seperate file
import csv

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
# calculate list of death terms, i_ext, and upper-bound on death term size
# for population to be viable
[di,i_ext,d_max] = myfun.get_death_terms_classes(parValue[parName['b']], \
                            parValue[parName['do']], parValue[parName['sa']])

# calculate mutation rates (running out of mutations model)
Uai = myfun.get_mutation_rates_RunOutMutModel(parValue[parName['UaMax']],i_ext)

# TEMP: plot di with classes to check values.
absFitClass = -np.asarray(range(0,i_ext+1))

# set root solving option for equilibrium densities
yi_option = 3   # select option to get equilibrium population density (=1,2,3)
subsamp = len(di)#/part
# set up names to evolution parameters at each class.
eq_y_i =  np.zeros(subsamp)     # equilibrium density of fitness class i
eq_N_i =  np.zeros(subsamp)    # equilibrium population size of fitness class i
eff_sr_i =  np.zeros(subsamp)  # selection coefficient of "c" trait beneficial mutation
eff_sr_k =  np.zeros(subsamp)
eff_sa_i =  np.zeros(subsamp)  # selection coefficient of "c" trait beneficial mutation

va_i =  np.zeros(subsamp)      # rate of adaptation in absolute fitness trait alone
vr_i =  np.zeros(subsamp)      # rate of adaptation in relative fitness trait alone
ve_i =  np.zeros(subsamp)      # rate of fitness decrease due to environmental degradation
ve_i_us =  np.zeros(subsamp) 
prob_incr_abs =  np.zeros(subsamp)   # transition probability i->i-1
prob_decr_abs =  np.zeros(subsamp)   # transition probability i->i+1

# calculate equilibrium population size,density, and rates of adaptation
for i in range(len(di)):
    # absolute fitness mutation rate, equilb.-density,equilb.-popsize,eff_sr
    eq_y_i[i] = myfun.get_eq_pop_density(parValue[parName['b']],di[i], \
                                        parValue[parName['sa']],yi_option)
    
    eq_N_i[i] = parValue[parName['T']]*eq_y_i[i]
    eff_sr_i[i] = myfun.get_c_selection_coefficient(parValue[parName['b']],di[i], \
                                                    eq_y_i[i],parValue[parName['cr']])
    eff_sr_k[i] = myfun.get_c_selection_coefficient_OLD(parValue[parName['b']], \
                                                    eq_y_i[i],parValue[parName['cr']])
    
    ve_i[i] = parValue[parName['sa']]*myfun.get_rate_env_change_scaled(parValue[parName['R']],di[i])
    ve_i_us[i] = parValue[parName['sa']] * parValue[parName['R']]
#coeff_dT_a = (eff_sa_i/(np.log(eff_sa_i/Ua_i)))**2
#coeff_dT_r = (eff_sr_i/(np.log(eff_sr_i/Ur)))**2

# ------------------ figure of Markov chain approximation --------------------

fig1,ax1 = plt.subplots(1,1,figsize=[7,7])
ax1.scatter(absFitClass,eff_sr_k,color="blue",linewidth=1.0,label=r'$Eff_Sr$')
ax1.scatter(absFitClass,eff_sr_i,color="red",linewidth=1.0,label=r'$Eff_Sr_Rescaled$')
#ax1.scatter([-1.89],[sa/103.0],color="none",s=100,edgecolor="black",linewidth=2.0,marker='o')
#ax1.scatter([-1.59],[sa/115],color="none",s=100,edgecolor="black",linewidth=2.0,marker='o')

#ax1.plot(d_i,ve_i,color="green",linewidth=2)
#ax1.plot(d_i,(130.0/103.0)*ve_i,color="black",linewidth=2)

ax1.set_xlim(-i_ext-1,0)
ax1.set_ylim(0,np.max(1.25*max(eff_sr_i)))

plt.tight_layout()
fig1.savefig('fig_AbsVsRel_MC_chain.pdf')

fig1,ax1 = plt.subplots(1,1,figsize=[7,7])
ax1.scatter(absFitClass,ve_i_us,color="blue",linewidth=1.0,label=r'$Eff_Sr$')
ax1.scatter(absFitClass,ve_i,color="red",linewidth=1.0,label=r'$Eff_Sr_Rescaled$')
#ax1.scatter([-1.89],[sa/103.0],color="none",s=100,edgecolor="black",linewidth=2.0,marker='o')
#ax1.scatter([-1.59],[sa/115],color="none",s=100,edgecolor="black",linewidth=2.0,marker='o')

#ax1.plot(d_i,ve_i,color="green",linewidth=2)
#ax1.plot(d_i,(130.0/103.0)*ve_i,color="black",linewidth=2)

ax1.set_xlim(-i_ext-1,0)
ax1.set_ylim(0,np.max(1.25*max(ve_i)))

#plt.tight_layout()
#fig1.savefig('fig_AbsVsRel_MC_chain.pdf')