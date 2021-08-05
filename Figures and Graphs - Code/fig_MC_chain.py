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

# import pfix values derived from simulations/interpolation from csv file (parallel_abs.csv),
# pfix values should be listed from first class 0 (at optimal death term) to last class (extinction).
pfixValsAbs = np.zeros(di.shape)
with open('parallel_abs.csv','r') as csvfile:
    pfixInput = csv.reader(csvfile)
    for (line,row) in enumerate(pfixInput):
        pfixValsAbs[line] = float(row[0])

pfixValsRel = np.zeros(di.shape)
with open('parallel_rel.csv','r') as csvfile:
    pfixInput = csv.reader(csvfile)
    for (line,row) in enumerate(pfixInput):
        pfixValsRel[line] = float(row[0])

# set root solving option for equilibrium densities
yi_option = 3   # select option to get equilibrium population density (=1,2,3)
subsamp = len(di)/part
# set up names to evolution parameters at each class.
eq_y_i =  np.zeros(subsamp)     # equilibrium density of fitness class i
eq_N_i =  np.zeros(subsamp)    # equilibrium population size of fitness class i
eff_sr_i =  np.zeros(subsamp)  # selection coefficient of "c" trait beneficial mutation
eff_sa_i =  np.zeros(subsamp)  # selection coefficient of "c" trait beneficial mutation

va_i =  np.zeros(subsamp)      # rate of adaptation in absolute fitness trait alone
vr_i =  np.zeros(subsamp)      # rate of adaptation in relative fitness trait alone
ve_i =  np.zeros(subsamp)      # rate of fitness decrease due to environmental degradation

prob_incr_abs =  np.zeros(subsamp)   # transition probability i->i-1
prob_decr_abs =  np.zeros(subsamp)   # transition probability i->i+1

# calculate equilibrium population size,density, and rates of adaptation
for i in range(subsamp):
    # absolute fitness mutation rate, equilb.-density,equilb.-popsize,eff_sr
    eq_y_i[i] = myfun.get_eq_pop_density(parValue[parName['b']],di[i*part], \
                                        parValue[parName['sa']],yi_option)
    
    eq_N_i[i] = parValue[parName['T']]*eq_y_i[i]
    eff_sr_i[i] = myfun.get_c_selection_coefficient(parValue[parName['b']],di[i*part], \
                                                    eq_y_i[i],parValue[parName['cr']])
    
    # rates of fitness change (relative scale)
    va_i[i] = myfun.get_vDF_pfix_vers(eq_N_i[i],parValue[parName['sa']],Uai[i*part],pfixValsAbs[i])
    vr_i[i] = myfun.get_vDF_pfix_vers(eq_N_i[i],eff_sr_i[i],parValue[parName['Ur']],pfixValsRel[i])
    ve_i[i] = parValue[parName['sa']]*myfun.get_rate_env_change_scaled(parValue[parName['R']],di[i*part])
    
    prob_incr_abs[i] = va_i[i]/parValue[parName['sa']]
    prob_decr_abs[i] = myfun.get_rate_env_change_scaled(parValue[parName['R']],di[i*part])
    
#coeff_dT_a = (eff_sa_i/(np.log(eff_sa_i/Ua_i)))**2
#coeff_dT_r = (eff_sr_i/(np.log(eff_sr_i/Ur)))**2

# ------------------ figure of Markov chain approximation --------------------

fig1,ax1 = plt.subplots(1,1,figsize=[7,7])
ax1.scatter(absFitClass,va_i,color="blue",linewidth=1.0,label=r'$v_a$')
ax1.scatter(absFitClass,vr_i,color="red",linewidth=1.0,label=r'$v_r$')
ax1.plot(absFitClass,ve_i,color="green",linewidth=2,label=r' $R$')
#ax1.scatter([-1.89],[sa/103.0],color="none",s=100,edgecolor="black",linewidth=2.0,marker='o')
#ax1.scatter([-1.59],[sa/115],color="none",s=100,edgecolor="black",linewidth=2.0,marker='o')

#ax1.plot(d_i,ve_i,color="green",linewidth=2)
#ax1.plot(d_i,(130.0/103.0)*ve_i,color="black",linewidth=2)

ax1.set_xlim(-i_ext-1,0)
ax1.set_ylim(0,np.max(1.25*vr_i))
ax1.set_xticks([-3.0+.5*i for i in range(0,5)])
ax1.set_xticklabels([3.0-.5*i for i in range(0,5)],fontsize=18)

ax1.set_yticks([10**(-4)*(0+0.2*i) for i in range(0,7)])
ax1.set_yticklabels([str(0+0.2*i) for i in range(0,7)],fontsize=18)

ax1.set_xlabel(r'Death term (d)',fontsize=22,labelpad=8)
ax1.set_ylabel(r'Rate of adaptation',fontsize=22,labelpad=8)
#ax1.ticklabel_format(style='sci',axis='y',scilimits=None)

handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,0,1]

ax1.legend([handles[idx] for idx in order],[labels[idx] for idx in order],fontsize = 20,ncol=2)
#plt.text(-20,1.2e-5 ,r'$i_{max} = -98$',fontsize=16)
#plt.text(-20,0.5e-5 ,r'$i_{ext} \approx -64$',fontsize=16)

#ax1.plot([-1.89,-1.89],[0,sa/104],c="black",linewidth=2,linestyle='--')
#ax1.annotate("", xy=(-1.92,sa/106), xytext=(-2.1, sa/106),arrowprops={'arrowstyle':'-|>','lw':3})
#ax1.annotate("", xy=(-1.86,sa/100), xytext=(-1.69, sa/100),arrowprops={'arrowstyle':'-|>','lw':3})
#
#ax1.plot([-1.59,-1.59],[0,sa/116],c="black",linewidth=2,linestyle='--')
#ax1.annotate("", xy=(-1.62,sa/115), xytext=(-1.79, sa/115),arrowprops={'arrowstyle':'-|>','lw':3})
#ax1.annotate("", xy=(-1.55,sa/115), xytext=(-1.37, sa/115),arrowprops={'arrowstyle':'-|>','lw':3})
sa = parValue[parName['sa']]
ax1.plot([-1.89,-1.89],[sa/1000,sa/104],c="black",linewidth=2,linestyle='--')
ax1.plot([-1.59,-1.59],[sa/300,sa/116],c="black",linewidth=2,linestyle='--')

ax1.annotate("", xy=(-1.92,sa/1000), xytext=(-2.1, sa/1000),arrowprops={'arrowstyle':'-|>','lw':3,'color':'blue'})
ax1.annotate("", xy=(-1.86,sa/1000), xytext=(-1.69, sa/1000),arrowprops={'arrowstyle':'-|>','lw':3})

ax1.annotate("", xy=(-1.62,sa/300), xytext=(-1.79, sa/300),arrowprops={'arrowstyle':'-|>','lw':3,'color':'blue'})
ax1.annotate("", xy=(-1.55,sa/300), xytext=(-1.37, sa/300),arrowprops={'arrowstyle':'-|>','lw':3,'color':'red'})

plt.text(-3.35,1.32*10**(-4),r'$\times 10^{-4}$', fontsize = 20)

# add extra legend to point out that solid line is for absolute fitness 
# rate of adaptation and dashed is for relative fitness rate of adaptation
#ax2 = ax1.twinx()
#ax2.set_xlim(-3.1,-0.8)
#ax2.set_ylim(0,np.max(1.25*Va))
#ax2.set_yticklabels([])
#
#plot_lines = []
#plt.hold(True)
#c = v_colors[0]
#l1, = plt.plot([0],[1], '-', color="Black")
#l2, = plt.plot([0],[2], '-', color="Green")
#
#plot_lines.append([l1, l2,])
#
#legend1 = plt.legend(plot_lines[0], [r'$R_{130}$', r'$R_{103}$'], loc=2,fontsize=16)
##plt.legend([l[0] for l in plot_lines], [1 2], loc=4)
#plt.gca().add_artist(legend1)

plt.tight_layout()
fig1.savefig('fig_AbsVsRel_MC_chain.pdf')

