# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 08:46:05 2021

@author: wyseq

This is a work of Icelandic fiction, be advised.

Code for modified adaptation rates graph generation.
"""
######Libre#######
import numpy as np
import matplotlib.pyplot as plt
import bisect
import scipy as sp
import matplotlib.pyplot as plt
# my functions in a seperate file
import scipy.optimize as opt
from scipy.optimize import fsolve
from scipy.optimize import root 
import fig_functions as myfun  
import pfix_sim_functions as pfix_sim
######Libre#######
######Datum#######
T = 1e9         # carrying capacity

b=2.0           # juvenile birth rate
do=100/98.0      # minimum death rate / death rate of optimal genotype
di=np.array([1.54, 1.576, 1.612, 1.648, 1.684, 1.72, 1.756, 1.792,  1.684, 1.684, 1.684, 1.684, 1.684, 1.684, 1.684, 1.684]);       # juvenile birth rate
#b = np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.82, 1.748, 1.676, 1.604, 1.532, 1.46, 1.388, 1.316]);

sa = 1e-2       # selection coefficient of beneficial mutation in 
sr = 0.175     # multiplicative increment to "c" is (1+sr)

yi_option = 3   # select option to get equilibrium population density (=1,2,3)

max_index = 400
de = b+1


sa = 1e-2       # selection coefficient of beneficial mutation in 
sr = 0.175     # multiplicative increment to "c" is (1+sr)

yi_option = 3   # select option to get equilibrium population density (=1,2,3)
######Datum#######
# we will be using modsimpop(d_Inc,c_Inc,samp,T,sr,b,dis,do,sa,de,yi_option)
# here 'dis' has been modified to take in 2 death rates, and we have removed the step size restriction
#pfixes = np.zeros(len(nsample))
#d_Inc, and c_Inc consider the inclusion of unuequal competition or death in the model
# 0 = "exclude", 1 = "include" 

# we need to set up the 'dis' to take out the first and next term so 
# dis[j:j+2]
for i in range(len(di)-1):
    print(di[i:(i+2)])
# as an examples
    
# there is no need to modify the 'dis', 'sa' or 'sr' terms as it is controlled by
# the inclusion terms
    
# in actuality our graphs go from di = 102/98 to 3 or so, but for now lets ignore that fact

# Vectors for pfixes in each type
d_pfixes = np.zeros(len(di)-1)
c_pfixes = np.zeros(len(di)-1)

# Now the sample number is questionable; how long are we willing and able to let it run?

samp = 1

# Pure Death
d_Inc = 1
c_Inc = 0
for i in range(len(di)-1):
    d_pfixes[i] = pfix_sim.modsimpop(d_Inc,c_Inc,samp,T,sr,b,di[i:(i+2)],do,sa,de,yi_option) #fix
    #array([0.024, 0.024, 0.008])

# Pure Competition
d_Inc = 0
c_Inc = 1
for i in range(len(di)-1):
    c_pfixes[i] = pfix_sim.modsimpop(d_Inc,c_Inc,samp,T,sr,b,di[i:(i+2)],do,sa,de,yi_option) #fix

