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


# calculate list of death terms, i_ext, and upper-bound on death term size
# for population to be viable
[di,i_ext,d_max] = myfun.get_death_terms_classes(parValue[parName['b']], \
                            parValue[parName['do']], parValue[parName['sa']])
di = np.array(list(reversed((di))))

pfix_a = np.zeros([(int(math.ceil(len(di)/part) - 3))])
with open('parallel_abs.txt','r') as csvfile:
    parInput = csv.reader(csvfile)
    for (line,row) in enumerate(parInput):
        pfix_a[line] = float(row[0])
        
pfix_r = np.zeros([(int(math.ceil(len(di)/part) - 3))])
with open('parallel_abs.txt','r') as csvfile:
    parInput = csv.reader(csvfile)
    for (line,row) in enumerate(parInput):
        pfix_r[line] = float(row[0])
        
pop_eq = np.zeros([(int(math.ceil(len(di)/part)-3))])
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






