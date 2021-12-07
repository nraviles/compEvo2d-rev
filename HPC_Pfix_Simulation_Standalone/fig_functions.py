# -*- coding: utf-8 -*-
"""
Created on Sat Feb 02 11:25:45 2019
Masel Lab
Project: Mutation-driven Adaptation
@author: Kevin Gomez

Description:
Defines the basics functions used in all scripts that process matlab
data and create figures in the mutation-driven adaptation manuscript.
"""
# libraries
import numpy as np
import scipy.optimize as opt
import bisect
import scipy.stats as st
# my functions in a seperate file
import math as math
import scipy.special

np.seterr(divide='raise')

# *****************************************************************************
# FUNCTIONS TO GET QUANTITIES FROM DESAI AND FISHER 2007
# *****************************************************************************


def get_vDF(N,s,U):
    # Calculates the rate of adaptation v, derived in Desai and Fisher 2007
    # for a given population size (N), selection coefficient (s) and beneficial
    # mutation rate (U)
    #
    # Inputs:
    # N - population size
    # s - selection coefficient
    # U - beneficial mutation rate
    #
    # Output: 
    # v - rate of adaptation
        
    v = s**2*(2*np.log(N*s)-np.log(s/U))/(np.log(s/U)**2)
    
    if np.isnan(v):
        return 0
    else:
        return v

# -----------------------------------------------------------------------------

def get_vDF_lnN_vers(lnN,s,U):
    # Calculates the rate of adaptation v, derived in Desai and Fisher 2007
    # for a given natural log of population size (lnN), selection coefficient (s) and beneficial
    # mutation rate (U)
    #
    # Inputs:
    # lnN - natural log of population size
    # s - selection coefficient
    # U - beneficial mutation rate
    #
    # Output: 
    # v - rate of adaptation
        
    v = s**2*(2*lnN + 2*np.log(s)-np.log(s/U))/(np.log(s/U)**2)
    
    if np.isnan(v):
        return 0
    else:
        return v

# -----------------------------------------------------------------------------

def get_vDF_pfix_vers(N,s,U,pfix):
    # Calculates the rate of adaptation v, derived in Desai and Fisher 2007
    # for a given natural log of population size (lnN), selection coefficient (s) and beneficial
    # mutation rate (U)
    #
    # Inputs:
    # lnN - natural log of population size
    # s - selection coefficient
    # U - beneficial mutation rate
    # pfix - numerically calcualted pfix value for class
    #
    # Output: 
    # v - rate of adaptation
        
    # need to split out into successional and multiple mutations regimes
    
    v = s**2*(2.0*np.log(N*pfix)-np.log(s/U))/(np.log(s/U)**2)
    
    if np.isnan(v):
        return 0
    else:
        return v
#------------------------------------------------------------------------------

def get_tau(s,U):
    # Calculates the rate of adaptation v, derived in Desai and Fisher 2007
    # for a given natural log of population size (lnN), selection coefficient (s) and beneficial
    # mutation rate (U)
    #
    # Inputs:
    # lnN - natural log of population size
    # s - selection coefficient
    # U - beneficial mutation rate
    # pfix - numerically calcualted pfix value for class
    #
    # Output: 
    # v - rate of adaptation
        
    # need to split out into successional and multiple mutations regimes
    
    tau = (1/s)*abs(np.log(s/U))
    
    if np.isnan(tau):
        return 0
    else:
        return tau
#------------------------------------------------------------------------------

def get_qDF(N,s,U):
    # Calculates the rate of adaptation v, derived in Desai and Fisher 2007
    # for a given population size (N), selection coefficient (s) and beneficial
    # mutation rate (U)
    #
    # Inputs:
    # N - population size
    # s - selection coefficient
    # U - beneficial mutation rate
    #
    # Output: 
    # v - rate of adaptation
        
    q = 2*np.log(N*s)/np.log(s/U)
    
    return q

#------------------------------------------------------------------------------
    
def get_death_terms_classes(b,do,sa):
    # function calculates the class for which population size has a negative 
    # growth rate in the Bertram & Masel 2019 lottery model
    #
    # inputs:
    # b - juvenile birth rate
    # do - death term of optimal genotype
    # sa - selection coefficient of beneficial mutations in "d" trait
    #
    # Output: 
    # di - set of death terms for absolute fitness classes
    # d_max - death term at which population is no longer viable
    # i_ext - smallest absolute fitness class with negative growth
    # 
    
    d_max = b + 2 
    di = [do]
    i_ext = 0
    
    while (di[-1] <= d_max):
        di = di+[di[-1]*(1+sa*(di[-1]-1))]
        i_ext = i_ext+1   
            
    di = np.asarray(di)
    return [di,i_ext,d_max]
#------------------------------------------------------------------------------

def get_mutation_rates_RunOutMutModel(UaMax,i_ext):
    # function calculates mutation rates for absolute fitness classes
    #
    # inputs:
    # UaMax - Maximum mutation rate at i_ext
    # i_ext - extinction class (last viable class)
    #
    # Output: 
    # Uai - set of death terms for absolute fitness classes
    #     
    
    Uai = np.zeros([i_ext+1]) 

    for i in range(0,i_ext+1):
        Uai[i] = UaMax*((1.0*i)/i_ext)
        
    return Uai

#------------------------------------------------------------------------------

def get_rate_env_change_scaled(R,di):
    # function calculates rate of environmental change scaled to match rates
    # of adaptation given in generation
    #
    # inputs:
    # R - Rate of environmental change (per unit of time for one iteration)
    # di - death term for class, which determines generation time of population
    #
    # Output: 
    # Ri - Rate of Env. change per generation (class i) 
    #     
    
    Ri = R/(di-1)
    
    return Ri

#------------------------------------------------------------------------------
    
def get_eq_pop_density(b,di,sa,option):
    # Calculate the equilibrium population size for the Bertram & Masel
    # variable density lottery model, single class case.
    #
    # Inputs:
    # b - juvenile birth rate
    # d - death rate
    # i - absolute fitness class (i beneficial mutations from optimal)
    # sa - absolute fitness selection coefficient
    #
    # Output: 
    # eq_density - equilibrium population density
    #
    # Remark: must have 1/sa*d > i >= 0 in general but for i = (1-d/(b+1))/sd
    # the population is in decline, unless it is rescued. Also, returns zero 
    # if formula gives negative density.
    
    def eq_popsize_err(y):    
        # used in numerical approach to obtain equilibrium population density
        return (1-y)*(1-np.exp(-b*y))-(di-1)*y

    if option == 1:
        # approximation near optimal gentotype
        eq_density = (1-np.exp(-b))/(di-np.exp(-b))+(di-1)/(di-np.exp(-b))*(np.exp(-b)-np.exp(-b*(1-np.exp(-b))/(di-np.exp(-b))))/(di-np.exp(-b*(1-np.exp(-b))/(di-np.exp(-b))))
        
        eq_density = np.max([eq_density,0]) # formula should 
        
    elif option == 2:
        # approximation near extinction genotype
        eq_density = (b+2)/(2*b)*(1-np.sqrt(1-8*(b-di+1)/(b+2)**2))
        
        eq_density = np.max([eq_density,0])
        
    else:
        # numerical solution to steady state population size equation
        eq_density = opt.broyden1(eq_popsize_err,[1], f_tol=1e-14)
        
        eq_density = np.max([eq_density,0])
        
    return eq_density

#------------------------------------------------------------------------------
    
def get_c_selection_coefficient(b,di,y,cr):
    # Calculate the "c" selection coefficient for the Bertram & Masel variable 
    # density lottery model for choice of ci = (1+sr)^i
    #
    # Inputs:
    # b - juvenile birth rate
    # di - death term associated with genotype
    # y - equilibrium population density (note that eff_sr has dependence on "d"
    # 	  through the equilibrium density as well).
    # cr - increase to ci from a single beneficial mutation is (1+sr)
    #
    # Output: 
    # eff_sr - selection coefficient of beneficial mutation in "c" trait
    #

    eff_r = cr*(1-y)*(1-(1+b*y)*np.exp(-b*y))/(y+(1-y)*(1-np.exp(-b)))*(b*y-1+np.exp(-b*y))/(b*y*(1-np.exp(-b*y))+cr*(1-(1+b*y)*np.exp(-b*y)))
    
    # scale growth rate to get growth rate per gen (1 gen = 1/(di-1) model iterations)
    eff_sr = eff_r/(di-1)
    
    if np.isnan(eff_sr):
        return 0
    else:
        return eff_sr

def get_c_selection_coefficient_OLD(b,y,cr):
    # Calculate the "c" selection coefficient for the Bertram & Masel variable 
    # density lottery model for choice of ci = (1+sr)^i
    #
    # Inputs:
    # b - juvenile birth rate
    # di - death term associated with genotype
    # y - equilibrium population density (note that eff_sr has dependence on "d"
    # 	  through the equilibrium density as well).
    # cr - increase to ci from a single beneficial mutation is (1+sr)
    #
    # Output: 
    # eff_sr - selection coefficient of beneficial mutation in "c" trait
    #

    eff_r = cr*(1-y)*(1-(1+b*y)*np.exp(-b*y))/(y+(1-y)*(1-np.exp(-b)))*(b*y-1+np.exp(-b*y))/(b*y*(1-np.exp(-b*y))+cr*(1-(1+b*y)*np.exp(-b*y)))
    
    # scale growth rate to get growth rate per gen (1 gen = 1/(di-1) model iterations)
    eff_sr = eff_r
    
    if np.isnan(eff_sr):
        return 0
    else:
        return eff_sr

'''
def deltnplussim(m,c,U): old
#    scatter=np.zeros([int(U),len(m)])
#    for i in range(len(m)):
#        for y in xrange(int(m[i])):
#            cell=int(int(U)*np.random.rand());
#            scatter[cell,i]=scatter[cell,i]+1;
    l=m/float(U)      
    #print(l * U)
    scatter1=np.random.poisson(lam=l[1],size=[U,1]) #dispersion of propugules 
    scatter1 = scatter1[scatter1[:,0]>0];
    scatter2 = np.random.poisson(lam=l[0],size=[len(scatter1),1])
    scatter = np.hstack([scatter2,scatter1])
    Up = len(scatter)
    #print(sum(scatter[:,1]))
    wins = np.zeros(len(m))
#    winsnocomp=np.zeros(len(m)); wins1=np.zeros(len(m)); wins2=np.zeros(len(m));
    comp=np.zeros(Up);
    for i in range(int(Up)):
        comp[i]=sum(scatter[i]) #total number competing per territory
        if comp[i]>0:            
            lotterycmf=np.cumsum(np.array(scatter[i])*c) # Sum mi ci / Sum mi, cbar n *c, n*c + m*c, ..., Sum mi ci is lotterycmf[-1]
            victor=bisect.bisect(lotterycmf,np.random.rand()*lotterycmf[-1]) #random.rand random between 0-1, [0 , c1 m1, c1 m1 + c2 m2] winner based on uniform
            
            wins[victor] = wins[victor] + 1
#            if scatter[i][victor]==1 and comp[i]==1:
#                winsnocomp[victor]=winsnocomp[victor]+1
#            elif scatter[i][victor]==1:
#                wins1[victor]=wins1[victor]+1
#            else: 
#                wins2[victor]=wins2[victor]+1
        # wins based on size: wins based on no comp, based on 1-on-1, 1-on-many
    #print('stochver')
    #print(Up)
    return wins
'''

def deltnplussim(m,c,U):
    l=m/float(U)      

    prob_neq_0 = st.poisson.sf(0, mu=l[1])
    comp_U = np.random.binomial(U,prob_neq_0)
    
    rng = np.arange(1,int(4*math.ceil(l[1]))+1)
    zt_poiss_prbs = (l[1]**rng)/((scipy.special.factorial(rng)*(np.exp(l[1]) - 1)))
    
    comp_mut = np.random.choice(rng,p=zt_poiss_prbs/sum(zt_poiss_prbs),size=[comp_U,1])
    comp_wld =np.random.poisson(lam=l[0],size=[comp_U,1])
    
    scatter = np.hstack([comp_wld,comp_mut])
    Up = len(scatter)
    
    wins = np.zeros(len(m))
    comp=np.zeros(Up);
    for i in range(int(Up)):
        comp[i]=sum(scatter[i]) #total number competing per territory
        if comp[i]>0:            
            lotterycmf=np.cumsum(np.array(scatter[i])*c) # Sum mi ci / Sum mi, cbar n *c, n*c + m*c, ..., Sum mi ci is lotterycmf[-1]
            victor=bisect.bisect(lotterycmf,np.random.rand()*lotterycmf[-1]) #random.rand random between 0-1, [0 , c1 m1, c1 m1 + c2 m2] winner based on uniform
            wins[victor] = wins[victor] + 1
    return wins

def R(m,c,U):
    l=m/float(U)
    L=sum(l)
    cbar=sum(m*c)/sum(m)
    out = l
    for i in range(len(l)):
        try:
            out[i]=cbar*np.exp(-l[i])*(1-np.exp(-(L-l[i])))\
                    /(c[i] + (L-1+np.exp(-L))/(1-(1+L)*np.exp(-L))*(cbar*L-c[i]*l[i])/(L-l[i]))
        except FloatingPointError:
            out[i]=np.exp(-l[i])*(1-np.exp(-(L-l[i])))\
                /(1 + (L-1+np.exp(-L))/(1-(1+L)*np.exp(-L)))
#    for i in range(len(out)):
#        if np.isnan(out)[i]: out[i]=0            
    return out

def A(m,c,U):
    l=m/float(U)
    L=sum(l)
    cbar=sum(m*c)/sum(m)    
    out = l
    
    for i in range(len(l)):    
        try:
            out[i]=cbar*(1-np.exp(-l[i]))\
                    /((1-np.exp(-l[i]))/(1-(1+l[i])*np.exp(-l[i]))*c[i]*l[i]\
                    +(L*(1-np.exp(-L))/(1-(1+L)*np.exp(-L))-l[i]*(1-np.exp(-l[i]))/(1-(1+l[i])*np.exp(-l[i])))/(L-l[i])*(cbar*L-c[i]*l[i]))
        except FloatingPointError:
            out[i]=cbar*(1-np.exp(-l[i]))\
                    /(c[i]+(L*(1-np.exp(-L))/(1-(1+L)*np.exp(-L))-1))
#    for i in range(len(out)):
#        if np.isnan(out)[i]: out[i]=-1
    return out

# This function is drawing divide by zero errors, most likely just numerical precision, but ask kevin
def deltnplus(m,c,U):
    if sum(m)>0 and U>0:
        L=sum(m)/float(U)
        cbar=sum(m*c)/sum(m)
        return m*(np.exp(-L)+(R(m,c,U)+A(m,c,U))*c/cbar)
    else:
        return np.zeros(len(m))
    
def popdeath(pop,di):
    # goal here is the sum of the win vectors
	# then we just add them to the pop
    npop = np.array([np.random.binomial(pop[0],1/di), np.random.binomial(pop[1],1/di)]); #is the number that survive
	# we would then output this number as the new initial population and iterated
	# the process some arbitrary amount of times
    return npop

def di_class(d,di,sa):
    return math.ceil(np.log(float(d) / float(di))/np.log(1+sa))

def simpop(samp,steps,T,sr,b,di,do,sa,de,yi_option): #and all other shit here

	# here we need to take an initial population, with all associated parameters
	# and run through Jasons model i.e.
    c = np.array([1,(1+sr)]);
    pfix = 0;
    print(get_eq_pop_density(b,di,sa,yi_option))
    for i in range(int(samp)):
        eq_pop = math.ceil(T * get_eq_pop_density(b,di,sa,yi_option));
        #print(eq_pop/T)
        pop = np.array([eq_pop-1, 1]);

        while ((pop[1] > 0) & (pop[1] < 1000)): # 4/18: set to 1000 for now, no clue though, != 0 1/pfix calculated via mathematica, grown larger than wild
            U = int(T - sum(pop));
            mi = pop * ((b * U) / T); #or something to generate propugules numbers
            
            deltan = deltnplussim(mi,c,U)
            wins = np.array([0,deltan[1]]) # 4/18: Change to deterministic growth 
            
            npop = pop+wins;
            pop = np.array([(1/di)*(pop[0] + deltnplus(mi,c,U)[0]), np.random.binomial(npop[1],1/di)]); # 4/18: changed to deterministic
            #print(pop[1])
            
        pfix = pfix + (pop[1] > 1)/float(samp);
        
    return pfix
# Everythings running well, only issue is that it spits out nonsensical probabilities sometimes, (how the hell do you get 0.152 from 1 sample?)
#possibility of poisson probabilites not working right
#sum of mi and li.... poiss(mi) =dist= sum_{to U}[poiss(li)]

def trackpop(samp,steps,T,sr,b,di,do,sa,de,yi_option): #and all other shit here

	# here we need to take an initial population, with all associated parameters
	# and run through Jasons model i.e.
    c = np.array([1,(1+sr)]);
    print(get_eq_pop_density(b,di,sa,yi_option))
    tracker = np.array(np.zeros(steps))
    for i in range(int(samp)):
        eq_pop = math.ceil(T * get_eq_pop_density(b,di,sa,yi_option));
        #print(eq_pop/T)
        pop = np.array([eq_pop, 1]);
        trk = np.array(np.zeros(steps))
        trk[0] = 1
        its = 1 
        while ((pop[1] > 0) & (pop[1] < 100) & (its < steps)): # pop2 > 1000, while < 1000 or != 0 1/pfix calculated via mathematica, grown larger than wild
            U = int(T - sum(pop));
            mi = pop * ((b * U) / T); #or something to generate propugules numbers
            #print(mi)
            deltan = deltnplussim(mi,c,U)
            wins = np.array([deltan[0],deltan[1]])
            npop = pop+wins;
            pop = np.array([np.random.binomial(npop[0],1/di), np.random.binomial(npop[1],1/di)]);
            trk[its] = pop[1]
            its = its + 1;
        
            
        tracker = np.vstack((tracker,trk))
        
    return tracker


def compwin(n,samp,T,sr,b,di,do,sa,de,yi_option): #and all other shit here

	# here we need to take an initial population, with all associated parameters
	# and run through Jasons model i.e.
    c = np.array([1,(1+sr)]);
    wins = np.array(np.zeros(samp))
    eq_pop = math.ceil(T * get_eq_pop_density(b,di,sa,yi_option));
        #print(eq_pop/T)
    pop = np.array([eq_pop-n, n]);
    U = int(T - sum(pop));
    mi = pop * ((b * U) / T); #or something to generate propugules numbers
            #print(mi)
    for i in range(samp):
        wins[i] = deltnplussim(mi,c,U)[1]
        
    return wins

#deterministic equations for wild, only consider competition when the mutant has a chance to win, everything else is a lost.  

# i.e. remove the poisson deterministic sampling for the one here: n -> m -> Bin(n + m, 1 - 1/di) = n (wrap)	
# when the approximations break down, assumptions regarding the poisson (when the bins work, binomial -> poisson) - go back through the paper.

def modsimpop(d_Inc,c_Inc,samp,T,sr,b,dis,do,sa,de,yi_option): #and all other shit here

	# here we need to take an initial population, with all associated parameters
	# and run through Jasons model i.e.
    c = np.array([1,(1+sr*c_Inc)]);
    pfix = 0;
    #print(get_eq_pop_density(b,dis[0],sa,yi_option))
    for i in range(int(samp)):
        eq_pop = int(math.ceil(T * get_eq_pop_density(b,dis[0],sa,yi_option)));
        #print(eq_pop/T)
        pop = np.array([eq_pop-1, 1]);

        while ((pop[1] > 0) & (pop[1] < 1000)): # 4/18: set to 1000 for now, no clue though, != 0 1/pfix calculated via mathematica, grown larger than wild
            U = int(T - sum(pop));
            mi = pop * ((b * U) / T); #or something to generate propugules numbers
            deltan = deltnplussim(mi,c,U)

            pop = np.array([(1/dis[0])*(pop[0] + deltnplus(mi,c,U)[0]), np.random.binomial(pop[1]+int(deltan[1]),1/dis[d_Inc])]); # 4/18: changed to deterministic
            # print(pop[1])
            
        pfix = pfix + (pop[1] > 1)/float(samp);
        
    return pfix
#possibility of poisson probabilites not working right
#sum of mi and li.... poiss(mi) =dist= sum_{to U}[poiss(li)]