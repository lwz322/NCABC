#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 03:22:51 2020

@author: lwz322
"""

import numpy as np
from cec2013.cec2013 import *
import matplotlib.pyplot as plt

# standard DE
def DE(popold, strategy=1, F=0.7, CR=0.2, dim=1, NP=20,func_n=2):
    fun = CEC2013(func_n)
    global FEs
    global MAX_FEs
    def out_bound(trail_vector):
        if np.random.rand(1) > 0.5:
            trail_vector=(trail_vector < lb) * lb + (trail_vector >= lb) * trail_vector
            trail_vector=(trail_vector > ub) * ub + (trail_vector <= ub) * trail_vector
        else:
            trail_vector=(trail_vector < lb) * (lb + np.random.rand(dim) * (ub-lb)) + (trail_vector >= lb) * trail_vector
            trail_vector=(trail_vector > ub) * (lb + np.random.rand(dim) * (ub-lb)) + (trail_vector <= ub) * trail_vector            
        
        return trail_vector
    
    ub = np.zeros(dim)
    lb = np.zeros(dim)
    # Get lower, upper bounds
    for k in range(dim):
        ub[k] = fun.get_ubound(k)
        lb[k] = fun.get_lbound(k)

    popnew=popold.copy()
    for x_index in range(NP):
        r_set = np.random.randint(NP, size=5)
        while len(set(r_set)) < 5 & x_index in r_set:
            r_set = np.random.randint(NP, size=5)

        x_r1 = popnew[r_set[0]]
        x_r2 = popnew[r_set[1]]
        x_r3 = popnew[r_set[2]]
        if strategy == 1:
            mutant_vector=x_r1+F*(x_r2-x_r3)

# Crossover
        crossover_bit_flag = np.random.rand(dim) < CR
        crossover_bit_mask=crossover_bit_flag<0.5
        trail_vector=popnew[x_index]*crossover_bit_mask+mutant_vector*crossover_bit_flag

# bound
        trail_vector=out_bound(trail_vector)
        
# Crowd
        distance=np.linalg.norm(popnew-trail_vector,ord=1,axis=1)
        distance[x_index]=max(distance)
        
# Selection
        trail_FE=fun.evaluate(trail_vector)
        FEs=FEs+1
        similar_FE=fun.evaluate(popnew[distance.argmin()])
        FEs=FEs+1
        if trail_FE >= similar_FE:
            popnew[distance.argmin()]=trail_vector
        
        if FEs%100000==0:
            print(FEs,"/",MAX_FEs)
        #if trail_FE < similar_FE:
            #print("evolve")

    return popnew, FEs

# Population init
evolve_plot=[0,0] # [evo process plot, final plot]
func_n=6
global FEs 
global MAX_FEs

f = CEC2013(func_n)
dim = f.get_dimension()
NP=300

X = np.zeros((NP, dim))
ub = np.zeros(dim)
lb = np.zeros(dim)
# Get lower, upper bounds
for k in range(dim):
    ub[k] = f.get_ubound(k)
    lb[k] = f.get_lbound(k)

# Create population within bounds
for j in range(NP):
    X[j] = lb + (ub - lb) * np.random.rand(1, dim)

if 1 in evolve_plot:
    x_plot = np.arange(lb, ub, 0.001)
    x_plot.resize(x_plot.size,dim)
    y_plot=np.zeros([x_plot.size,dim])

    for i in range(x_plot.size):
        y_plot[i]=f.evaluate(np.array(x_plot[i]))
    
# evolve
#MAX_Gen=300
MAX_FEs=f.get_maxfes()
N_OPT=f.get_no_goptima()
FEs=0
#for i in range(MAX_Gen):
while FEs<MAX_FEs:
    pop_new=DE(X, 1, 0.7, 0.2, dim, NP,func_n)
    #if sum(pop_new-X)!=0:
    if evolve_plot[0]==True:
        Y = np.zeros((NP, dim))
        for i in range(NP):
            Y[i]=f.evaluate(X[i])
        
        new_Y = np.zeros((NP, dim))
        for i in range(NP):
            new_Y[i]=f.evaluate(pop_new[i])
            
        plt.figure()
        plt.plot(x_plot, y_plot)
        plt.scatter(X,Y,alpha=0.6)
        plt.plot(pop_new,new_Y,'rx')
        plt.show()
        
    X=pop_new.copy()

# Calculate how many global optima are in the population
accuracy = 0.001
count, seeds = how_many_goptima(X, f, accuracy)
print("Fun",func_n,"Fuond global optimizers",count,"/",N_OPT)
#print("In the current population there exist",
      count, "global optimizers.")
# print("Global optimizers:", seeds)

if evolve_plot[1]==True:
    Y = np.zeros((NP, dim))
    for i in range(NP):
        Y[i]=f.evaluate(X[i])
    
    new_Y = np.zeros((NP, dim))
    for i in range(NP):
        new_Y[i]=f.evaluate(pop_new[i])
        
    plt.figure()
    plt.plot(x_plot, y_plot)
    plt.scatter(X,Y,alpha=0.6)
    plt.plot(pop_new,new_Y,'rx')
    plt.show()