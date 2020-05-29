#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from cec2013.cec2013 import *
import matplotlib.pyplot as plt

# control parameters
NP=300
FoodNumber=int(NP/2)
limit=100
maxCycle=25
evolve_plot=[0,0] # [evo process plot, final plot]
func_n=4

f = CEC2013(func_n)
dim = f.get_dimension()

global FEs
global MAX_FEs

# Fitness
def calculateFitness(FEs):
    ind_positive=FEs>=0
    ind_negative=FEs<0
    Fitness=(abs(FEs*ind_negative)+1)*ind_negative+1/((FEs*ind_positive)+1)*ind_positive
    return Fitness

# Population init
Foods = np.zeros((FoodNumber,dim))
FEs = np.zeros(FoodNumber)

ub = np.zeros(dim)
lb = np.zeros(dim)
# Get lower, upper bounds
for k in range(dim):
    ub[k] = f.get_ubound(k)
    lb[k] = f.get_lbound(k)

if 1 in evolve_plot:
    x_plot = np.arange(lb, ub, 0.001)
    x_plot.resize(x_plot.size,dim)
    y_plot=np.zeros([x_plot.size,dim])

    for i in range(x_plot.size):
        y_plot[i]=f.evaluate(np.array(x_plot[i]))
        
# Create population within bounds
# init foodsource, FE on Food source
for j in range(FoodNumber):
    Foods[j] = lb + (ub - lb) * np.random.rand(dim)
    FEs[j]=-f.evaluate(Foods[j])

Fitness=calculateFitness(FEs)

trail=np.zeros(FoodNumber)

BestInd=FEs.argmin()
GlobalMin=FEs[BestInd]
GlobalParams=Foods[BestInd]

Cycle=0
while Cycle <=maxCycle:
    for i in range(FoodNumber):
        Param2Change=np.random.randint(dim)
        neigbhbour=np.random.randint(FoodNumber)
        while neigbhbour==i:
            neigbhbour=np.random.randint(FoodNumber)

        sol=Foods[i].copy()
        sol[Param2Change]=Foods[i,Param2Change]+(Foods[i,Param2Change]-Foods[neigbhbour,Param2Change])*(np.random.rand()-0.5)*2
        sol[sol<lb]=lb[sol<lb]
        sol[sol>ub]=ub[sol>ub]

        FE_sol=-f.evaluate(sol)
        Fitness_sol=calculateFitness(FE_sol)

        if Fitness_sol > Fitness[i]:
            Foods[i]=sol
            Fitness[i]=Fitness_sol
            FEs[i]=FE_sol
            trail[i]=0
        else:
            trail[i]=trail[i]+1

        #  Calculate Probabilities          
        prob=(0.9*Fitness/max(Fitness))+0.1

        # onlooker
        i=0
        t=0
        while t<FoodNumber:
            if np.random.rand() < prob[i]:
                t=t+1
                Param2Change=np.random.randint(dim)
                neigbhbour=np.random.randint(FoodNumber)
                while neigbhbour==i:
                    neigbhbour=np.random.randint(FoodNumber)

                sol=Foods[i].copy()
                sol[Param2Change]=Foods[i,Param2Change]+(Foods[i,Param2Change]-Foods[neigbhbour,Param2Change])*(np.random.rand()-0.5)*2
                sol[sol<lb]=lb[sol<lb]
                sol[sol>ub]=ub[sol>ub]

                FE_sol=-f.evaluate(sol)
                Fitness_sol=calculateFitness(FE_sol)

                if Fitness_sol > Fitness[i]:
                    Foods[i]=sol
                    Fitness[i]=Fitness_sol
                    FEs[i]=FE_sol
                    trail[i]=0
                else:
                    trail[i]=trail[i]+1

                i=i+1
                if i==FoodNumber:
                    i=1

        BestInd=FEs.argmin()
        if FEs[BestInd]<GlobalMin:
            GlobalMin=FEs[BestInd]
            GlobalParams=Foods[BestInd]

        # scout bee
        ind=trail.argmax()
        if trail[ind]>limit:
            trail[ind]=0
            sol=(ub-lb)*np.random.rand(dim)+lb
            FE_sol=-f.evaluate(sol)
            Fitness_sol=calculateFitness(FE_sol)
            Foods[ind]=sol
            Fitness[ind]=Fitness_sol
            FEs[ind]=FE_sol
        
    if Cycle%5==6:
        print("Cycle:",Cycle,"FE",GlobalMin)
        X=Foods.copy()
        Y = np.zeros((FoodNumber, dim))
        for i in range(FoodNumber):
            Y[i]=f.evaluate(X[i])
                            
        plt.figure()
        plt.plot(x_plot, y_plot)
        plt.scatter(X,Y,alpha=0.3)
        plt.show()
            
    Cycle=Cycle+1

accuracy = 0.0001
count, seeds = how_many_goptima(Foods, f, accuracy)
N_OPT=f.get_no_goptima()
print("Fun",func_n,"Fuond global optimizers",count,"/",N_OPT)
#print("In the current population there exist",count, "global optimizers.")
# print("Global optimizers:", seeds)