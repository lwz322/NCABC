#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from cec2013.cec2013 import *
import matplotlib.pyplot as plt

from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing

import xlwt
# standard DE
def DE(popold, total_FEs, strategy=1, F=0.9, CR=0.1, dim=1, NP=20,func_n=2):
    fun = CEC2013(func_n)
    np.random.seed()

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

        x_r1 = popold[r_set[0]]
        x_r2 = popold[r_set[1]]
        x_r3 = popold[r_set[2]]
        if strategy == 1:
            mutant_vector=x_r1+F*(x_r2-x_r3)

    # Crossover
        crossover_bit_flag = np.random.rand(dim) < CR
        crossover_bit_mask=crossover_bit_flag<0.5
        trail_vector=popold[x_index]*crossover_bit_mask+mutant_vector*crossover_bit_flag

    # bound
        trail_vector=out_bound(trail_vector)
        
    # Selection
        trail_FE=fun.evaluate(trail_vector)
        old_FE=fun.evaluate(popold[x_index])
        total_FEs=total_FEs+2
        if trail_FE >= old_FE:
            popnew[x_index]=trail_vector
        
        #if total_FEs%100000==0:
        #    print(total_FEs,"/",MAX_FEs)
        #if trail_FE < similar_FE:
            #print("evolve")

    return popnew, total_FEs

# standard CDE
def CDE(popold, total_FEs, strategy=1, F=0.5, CR=0.9, dim=1, NP=20,func_n=2):
    fun = CEC2013(func_n)
    np.random.seed()

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
        total_FEs=total_FEs+1
        similar_FE=fun.evaluate(popnew[distance.argmin()])
        total_FEs=total_FEs+1
        if trail_FE >= similar_FE:
            popnew[distance.argmin()]=trail_vector
        
        #if total_FEs%100000==0:
        #    print(total_FEs,"/",MAX_FEs)
        #if trail_FE < similar_FE:
            #print("evolve")

    return popnew, total_FEs

# NCDE
def NCDE(popold, total_FEs, strategy=1, F=0.5, CR=0.9, dim=1, NP=20,func_n=2):
    fun = CEC2013(func_n)
    np.random.seed()

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
        distance=np.linalg.norm(popnew-popnew[x_index],ord=1,axis=1)
        order_distance=np.argsort(distance)
        sorted_pop = popnew[order_distance, :]
        subpop=sorted_pop[0:int(NP/5),:]
        
        r_set = np.random.randint(len(subpop), size=3)
        while len(set(r_set)) < 3 & x_index in r_set:
            r_set = np.random.randint(len(subpop), size=3)

        x_r1 = subpop[r_set[0]]
        x_r2 = subpop[r_set[1]]
        x_r3 = subpop[r_set[2]]
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
        total_FEs=total_FEs+1
        similar_FE=fun.evaluate(popnew[distance.argmin()])
        total_FEs=total_FEs+1
        if trail_FE >= similar_FE:
            popnew[distance.argmin()]=trail_vector
        
        #if total_FEs%100000==0:
        #    print(total_FEs,"/",MAX_FEs)
        #if trail_FE < similar_FE:
            #print("evolve")

    return popnew, total_FEs

def calculateFitness(FEs):
    ind_positive=FEs>=0
    ind_negative=FEs<0
    Fitness=(abs(FEs*ind_negative)+1)*ind_negative+1/((FEs*ind_positive)+1)*ind_positive
    return Fitness

def ABC(Foods,func_n,limit=100):
    FoodNumber=len(Foods)
    f = CEC2013(func_n)
    dim = f.get_dimension()
    MAX_FEs=f.get_maxfes()
    N_OPT=f.get_no_goptima()
    np.random.seed()

    ub = np.zeros(dim)
    lb = np.zeros(dim)
    # Get lower, upper bounds
    for k in range(dim):
        ub[k] = f.get_ubound(k)
        lb[k] = f.get_lbound(k)

    # Fitness
    FEs = np.zeros(FoodNumber)
    for j in range(FoodNumber):
        FEs[j]=-f.evaluate(Foods[j])

    Fitness=calculateFitness(FEs)
    total_FEs=len(Foods)

    trail=np.zeros(FoodNumber)

    BestInd=FEs.argmin()
    GlobalMin=FEs[BestInd]
    GlobalParams=Foods[BestInd]

    result_FEs=0
    result_N_opt=0
    while total_FEs<MAX_FEs:
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
            total_FEs=total_FEs+1
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
                    total_FEs=total_FEs+1
                    Fitness_sol=calculateFitness(FE_sol)

                    if Fitness_sol > Fitness[i]:
                        Foods[i]=sol
                        Fitness[i]=Fitness_sol
                        FEs[i]=FE_sol
                        trail[i]=0
                    else:
                        trail[i]=trail[i]+1

                    i=i+1
                    if i==FoodNumber-1:
                        i=0

            # scout bee
            ind=trail.argmax()
            if trail[ind]>limit:
                trail[ind]=0
                sol=(ub-lb)*np.random.rand(dim)+lb
                FE_sol=-f.evaluate(sol)
                total_FEs=total_FEs+1
                Fitness_sol=calculateFitness(FE_sol)
                Foods[ind]=sol
                Fitness[ind]=Fitness_sol
                FEs[ind]=FE_sol
            
        count, seeds = how_many_goptima(Foods, f, accuracy)
        if count==N_OPT :
            if result_FEs==0:
                result_FEs=total_FEs
                result_N_opt=count
                break

        if total_FEs>=MAX_FEs:
            if result_FEs==0:
                result_FEs=MAX_FEs
                
            result_N_opt=count

    print("Fun",func_n,":","Fuond global optimizers",count,"/",N_OPT)
    return result_N_opt,result_FEs

def MMOPABC(Foods,func_n,limit=100):
    FoodNumber=len(Foods)
    f = CEC2013(func_n)
    dim = f.get_dimension()
    MAX_FEs=f.get_maxfes()
    N_OPT=f.get_no_goptima()
    np.random.seed()

    ub = np.zeros(dim)
    lb = np.zeros(dim)
    # Get lower, upper bounds
    for k in range(dim):
        ub[k] = f.get_ubound(k)
        lb[k] = f.get_lbound(k)

    # Fitness
    FEs = np.zeros(FoodNumber)
    for j in range(FoodNumber):
        FEs[j]=-f.evaluate(Foods[j])

    Fitness=calculateFitness(FEs)
    total_FEs=len(Foods)

    trail=np.zeros(FoodNumber)

    BestInd=FEs.argmin()
    GlobalMin=FEs[BestInd]
    GlobalParams=Foods[BestInd]

    result_FEs=0
    result_N_opt=0
    vice_Foods=Foods.copy()
    vice_Fitness=Fitness.copy() 
    while total_FEs<MAX_FEs:
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
            total_FEs=total_FEs+1
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
                total_FEs=total_FEs+1
                Fitness_sol=calculateFitness(FE_sol)

                if Fitness_sol > Fitness[i]:
                    Foods[i]=sol
                    Fitness[i]=Fitness_sol
                    FEs[i]=FE_sol
                    trail[i]=0
                else:
                    trail[i]=trail[i]+1

                i=i+1
                if i==FoodNumber-1:
                    i=0

        # scout bee
        while trail.max()>limit:
            ind=trail.argmax()
        #if trail[ind]>limit:
            trail[ind]=0
            gf_min_ind=vice_Fitness.argmin()
            vice_Foods[gf_min_ind]=Foods[ind]
            vice_Fitness[gf_min_ind]=Fitness[ind]
            sol=(ub-lb)*np.random.rand(dim)+lb
            FE_sol=-f.evaluate(sol)
            total_FEs=total_FEs+1
            Fitness_sol=calculateFitness(FE_sol)
            Foods[ind]=sol
            Fitness[ind]=Fitness_sol
            FEs[ind]=FE_sol

            
        final_Foods=np.append(Foods,vice_Foods)
        final_Foods.resize(len(Foods),dim)
        count, seeds = how_many_goptima(final_Foods, f, accuracy)
        if count==N_OPT :
            if result_FEs==0:
                result_FEs=total_FEs
                result_N_opt=count
                break

        if total_FEs>=MAX_FEs:
            if result_FEs==0:
                result_FEs=MAX_FEs
                
            result_N_opt=count

    print("Fun",func_n,":","Fuond global optimizers",count,"/",N_OPT)
    return result_N_opt,result_FEs

def CABC(Foods,func_n,limit=100):
    FoodNumber=len(Foods)
    f = CEC2013(func_n)
    dim = f.get_dimension()
    MAX_FEs=f.get_maxfes()
    N_OPT=f.get_no_goptima()

    np.random.seed()
    ub = np.zeros(dim)
    lb = np.zeros(dim)
    # Get lower, upper bounds
    for k in range(dim):
        ub[k] = f.get_ubound(k)
        lb[k] = f.get_lbound(k)
    # Fitness
    FEs = np.zeros(FoodNumber)
    for j in range(FoodNumber):
        FEs[j]=-f.evaluate(Foods[j])

    Fitness=calculateFitness(FEs)
    total_FEs=len(Foods)

    trail=np.zeros(FoodNumber)

    result_FEs=0
    result_N_opt=0
    while total_FEs<MAX_FEs:
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
            total_FEs=total_FEs+1
            Fitness_sol=calculateFitness(FE_sol)

            # Crowd
            distance=np.linalg.norm(Foods-sol,ord=1,axis=1)
            distance[i]=max(distance)
            index_similar=distance.argmin()
            FE_similar=-f.evaluate(Foods[index_similar])
            Fitness_similar=calculateFitness(FE_similar)

            if Fitness_sol > Fitness_similar:
                Foods[index_similar]=sol
                Fitness[index_similar]=Fitness_sol
                FEs[index_similar]=FE_sol
                trail[index_similar]=0
            else:
                trail[index_similar]=trail[index_similar]+1

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
                total_FEs=total_FEs+1
                Fitness_sol=calculateFitness(FE_sol)

                # Crowd
                distance=np.linalg.norm(Foods-sol,ord=1,axis=1)
                distance[i]=max(distance)
                index_similar=distance.argmin()
                Fitness_similar=calculateFitness(FE_similar)

                if Fitness_sol > Fitness_similar:
                    Foods[index_similar]=sol
                    Fitness[index_similar]=Fitness_sol
                    FEs[index_similar]=FE_sol
                    trail[index_similar]=0
                else:
                    trail[index_similar]=trail[index_similar]+1

                if Fitness_sol > Fitness[i]:
                    Foods[i]=sol
                    Fitness[i]=Fitness_sol
                    FEs[i]=FE_sol
                    trail[i]=0
                else:
                    trail[i]=trail[i]+1    

                i=i+1
                if i==FoodNumber-1:
                    i=0

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

        count, seeds = how_many_goptima(Foods, f, accuracy)
        if count==N_OPT :
            if result_FEs==0:
                result_FEs=total_FEs
                result_N_opt=count
                break

        if total_FEs>=MAX_FEs:
            if result_FEs==0:
                result_FEs=MAX_FEs
                
            result_N_opt=count

    print("Fun",func_n,":","Fuond global optimizers",count,"/",N_OPT)
    return result_N_opt,result_FEs

def NCABC(Foods,func_n):
    FoodNumber=len(Foods)
    f = CEC2013(func_n)
    dim = f.get_dimension()
    MAX_FEs=f.get_maxfes()
    N_OPT=f.get_no_goptima()
    accuracy=0.0001

    np.random.seed()
    ub = np.zeros(dim)
    lb = np.zeros(dim)
    # Get lower, upper bounds
    for k in range(dim):
        ub[k] = f.get_ubound(k)
        lb[k] = f.get_lbound(k)
    # Fitness
    FEs = np.zeros(FoodNumber)
    for j in range(FoodNumber):
        FEs[j]=-f.evaluate(Foods[j])

    Fitness=calculateFitness(FEs)
    total_FEs=len(Foods)

    trail=np.zeros(FoodNumber)

    result_FEs=0
    result_N_opt=0
    limit=FoodNumber*dim
        
    vice_Foods=Foods.copy()
    vice_Fitness=Fitness.copy()   
    while total_FEs<MAX_FEs:
        for i in range(FoodNumber):
            Param2Change=np.random.randint(dim)
            
            distance=np.linalg.norm(Foods-Foods[i],ord=1,axis=1)
            order_distance=np.argsort(distance)
            sorted_pop = Foods[order_distance, :]
            subpop=sorted_pop[0:int(FoodNumber/5),:]
            
            neigbhbour=np.random.randint(len(subpop))
            
            while neigbhbour==i:
                neigbhbour=np.random.randint(len(subpop))

            sol=Foods[i].copy()
            sol[Param2Change]=Foods[i,Param2Change]+(Foods[i,Param2Change]-subpop[neigbhbour,Param2Change])*(np.random.rand()-0.5)*2
            sol[sol<lb]=lb[sol<lb]
            sol[sol>ub]=ub[sol>ub]

            FE_sol=-f.evaluate(sol)
            total_FEs=total_FEs+1
            Fitness_sol=calculateFitness(FE_sol)

            # Crowd
            distance=np.linalg.norm(Foods-sol,ord=1,axis=1)
            distance[i]=max(distance)
            index_similar=distance.argmin()
            FE_similar=-f.evaluate(Foods[index_similar])
            Fitness_similar=calculateFitness(FE_similar)

            if Fitness_sol > Fitness_similar:
                Foods[index_similar]=sol
                Fitness[index_similar]=Fitness_sol
                FEs[index_similar]=FE_sol
                trail[index_similar]=0

            elif Fitness_sol > Fitness[i]:
                Foods[i]=sol
                Fitness[i]=Fitness_sol
                FEs[i]=FE_sol
                trail[i]=0
                trail[index_similar]=trail[index_similar]+1
                
            else:
                trail[i]=trail[i]+1
                trail[index_similar]=trail[index_similar]+1

        #  Calculate Probabilities          
        prob=np.zeros(len(Foods))
        for i in range(FoodNumber):
            distance=np.linalg.norm(Foods-Foods[i],ord=1,axis=1)
            order_distance=np.argsort(distance)
            sorted_pop_Fitness = Fitness[order_distance,]
            subpop_Fitness=sorted_pop_Fitness[0:int(FoodNumber/5)]
            prob[i]=(0.9*Fitness[i]/max(subpop_Fitness))+0.1

        # onlooker
        i=0
        t=0
        while t<FoodNumber:
            if np.random.rand() < prob[i]:
                t=t+1
                Param2Change=np.random.randint(dim)
                distance=np.linalg.norm(Foods-Foods[i],ord=1,axis=1)
                order_distance=np.argsort(distance)
                sorted_pop = Foods[order_distance, :]
                subpop=sorted_pop[0:int(FoodNumber/5),:]
                
                neigbhbour=np.random.randint(len(subpop))
                
                while neigbhbour==i:
                    neigbhbour=np.random.randint(len(subpop))
    
                sol=Foods[i].copy()
                sol[Param2Change]=Foods[i,Param2Change]+(Foods[i,Param2Change]-subpop[neigbhbour,Param2Change])*(np.random.rand()-0.5)*2
                sol[sol<lb]=lb[sol<lb]
                sol[sol>ub]=ub[sol>ub]

                FE_sol=-f.evaluate(sol)
                total_FEs=total_FEs+1
                Fitness_sol=calculateFitness(FE_sol)

                # Crowd
                distance=np.linalg.norm(Foods-sol,ord=1,axis=1)
                distance[i]=max(distance)
                index_similar=distance.argmin()
                Fitness_similar=calculateFitness(FE_similar)

                if Fitness_sol > Fitness_similar:
                    Foods[index_similar]=sol
                    Fitness[index_similar]=Fitness_sol
                    FEs[index_similar]=FE_sol
                    trail[index_similar]=0
    
                elif Fitness_sol > Fitness[i]:
                    Foods[i]=sol
                    Fitness[i]=Fitness_sol
                    FEs[i]=FE_sol
                    trail[i]=0
                    trail[index_similar]=trail[index_similar]+1
                    
                else:
                    trail[i]=trail[i]+1
                    trail[index_similar]=trail[index_similar]+1

                i=i+1
                if i==FoodNumber-1:
                    i=0

        # scout bee
        if trail.max()>limit:
            ind=trail.argmax()
            trail[ind]=0
            distance=np.linalg.norm(Foods-Foods[ind],ord=1,axis=1)
            order_distance=np.argsort(distance)
            #guard_bee_order=Fitness[order_distance[0:int(FoodNumber/20)]].argmax()
            global_min_order=Fitness[order_distance[0:int(FoodNumber/5)]].argmin()
            global_max_order=Fitness[order_distance[0:int(FoodNumber/5)]].argmax()
            local_max_order=Fitness[order_distance[0:int(FoodNumber/20)]].argmax()
            global_max_ind=order_distance[global_max_order]
            global_min_ind=order_distance[global_min_order]
            local_max_ind=order_distance[local_max_order]
            if local_max_ind==global_min_order:
                for i in range(1,int(FoodNumber/5)):
                    neighbor_max_order=Fitness[order_distance[0:int(FoodNumber/5)-i]].argmax() 
                    if ind==order_distance[neighbor_max_order]:
                        for reset_ind in order_distance[0:int(FoodNumber/5)-i]:
                            gf_min_ind=vice_Fitness.argmin()
                            vice_Foods[gf_min_ind]=Foods[reset_ind]
                            vice_Fitness[gf_min_ind]=Fitness[reset_ind]
                            trail[reset_ind]=0
                            sol=(ub-lb)*np.random.rand(dim)+lb
                            FE_sol=-f.evaluate(sol)
                            total_FEs=total_FEs+1
                            Fitness_sol=calculateFitness(FE_sol)
                            Foods[reset_ind]=sol
                            Fitness[reset_ind]=Fitness_sol
                            FEs[reset_ind]=FE_sol

            if False:
                trail[ind]=0
                gf_min_ind=vice_Fitness.argmin()
                vice_Foods[gf_min_ind]=Foods[ind]
                vice_Fitness[gf_min_ind]=Fitness[ind]

            if False:
                gf_min_ind=vice_Fitness.argmin()
                vice_Foods[gf_min_ind]=Foods[global_max_ind]
                vice_Fitness[gf_min_ind]=Fitness[global_max_ind]
                trail[global_max_ind]=0
                sol=(ub-lb)*np.random.rand(dim)+lb
                FE_sol=-f.evaluate(sol)
                total_FEs=total_FEs+1
                Fitness_sol=calculateFitness(FE_sol)
                Foods[global_max_ind]=sol
                Fitness[global_max_ind]=Fitness_sol


        final_Foods=np.append(Foods,vice_Foods)
        final_Foods.resize(len(Foods),dim)
        count, seeds = how_many_goptima(final_Foods, f, accuracy)

        if count==N_OPT :
            if result_FEs==0:
                result_FEs=total_FEs
                result_N_opt=count
                break

        if total_FEs>=MAX_FEs:
            if result_FEs==0:
                result_FEs=MAX_FEs
                
            result_N_opt=count
    
    print("Fun",func_n,":","Fuond global optimizers",count,"/",N_OPT,total_FEs)
    return result_N_opt,result_FEs

def par_ABC_run(runtime):
    global func_n
    global NP
    FoodNumber=int(NP/2)

    f = CEC2013(func_n)
    dim = f.get_dimension()
    limit=FoodNumber*dim 
    np.random.seed()

    # Population init
    Foods = np.zeros((FoodNumber,dim))

    ub = np.zeros(dim)
    lb = np.zeros(dim)
    # Get lower, upper bounds
    for k in range(dim):
        ub[k] = f.get_ubound(k)
        lb[k] = f.get_lbound(k)

    # Create population within bounds
    # init foodsource, FE on Food source
    for j in range(FoodNumber):
        Foods[j] = lb + (ub - lb) * np.random.rand(dim)

    result_N_opt,result_FEs=MMOPABC(Foods,func_n)
    return result_N_opt,result_FEs

def par_CDE_run(runtime):
    result_N_opt=0
    result_FEs=0
    global func_n
    global NP
    # init
    np.random.seed()
    f = CEC2013(func_n)
    dim = f.get_dimension()
    MAX_FEs=f.get_maxfes()
    N_OPT=f.get_no_goptima()
    
    # pop init
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
    
    total_FEs=0
    while total_FEs<MAX_FEs:
        pop_new, total_FEs = NCDE(X, total_FEs, 1, 0.5, 0.9, dim, NP, func_n)
        X=pop_new.copy()
        count, seeds = how_many_goptima(X, f, accuracy)
        if count==N_OPT :
            if result_FEs==0:
                result_FEs=total_FEs
                result_N_opt=count
                break

        if total_FEs>=MAX_FEs:
            if result_FEs==0:
                result_FEs=MAX_FEs
                
            result_N_opt=count
            
    print("Fun",func_n,":",runtime+1,"RUN:","Fuond global optimizers",count,"/",N_OPT)
    return result_N_opt,result_FEs        

def CEC_benchmark(benchmark_fun_set=[1]):
    global func_n
    global runtimes
    global accuracy
    global MAX_FEs
    global NP
    accuracy = 0.0001
    runtimes=50
    
    book = xlwt.Workbook(encoding='utf-8',style_compression=0)
    sheet = book.add_sheet('CDE',cell_overwrite_ok=True)
    result_N_opt=np.zeros(runtimes)
    result_FEs=np.zeros(runtimes)
    for func_n in benchmark_fun_set:
        f = CEC2013(func_n)
        N_OPT=f.get_no_goptima()
        MAX_FEs=f.get_maxfes()
        #NP=int(0.5*f.get_popsize())
        NP=f.get_popsize()
        #NP=120
        par_input=range(runtimes)

        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=10)
        #pool = ThreadPool()
        index=0
        for A,B in pool.map(par_ABC_run,par_input):
            result_N_opt[index]=A
            result_FEs[index]=B
            index=index+1

        pool.close()
        pool.join()    
            
        #print(result_N_opt)
        result_PR=sum(result_N_opt)/(N_OPT*runtimes)
        result_SR=sum(result_FEs<MAX_FEs)/runtimes
        AveFEs=sum(result_FEs)/runtimes
        print(result_PR,result_SR,AveFEs)   

        sheet.write(func_n-1,0,result_PR)
        sheet.write(func_n-1,1,result_SR)
        sheet.write(func_n-1,2,AveFEs)
        sheet.write(func_n-1,3,np.std(result_FEs))
        
    book.save('./MMOPCABC_all_NP_final.xls')
    
CEC_benchmark(benchmark_fun_set=range(1,21))