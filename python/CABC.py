#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from cec2013.cec2013 import *
import matplotlib.pyplot as plt

def calculateFitness(FEs):
    ind_positive=FEs>=0
    ind_negative=FEs<0
    Fitness=(abs(FEs*ind_negative)+1)*ind_negative+1/((FEs*ind_positive)+1)*ind_positive
    return Fitness

def CABC(Foods,func_n):
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
        
    guard_Foods=Foods.copy()
    guard_Fitness=Fitness.copy()
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
            
            ind=trail.argmax()
            if trail[ind]>limit:
                trail[ind]=0
                distance=np.linalg.norm(Foods-Foods[ind],ord=1,axis=1)
                order_distance=np.argsort(distance)
                #guard_bee_order=Fitness[order_distance[0:int(FoodNumber/20)]].argmax()
                guard_bee_order=Fitness[order_distance[0:int(FoodNumber/10)]].argmin()
                guard_bee_ind=order_distance[guard_bee_order]
                if guard_bee_ind==ind:
                    trail[guard_bee_ind]=0
                    sol=(ub-lb)*np.random.rand(dim)+lb
                    FE_sol=-f.evaluate(sol)
                    total_FEs=total_FEs+1
                    Fitness_sol=calculateFitness(FE_sol)
                    Foods[guard_bee_ind]=sol
                    Fitness[guard_bee_ind]=Fitness_sol
                    FEs[guard_bee_ind]=FE_sol
                #gf_min_ind=guard_Fitness.argmin()
                #guard_Foods[gf_min_ind]=guard_bee
                #guard_Fitness[gf_min_ind]=Fitness[guard_bee_ind]
                
                #for reset_ind in order_distance[0:int(FoodNumber/20)]:
                    
                #    if reset_ind==guard_bee_ind:
                #        trail[reset_ind]=limit-5
                #        continue
                    
                #    trail[reset_ind]=0
                #    sol=(ub-lb)*np.random.rand(dim)+lb
                #    FE_sol=-f.evaluate(sol)
                #    total_FEs=total_FEs+1
                #    Fitness_sol=calculateFitness(FE_sol)
                #    Foods[reset_ind]=sol
                #    Fitness[reset_ind]=Fitness_sol
                #    FEs[reset_ind]=FE_sol
                    
                    
        final_Foods=np.append(Foods,guard_Foods)
        final_Foods.resize(len(Foods),dim)
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
    
    
    pop_plot(Foods,func_n)    
    print("Fun",func_n,":","Fuond global optimizers",count,"/",N_OPT,total_FEs)
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
    limit=FoodNumber
        
    vice_Foods=Foods.copy()
    flag=1
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
            sorted_pop_Fitness = Fitness[order_distance]
            subpop_Fitness=sorted_pop_Fitness[0:int(FoodNumber/10)]
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
            gf_min_ind=vice_Fitness.argmin()
            vice_Foods[gf_min_ind]=Foods[ind]
            vice_Fitness[gf_min_ind]=Fitness[ind]
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
                trail[ind]=0

        final_Foods=np.append(Foods,vice_Foods)
        final_Foods.resize(len(Foods),dim)
        count, seeds = how_many_goptima(final_Foods, f, accuracy)
        
        
        if total_FEs/40000>flag:
            flag=flag+1
            pop_plot(final_Foods,func_n) 
            
        if count==N_OPT :
            if result_FEs==0:
                result_FEs=total_FEs
                result_N_opt=count
                break

        if total_FEs>=MAX_FEs:
            if result_FEs==0:
                result_FEs=MAX_FEs
                
            result_N_opt=count
    
    pop_plot(final_Foods,func_n)  
    print("Fun",func_n,":","Fuond global optimizers",count,"/",N_OPT,total_FEs)
    return result_N_opt,result_FEs

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
            
        count, seeds = how_many_goptima(Foods, f, 0.0001)
        if count==N_OPT :
            if result_FEs==0:
                result_FEs=total_FEs
                result_N_opt=count
                break

        if total_FEs>=MAX_FEs:
            if result_FEs==0:
                result_FEs=MAX_FEs
                
            result_N_opt=count

    pop_plot(Foods,func_n)    
    print("Fun",func_n,":","Fuond global optimizers",count,"/",N_OPT)
    return result_N_opt,result_FEs

def par_ABC_run(runtime):
    global func_n
    global NP
    FoodNumber=int(NP/2)

    f = CEC2013(func_n)
    dim = f.get_dimension()

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

    result_N_opt,result_FEs=CABC(Foods,func_n)
    return result_N_opt,result_FEs

def pop_plot(X,func_n):
    f = CEC2013(func_n)
    dim = f.get_dimension()
    ub = np.zeros(dim)
    lb = np.zeros(dim)
    # Get lower, upper bounds
    for k in range(dim):
        ub[k] = f.get_ubound(k)
        lb[k] = f.get_lbound(k)


    if dim==1:
        x_plot = np.arange(lb, ub, 0.001)
        x_plot.resize(x_plot.size,dim)
        y_plot=np.zeros([x_plot.size,dim])
    
        for i in range(x_plot.size):
            y_plot[i]=f.evaluate(np.array(x_plot[i]))
        
        Y = np.zeros((len(X), dim))
        for i in range(len(X)):
            Y[i]=f.evaluate(X[i])
                            
        plt.figure()
        plt.plot(x_plot, y_plot)
        plt.xlabel("X")
        plt.ylabel("Fun_eval")
        plt.scatter(X,Y,alpha=0.1)
        plt.show()   
     
        
    if dim==2:
        x_plot = np.arange(lb[0], ub[0], 0.1)
        y_plot = np.arange(lb[1], ub[1], 0.1)
        x_plot.resize(x_plot.size,1)
        y_plot.resize(y_plot.size,1)
        A, B = np.meshgrid(x_plot, y_plot)
        plt.figure()
        #C, D = np.meshgrid(Foods[:,0],Foods[:,1])
        #grid=np.concatenate((x_plot,y_plot),axis=1)
        plt.xlim(lb[0],ub[0])
        plt.ylim(lb[1],ub[1])
        plt.scatter(X[:,0],X[:,1],alpha=0.1)
        
        if False:
            from mpl_toolkits.mplot3d import Axes3D
            Y = np.zeros((len(x_plot),len(x_plot)))
            for i in range(len(x_plot)):
                for j in range(len(y_plot)):
                    Y[i,j]=f.evaluate(np.append(x_plot[i],y_plot[j]))
    
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.plot_surface(A, B, Y, rstride=1, cstride=1, cmap='rainbow')
            plt.show()
            
    if dim==3:
        from mpl_toolkits.mplot3d import Axes3D
        x = X[:, 0]  # [ 0  3  6  9 12 15 18 21]
        y = X[:, 1]  # [ 1  4  7 10 13 16 19 22]
        z = X[:, 2]  # [ 2  5  8 11 14 17 20 23]
         
         
        # 绘制散点图
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(x, y, z)
        plt.show()
            
        
func_n=13
for func_n in range(12,13):
    #NP=400    
    f = CEC2013(func_n)
    NP=2*f.get_popsize()
    FoodNumber=int(NP/2)
    dim = f.get_dimension()
    
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
    
    result_N_opt,result_FEs=NCABC(Foods,func_n)