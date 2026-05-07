
import random
import math
import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.optimize import curve_fit
from joblib import Parallel, delayed
from kappa import rede
import KMC


#loading lattice
kappamean, d = rede(1000,filename="kappa.txt")




R_0   = 2    # R_F = k^2 * R_0
r     = 1      # distância entre as moléculas 
Tau   = 100
N     = 5000   #Número máximo de éxcitons 





resultados = Parallel(n_jobs=-1,verbose=10)(delayed(KMC.L_D_calc)(d,R_0,r,Tau,i) for i in range(1,N,int(N/10)))


lista_L_D_sig    = np.array(resultados)[:,0]#[resultados[i][0] for i in range(N-1)]
lista_L_D_steps  = np.array(resultados)[:,1]# [resultados[i][1] for i in range(N-1)]
lista_n          = np.arange(len(resultados))#[i for i in range(N-1)]



#SIMPLE TEST IN A SINGLE RUN EVENT
#STEPS    = np.mean(lista_L_D_steps)
#ANALYTIC = KMC.L_D_an(R_0,r)
#print(STEPS,ANALYTIC)
##################


#lista_L_D_sig = []
#lista_L_D_steps = []
#lista_n = []
#for n in range(1,N+1):
#    L_D = L_D_calc(d,R_0,r,Tau,n)
#    lista_L_D_sig.append(L_D[0])
#    lista_L_D_steps.append(L_D[1])
#    lista_n.append(n)

#print(lista_L_D_sig)
#print(lista_L_D_steps)

d = {"n_exc": lista_n.tolist(), "L_D_sig": lista_L_D_sig.tolist(), "L_D_steps": lista_L_D_steps.tolist()}

with open("dados_5.txt", "w") as arquivo:
    json.dump(d, arquivo)
arquivo.close()
