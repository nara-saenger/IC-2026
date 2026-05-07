
import random
import math
import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.optimize import curve_fit



def rede(N,filename=None):
    lista = []
    for i in range(int(N)):
        x = np.random.uniform(-1,1)
        y = np.random.uniform(-1,1)
        z = np.random.uniform(-1,1)
        norm = (x**2 + y**2 + z**2)**0.5
        x = x/norm
        y = y/norm
        z = z/norm
        lista.append([round(x,5),round(y,5),round(z,5)])
        
    kappas=[]
    array_mu = np.array(lista)
    for i in range(N-1):
        mu_i = array_mu[i]
        mu_iplus1= array_mu[i+1]
        kappa = np.dot(mu_i,mu_iplus1) - 3*(mu_i[0])*(mu_iplus1[0])
        kapppa2=kappa**2
        kappas.append(kapppa2)
    kappamean = np.mean(kappas)
    print(f'For this sample, <k^2> = {np.mean(kappas):.4f}, std = {np.std(kappas):.4f}')
    
    if filename is not None:
        with open(filename, "w") as arquivo: #saving 
            json.dump(lista, arquivo)
    return kappamean,lista



N_rede = 1000 #(tamanho da rede)





def testing():
  #testing convergence:
  means = []
  N = 1000
  for _ in range(N):
      kappamean,d = rede(N_rede)
      means.append(kappamean)
  plt.plot(np.arange(N),means)
  plt.savefig("kappa.png")
