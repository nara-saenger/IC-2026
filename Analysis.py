
import random
import math
import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.optimize import curve_fit
from joblib import Parallel, delayed
from kappa import rede



loadingfile = "dados_5.txt"



N     = 400 
with open(loadingfile, "r") as arquivo:
    d = json.load(arquivo)
    
    
lista_L_D_steps = d["L_D_steps"]
lista_L_D_sig = d["L_D_sig"]
lista_n = d["n_exc"]

# Comprimento de difusão calculado pelo desvio padrão vs número de éxcitons

plt.figure()
plt.plot(lista_n, lista_L_D_sig, "o")
plt.ylabel(r"$L_D$ pelo desvio padrão")
plt.xlabel(r"Número de éxcitons")
plt.xlim(0,N+1)
plt.ylim(0,max(lista_L_D_sig)+1)
plt.show()

# Comprimento de difusão calculado pels passos vs número de éxcitons

plt.figure()
plt.plot(lista_n, lista_L_D_steps, "o")
plt.ylabel(r"$L_D$ pelo número de passos")
plt.xlabel(r"Número de éxcitons")
plt.xlim(0,N+1)
plt.ylim(0,max(lista_L_D_steps)+1)
plt.show()

# Comparação entre os comprimentos calculados 

plt.figure()
plt.plot(lista_L_D_sig, lista_L_D_steps, "o")
plt.ylabel(r"$L_D$ pelo número de passos")
plt.xlabel(r"$L_D$ pelo desvio padrão")
plt.xlim(0,max(lista_L_D_sig)+1)
plt.ylim(0,max(lista_L_D_steps)+1)
plt.show()
