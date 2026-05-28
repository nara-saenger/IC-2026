
import random
import math
import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.optimize import curve_fit


def calc_dt(k):
    u = random.random()
    return -np.log(u)/sum(k)

def K_to_prob(k):
    K_T = sum(k)
    return [k[i]/K_T for i in range(len(k))]

def which_event(p):
    s = np.cumsum(p)
    s = s/s[-1]
    u = random.random()
    if u < s[0]:
        return 0
    if u < s[1]:
        return 1
    if u <= s[2]:
        return 2 

def calc_kappa2(mu_d, mu_a):
    mu_d = np.array(mu_d)
    mu_a = np.array(mu_a)
    kappa = np.dot(mu_d,mu_a) - 3*(mu_d[0])*(mu_a[0])  # Rede na coordenada x
    return kappa**2

def calc_K(kappa2,R_0,r,Tau):
    RF = ((kappa2)**(1/6))*R_0
    K = ((RF/r)**6)/Tau
    return K
    
def cond_contorno(N,X_0):
    if X_0 > N:
        return 0
    if X_0 < 0:
        return N
    return X_0


#Função principal

def traj(d,R_0,r,Tau):
    N          = len(d)     #Tamanho da rede
    X_0        = random.randint(0,N-1) #posicao inicial
    alive      = True
    t          = 0
    n          = 0
    traj       = [X_0]
    x_rel      = 0 #distancia relativa
    pos        = X_0
    contagem   = 0

    while alive:
        x_E_con = cond_contorno(N-1,pos -1)
        x_D_con = cond_contorno(N-1,pos +1)
        s_i = d[pos]
        s_E = d[x_E_con]
        s_D = d[x_D_con]
        kappa2_E = calc_kappa2(s_i,s_E)
        kappa2_D= calc_kappa2(s_i,s_D)
        
        kappa2_E = 2/3 # se forcar o 2/3, LD sigma concorda
        kappa2_D = 2/3 # se forcar o 2/3, LD sigma concorda
        K_Emi = 1/Tau
        K_Es  = calc_K(kappa2_E,R_0,r,Tau)
        K_Dir = calc_K(kappa2_D,R_0,r,Tau)

        #print(np.sqrt(2*K_Es/K_Emi)*r)
        list_K = [K_Emi, K_Dir, K_Es]

        P = K_to_prob(list_K) 
        dt = calc_dt(list_K)
        t += dt
        event = which_event(P)

        if event ==0:
            alive = False
        if event == 1:    #Pulo para a direita
            pos        = x_D_con
            x_rel     += r
            contagem  += 1
            traj.append(pos)
        if event ==2:     #Pulo para a esquerda 
            pos       = x_E_con
            x_rel    -= r #+=r#-= r 
            contagem += 1
            traj.append(pos)

    return(traj,t,contagem,x_rel)


# Funções do comprimento de difusão 
def L_D_an(R_0,r):
    RF = ((2/3)**(1/6))*R_0
    L_D = (np.sqrt(2))*((RF)**3)/(r**2)  # np.sqrt(2) para contagem de passos dos dois lados 
    return L_D
    
def L_D_calc(d,R_0,r,Tau,n):     # para n éxcitons
    s = []
    s_N = []
    dx = []
    dx2 = []
    for i in range(n):
        t = traj(d,R_0,r,Tau)
        trajetoria = t[0] 
        passos = t[2]
        deltaX = t[3]
        #print(deltaX/r,passos, flush=True)
        s.append(deltaX)
        s_N.append(passos)
        dx.append(abs(deltaX))
        dx2.append(deltaX**2)
        
    dx  = np.array(dx)
    dx2 = np.array(dx2)
    dx_mean = np.mean(dx)
    dx2_mean  = np.mean(dx2)
    
    N         = (sum(s_N)/n)
    L_D_steps = r*(N**0.5)
    L_D_sig    = np.sqrt(2)*np.sqrt(dx2_mean-dx_mean**2)
    L_D_sig2   = np.sqrt(2)*np.sqrt(np.var(dx))
    L_D_ana     = L_D_an(R_0,r)
    

    #print(L_D_steps,L_D_ana,L_D_sig,L_D_sig2,flush=True)
    return L_D_sig, L_D_steps






