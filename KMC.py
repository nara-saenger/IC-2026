
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
    s = [sum(p[:i+1]) for i in range(len(p))]
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
        kappa2_D = calc_kappa2(s_i,s_D)
        

        K_Emi = 1/Tau
        K_Es  = calc_K(kappa2_E,R_0,r,Tau)
        K_Dir = calc_K(kappa2_D,R_0,r,Tau)

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
            x_rel    -= r 
            contagem += 1
            traj.append(pos)
            
    # nota para a Nara:
    # x_rel mede a distância efetiva que a quasipartícula viaja em relação à sua pos.
    # Se for aplicar esse código para redes não uniformes, alguns cuidados:
    # Crie uma rede de tamanho N+1, onde o sítio N atuará como um sítio fantasma para intermediar os pulos de fronteira.
    # Ele não representará um sítio real. Será um artifício para induzir a condição de contorno.
    # Na prática, ele terá o seguinte comportamento:
    # - representa o sítio 0 quando o pulo for do N-1 para o 0
    # - representa o sítio N-1 quando o pulo for do 0 para o N-1
    #
    #
    #  . . .. . . . .    *
    #                    ^----- fronteira
    #  0 1 ....    N-1   N
    #
    #
    #  pulo do N-1 para o N vai levar para o sítio 0 com a distância   pos(N) - pos(N-1)              
    #  pulo do 0 para o N vai levar para o sítio N-1 com a distânica  pos(N-1)- pos(N)
    #
    #  você terá de tratar explicitamente essas duas linhas de condição de contorno. O resto seguirá a lógica usual.
    # 
    # Seja i = sítio antes do pulo, j = depois
    # x_rel += dist onde dist = boundary_cond_shift(i,j)
    #
    # boundary_cond_shift(i,j):
    #   se i==N-1 e j==N: <--- pulo na froteira
    #     shift = pos(N) - pos(N-1)
    #     x     = 0
    #   se i==0 e j==-1: <--- pulo na fronteira
    #     shift = pos(N-1) - pos(N)
    #     x     = N-1
    #   caso o contrário: <--- se o pulo for fora da fronteira
    #     shift = pos(j) - pos(i)
    #     x   = j  
    #  retorna shift,x
    #
    #
    # shift será para medir a distância efetiva que o exciton vai percorrer até sua fluorescência
    # x será a posição na rede (importante para localizar num eventual gráfico ou se a rede tiver regioês interessantes)
    #
    #
    #
    #
    # Aqui eu fiz com as distâncias, mas note que será necessário um remendo para o kappa.
    # você precisa forçar, na hora de calcular o kappa, a seguinte condição:
    # se o pulo for i==0 e j==N-1: dipolo_vizinho = dipolo(N-1)
    # se o pulo for i==N-1 e j==N: dipolo_vizinho = dipolo(0)
    # 
    # aí sim você pode chamar o calc_kappa2(s_i,dipolo_vizinho)

    return(traj,t,contagem,x_rel)


# Funções do comprimento de difusão 

def L_D_calc(d,R_0,r,Tau,n):     # para n éxcitons
    s = []
    s_N = []
    for i in range(n):
        t = traj(d,R_0,r,Tau)
        trajetoria = t[0] 
        passos = t[2]
        deltaX = t[3]
        
        s.append(deltaX)
        s_N.append(passos)

    
    N         = (sum(s_N)/n)
    L_D_steps = r*(N**0.5)
    # só é x**2 pq até o momento x_rel é medido partindo do 0. Se fosse Xdepois-Xantes, teria de compensar a distância Xdepois a cada volta na rede
    L_D_sig   = np.sqrt(np.mean([x**2 for x in s]))

    return L_D_sig, L_D_steps



def L_D_an(R_0,r):
    RF = ((2/3)**(1/6))*R_0
    L_D = (np.sqrt(2)*(RF)**3)/(r**2)  #para contagem de passos dos dois lados 
    return L_D



