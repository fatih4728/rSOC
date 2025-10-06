# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 10:08:54 2025

Matrix according to the Kogekar Paper

@author: smfadurm
"""

import numpy as np
import cantera as ct
from scipy.optimize import fsolve
from numpy.linalg import inv

N = 2

T = 700 + 273.15
P = ct.one_atm
R = ct.gas_constant*1.e-3
M = np.array([2e-3, 18e-3])
F = ct.faraday*1.e-3

XT = P/R/T

Bg = 5.37E-17
wH2 = 0.0314185
w = np.array([wH2, 1-wH2])
x = w / M / np.sum(w / M)


X_ch = x * XT

D_knudsenEff = np.array([2.8e-5, 9.4e-6])
# D_knudsenEff = np.array([9.4e-6, 2.8e-5])
D_binaryEff = np.array([[0., 8.e-4], [8.e-4, 0.]])

i = 24999.7155
J = np.array([i / 2 / F, - i/2/F])


Akk = np.zeros((N, N))
Akl = np.zeros((N, N))

# muss neu berechnet werden

def calcAJ(X_ch, J, D_binaryEff, D_knudsenEff, XT):
    sumTerm = 0
    for k in range(N):
        for l in range(N):
            if l != k:
                sumTerm += X_ch[l] / D_binaryEff[k, l]  
                Akl[k, l] = - X_ch[k] / XT / D_binaryEff[k, l]
            else:
                Akk[k, k] = 1 / D_knudsenEff[k] + 1/XT * sumTerm
                sumTerm = 0.
    A = Akk + Akl
    return (A @ J)
AJ_test = calcAJ(x, J, D_binaryEff, D_knudsenEff, XT)

# %% Expression from Zhu Paper
H = np.zeros((N, N))
sumTerm = 0.
# calculate the D_DGM Matrix, by calculating H and inversing it
for k in range(N):
    for l in range(N):
        for j in range(N):
            if j != k:
                sumTerm += x[j] / D_binaryEff[k, j]
                # print(sumTerm)
        if k == l:
            H[k, l] = ((1/D_knudsenEff[k] + sumTerm))
        else:
            H[k, l] = -x[k] / D_binaryEff[k, l]
        sumTerm = 0.
DGM_kl = inv(H)
mu = np.array([1.84e-5, 3.26e-5])    # dynamic visosities at 600°C
muMix = (X_ch*mu).sum()                # dynamic for the mix

# define the functioin
def funcZhou(x, X_ch, dz, Bg, P, DGM_kl, muMix):
    X1_tpb, X2_tpb, dp = x
    # X1_ch, X2_ch = X_ch
    f1 = ( -J[0] - (DGM_kl[0, 0] * (X1_tpb - X_ch[0]) / dz 
                    + DGM_kl[0, 1] * (X2_tpb - X_ch[1]) / dz) -
          (DGM_kl[0, 0]*X_ch[0]/D_knudsenEff[0] + 
           DGM_kl[0, 1]*X_ch[1]/D_knudsenEff[1])
          * Bg / muMix * dp/dz)
    f2 = ( -J[1] - (DGM_kl[1, 0] * (X1_tpb - X_ch[0]) / dz 
                    + DGM_kl[1, 1] * (X2_tpb - X_ch[1]) / dz) -
          (DGM_kl[1, 0]*X_ch[0]/D_knudsenEff[0] + 
           DGM_kl[1, 1]*X_ch[1]/D_knudsenEff[1])
          * Bg / muMix * dp/dz)
    f3 = (P + dp)/R/T  - X1_tpb - X2_tpb
    # print("f1:", f1, type(f1))
    return np.array([f1, f2, f3])

z = 200e-6
X1, X2, dp_zhou = fsolve(funcZhou, [XT*0.1, XT*0.9, 10000], 
                          args = (X_ch, z, Bg, P, DGM_kl, muMix))
x_zhou = np.array([X1, X2]) / XT
w_zhou = x_zhou * M / np.sum(x_zhou * M)
P_zhou = P + dp_zhou
# print(f'Zhou Compostion by weight @tpb     \t[H2, H20] is {w_zhou}')
# print(P_zhou*1e-5)

# %%



def func(x, X_ch, dz, Bg, P):
    X1_tpb, X2_tpb, dp = x
    X1_ch, X2_ch = X_ch
    AJ = calcAJ(X_ch, J, D_binaryEff, D_knudsenEff, (P + dp)/R/T)
    f1 = AJ[0] - (-(X1_tpb-X1_ch)/dz - X1_ch*Bg/D_knudsenEff[0]*dp/dz)
    f2 = AJ[1] - (-(X2_tpb-X2_ch)/dz - X2_ch*Bg/D_knudsenEff[1]*dp/dz)
    f3 = (P + dp)/R/T  - X1_tpb - X2_tpb
    return [f1, f2, f3]


# dz = np.array([25e-6, 150e-6, 25e-6])
L = 200e-6
i = 10 
dz = np.zeros((i)) + L / i
XX = X_ch
dp_total = 0
for z in dz:
    XX[0], XX[1], dp = fsolve(func, [XT*0.1, XT*0.9, 10000], 
                              args = (XX, z, Bg, P+dp_total))
    dp_total += dp
X1_tpb, X2_tpb = XX

    
x_tpb = np.array([X1_tpb, X2_tpb]) / XT
w_tpb = x_tpb * M / np.sum(x_tpb * M)

# %%
print(f'Compostion by weight @channel \t[H2, H20] is {w}')
print(f'Compostion by weight @tpb     \t[H2, H20] is {w_tpb}')
# print(f'The pressure difference is {dp_total:.5} Pa,')
print(f'the total pressure is {(P + dp)*1e-5:.5} bar')
    
    
    
    
    
    
    
    
    
    
    
    
    
    







