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


class DustyGasModelZhou:
    def __init__(self, Bg, c_ch, M, mu, i, L, T,
                 D_binaryEff, D_knudsenEff):
        # self.epsilon = epsilon
        # self.tau = tau
        # self.rp = rp
        self.c_ch = c_ch
        self.M = M
        self.muMix = (c_ch*mu).sum()
        F = ct.faraday*1.e-3
        self.J = np.array([i / 2 / F, - i/2/F])
        self.z = L
        self.T = T
        self.cT = c_ch.sum()
        self.x = c_ch/c_ch.sum()
        self.N = 2
        self.D_binaryEff = D_binaryEff
        self.D_knudsenEff = D_knudsenEff
        self.Bg = Bg
                

    def calculate_D_DGM(self):
        H = np.zeros((self.N, self.N))
        sumTerm = 0.
        for k in range(self.N):
            for l in range(self.N):
                for j in range(self.N):
                    if j != k:
                        sumTerm += self.x[j] / self.D_binaryEff[k, j]
                if k == l:
                    H[k, l] = ((1/self.D_knudsenEff[k] + sumTerm))
                else:
                    H[k, l] = -self.x[k] / self.D_binaryEff[k, l]
                sumTerm = 0.
        return inv(H)
    
    # define the functioin
    def funcZhou(self, x, P, DGM_kl):
        X1_tpb, X2_tpb, dp = x
        D_knudsenEff = self.D_knudsenEff
        X_ch = self.c_ch
        dz = self.z
        Bg = self.Bg
        muMix = self.muMix
        J = self.J
        
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
        return np.array([f1, f2, f3])


    def solveDGM(self):
        DGM_kl = self.calculate_D_DGM()
        P = ct.one_atm
        c1, c2, dp_zhou = fsolve(self.funcZhou, 
                                 [cT*0.1, cT*0.9, 10000], 
                                  args = (P, DGM_kl)
                                  )
        x_zhou = np.array([c1, c2]) / cT
        w_zhou = x_zhou * M / np.sum(x_zhou * M)
        P_zhou = P + dp_zhou
        print(f'Zhou Compostion by weight @tpb \t[H2, H20] is {w_zhou}')
        print(f'The total pressure is {P_zhou*1e-5:.4} bar')
        return x_zhou, P_zhou


def permeabilityFactorBg(epsilon, tau, rp):
    """
    This calculates the permeability factor commonly noted as Bg
    Here the Kozeny-Carman relationship is used, according to the
    Zhou paper.
    
    """
    # 5.37e-17
    return  (epsilon**3 * (2*rp)**2 / 72 / tau / (1 - epsilon)**2)

if __name__ =="__main__":

    # operational parameters
    D_knudsenEff = np.array([2.8e-5, 9.4e-6])
    D_binaryEff = np.array([[0., 8.e-4], [8.e-4, 0.]])
    i = 24999.998
    w = np.array([wH2 := 0.03148458, 1 - wH2])


    # universal constants
    T = 700 + 273.15
    P = ct.one_atm
    R = ct.gas_constant*1.e-3
    F = ct.faraday*1.e-3
    cT = P/R/T

    # geometry
    epsilon = 0.3       # porosity of the electrode (affects mass balance; smaller -> worse)
    tau = 5.            # tortuosity of the electrode
    rp = 5.e-7          # radius of the pores (affects mass balance; smaller -> worse)
    dEle = 200e-6       # thickness of the electrode (affects mass balance a lot; bigger -> worse)
    A = 100             # area in cm²
 
    # gas variables
    M = np.array([2e-3, 18e-3])
    x = w / M / np.sum(w / M)
    c_ch = x * cT
    mu = np.array([1.84e-5, 3.26e-5])    # dynamic visosities at 600°C
    muMix = (c_ch*mu).sum()                # dynamic for the mix

    Bg = permeabilityFactorBg(epsilon, tau, rp)
    dgm = DustyGasModelZhou(Bg, c_ch, M, mu, i, dEle, T, 
                            D_binaryEff, D_knudsenEff)
    
    DGM_kl = dgm.calculate_D_DGM()
    x_tpb, P_tpb = dgm.solveDGM()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# muss neu berechnet werden

    # def calcAJ(X_ch, J, D_binaryEff, D_knudsenEff, XT):
    #     Akk = np.zeros((N, N))
    #     Akl = np.zeros((N, N))
    #     sumTerm = 0
    #     for k in range(N):
    #         for l in range(N):
    #             if l != k:
    #                 sumTerm += X_ch[l] / D_binaryEff[k, l]  
    #                 Akl[k, l] = - X_ch[k] / XT / D_binaryEff[k, l]
    #             else:
    #                 Akk[k, k] = 1 / D_knudsenEff[k] + 1/XT * sumTerm
    #                 sumTerm = 0.
    #     A = Akk + Akl
    #     return (A @ J)
    
    # %% Expression from Zhu Paper
    
        
    # %%
    
    
    
    # def func(x, X_ch, dz, Bg, P):
    #     X1_tpb, X2_tpb, dp = x
    #     X1_ch, X2_ch = X_ch
    #     AJ = calcAJ(X_ch, J, D_binaryEff, D_knudsenEff, (P + dp)/R/T)
    #     f1 = AJ[0] - (-(X1_tpb-X1_ch)/dz - X1_ch*Bg/D_knudsenEff[0]*dp/dz)
    #     f2 = AJ[1] - (-(X2_tpb-X2_ch)/dz - X2_ch*Bg/D_knudsenEff[1]*dp/dz)
    #     f3 = (P + dp)/R/T  - X1_tpb - X2_tpb
    #     return [f1, f2, f3]
    
    
    # dz = np.array([25e-6, 150e-6, 25e-6])
    # L = 200e-6
    # i = 10 
    # dz = np.zeros((i)) + L / i
    # XX = X_ch
    # dp_total = 0
    # for z in dz:
    #     XX[0], XX[1], dp = fsolve(func, [XT*0.1, XT*0.9, 10000], 
    #                               args = (XX, z, Bg, P+dp_total))
    #     dp_total += dp
    # X1_tpb, X2_tpb = XX
    
        
    # x_tpb = np.array([X1_tpb, X2_tpb]) / XT
    # w_tpb = x_tpb * M / np.sum(x_tpb * M)
    
    # %%
    # print(f'Compostion by weight @channel \t[H2, H20] is {w}')
    # print(f'Compostion by weight @tpb     \t[H2, H20] is {w_tpb}')
    # # print(f'The pressure difference is {dp_total:.5} Pa,')
    # print(f'the total pressure is {(P + dp)*1e-5:.5} bar')
        
        
        
        
        
        
        
        
        
        
        
    
    
