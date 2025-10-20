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
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI



class DustyGasModelZhou:
    def __init__(self, Bg, c_ch, M, mu, i, j, L, T,
                 D_binaryEff, D_knudsenEff):
        # self.epsilon = epsilon
        # self.tau = tau
        # self.rp = rp
        self.c_ch = c_ch
        self.M = M
        self.muMix = (c_ch*mu).sum()
        # F = ct.faraday*1e-3
        # self.J = np.array([i / 2 / F, - i/2/F])
        self.J = j
        self.z = L
        self.T = T
        self.cT = c_ch.sum()
        self.x_ch = c_ch/c_ch.sum()
        self.N = 2
        self.D_binaryEff = D_binaryEff
        self.D_knudsenEff = D_knudsenEff
        self.Bg = Bg
        self.R = ct.gas_constant*1.e-3
        self.P = ct.one_atm
                

    def calculate_D_DGM(self, x_ch):
        """
        Generate the DGM Matrix

        Parameters
        ----------
        x_ch : float
            Molar concentration of the species at channel side.

        Returns
        -------
        matrix
            The DGM matrix, which is the inverse of H.

        """
        H = np.zeros((self.N, self.N))
        sumTerm = 0.
        for k in range(self.N):
            for l in range(self.N):
                for j in range(self.N):
                    if j != k:
                        sumTerm += x_ch[j] / self.D_binaryEff[k, j]
                if k == l:
                    H[k, l] = ((1/self.D_knudsenEff[k] + sumTerm))
                else:
                    H[k, l] = -x_ch[k] / self.D_binaryEff[k, l]
                sumTerm = 0.
        return inv(H)
    
    # define the functioin
    def funcZhou(self, x, P, DGM_kl, dz, c_ch):
        """
        This is the function that will be solved. The equation for the
        DGM so to speak. The three boundary conditions are (1) and (2) 
        that the flux of H2 and H2O are respected and (3) that the gasses
        add up to a total concentration.

        Parameters
        ----------
        x : Array
            This is, what will be solved for.
        P : float
            Total pressure.
        DGM_kl : Matrix
            The matrix from the function.
        dz : float
            the length orthogonal to the electrolyte, that will be passed.
        c_ch : float
            concentration in the channel side increment.

        Returns
        -------
        array
            the electrolyte side concentration and the pressure.

        """
        c1_tpb, c2_tpb, dp = x
        D_knudsenEff = self.D_knudsenEff
        Bg = self.Bg
        muMix = self.muMix
        J = self.J
        
        f1 = ( -J[0] - (DGM_kl[0, 0] * (c1_tpb - c_ch[0]) / dz 
                        + DGM_kl[0, 1] * (c2_tpb - c_ch[1]) / dz) -
              (DGM_kl[0, 0]*c_ch[0]/D_knudsenEff[0] + 
               DGM_kl[0, 1]*c_ch[1]/D_knudsenEff[1])
              * Bg / muMix * dp/dz)
        f2 = ( -J[1] - (DGM_kl[1, 0] * (c1_tpb - c_ch[0]) / dz 
                        + DGM_kl[1, 1] * (c2_tpb - c_ch[1]) / dz) -
              (DGM_kl[1, 0]*c_ch[0]/D_knudsenEff[0] + 
               DGM_kl[1, 1]*c_ch[1]/D_knudsenEff[1])
              * Bg / muMix * dp/dz)
        f3 = (P + dp*1.)/self.R/self.T  - c1_tpb - c2_tpb
        return np.array([f1, f2, f3])


    def solveDGM(self):
        x_ch = self.x_ch
        DGM_kl = self.calculate_D_DGM(x_ch)
        P_zhou = ct.one_atm
        dp_sum = 0.
        i = 1000
        delta_z = np.zeros((i)) + self.z / i
        c_ch = self.c_ch
        for dz in delta_z:
            c1, c2, dp_zhou = fsolve(self.funcZhou, 
                                     [self.cT*0.1, self.cT*0.9, 1000], 
                                      args = (P_zhou, DGM_kl, dz, c_ch)
                                      )
            c_ch = np.array([c1, c2])
            P_zhou += dp_zhou*1.
            dp_sum += dp_zhou 
            # x_ch = c_ch / P_zhou * self.R * self.T
            x_ch = c_ch / ct.one_atm * self.R * self.T
            DGM_kl = self.calculate_D_DGM(x_ch)

        x_zhou = np.array([c1, c2]) / P_zhou * self.R*self.T
        
        return x_zhou, P_zhou

class ControlVolume:
    def __init__(self, mDotIn, mDotDiff, wIn):
        self.mDotIn = mDotIn
        self.mDotDiff = mDotDiff
        self.wIn = wIn
    
    def massBalance(self):
        mDotOut = self.mDotDiff + self.mDotIn
        xInH2, xInH2O = self.wIn
        if mDotOut.any() < 0:
            print("the outlet is not outletting")
        mDotH2_out = abs(self.mDotIn[0]) - abs(self.mDotDiff[0])
        yOutH2 = mDotH2_out/mDotOut.sum()
        
        Uf = 1 - yOutH2/self.wIn[0]
        if Uf.any() > 1.0:
            print("Not enough H2 to support operation")
            print("increase volume flow, or decrease current density or area")
        return mDotOut, np.array([yOutH2, 1-yOutH2]), Uf


def getDensities(components, P, T):
    rho = []
    for i in components:
        rho_x = PropsSI('D', 'P', P, 'T', T, i)
        rho.append(rho_x)
    return np.array(rho)

def getViscosities(components, P, T):
    mu = []
    for i in components:
        mu_x = PropsSI('V', 'P', P, 'T', T, i)
        mu.append(mu_x)
    return np.array(mu)

def permeabilityFactorBg(epsilon, tau, rp):
    """
    This calculates the permeability factor commonly noted as Bg
    Here the Kozeny-Carman relationship is used, according to the
    Zhou paper.
    
    """
    # 5.37e-17
    return  (epsilon**3 * (2*rp)**2 / 72 / tau / (1 - epsilon)**2)

def calculateMolarFlux(i, A):
    """
    This function calculates the anode molar flux of a SOC
    Only valid for pure hydrogen and steam

    Parameters
    ----------
    i : float
        current density in A/cm².
    A : Area
        Area in cm².

    Returns
    -------
    J : Array
        molar flux of H2 and H2O.

    """
    F = 96485.33212331001          # faraday constant
    ndotH2 = i/F/2*A
    ndotH2O = ndotH2
    return np.array((-ndotH2, ndotH2O))

def x2w(x, M):
    return x * M / np.sum(x*M)

def w2x(w, M):
    return w / M / np.sum(w / M)

def D_Knudsen(rp, T, M):
    """
    Source is Transport Phenomena in materials processing.
    !The equation works if rp is given in m. Might be a mistake hidden!

    Parameters
    ----------
    rp : float
        pore radius / cm.
    T : float
        temperature / K.
    M : float
        molar mass / g/mol.

    Returns
    -------
    float
        Knudsen Diffusion Coefficient.

    """
    return 9700*rp*(T/M)**0.5 *1e-2

def D_Fuller(T, P = 1.):
    """    
    Fuller equation for the binary diffusion coefficient

    Parameters
    ----------
    T : float
        Temperature / K.
    P : float, optional
        Pressure / atm. The default is 1..

    Returns
    -------
    Binary diffusion coefficient.

    """
    M_A = 2.01588
    M_B = 18.01528
    v_A = 7.07
    v_B = 12.7
    return (1e-3*T**1.75*(1/M_A+1/M_B)**0.5/
            (P*(v_A**(1/3)+v_B**(1/3))**2)) * 1e-4

if __name__ =="__main__":

    # operational parameters
    D_knudsenEff = np.array([2.8e-5, 9.4e-6])
    D_binaryEff = np.array([[0., 8.e-4], [8.e-4, 0.]])
    i = 24999.998
    w = np.array([wH2 := 0.031418, 1 - wH2])


    # universal constants
    T = 700 + 273.15
    P = ct.one_atm
    R = ct.gas_constant*1.e-3
    F = ct.faraday*1.e-3
    cT = P/R/T

    # geometry
    epsilon = 0.54       # porosity of the electrode (affects mass balance; smaller -> worse)
    tau = 4.81            # tortuosity of the electrode
    rp = 5.e-7          # radius of the pores (affects mass balance; smaller -> worse)
    dEle = 200e-6       # thickness of the electrode (affects mass balance a lot; bigger -> worse)
    A = 100             # area in cm²
 
    # gas variables
    M = np.array([2e-3, 18e-3])
    x = w2x(w, M)
    c_ch = x * cT
    mu = np.array([1.84e-5, 3.26e-5])    # dynamic visosities at 600°C

    Bg = permeabilityFactorBg(epsilon, tau, rp)
    Bg = 5.3e-17
    dgm = DustyGasModelZhou(Bg, c_ch, M, mu, i, dEle, T, 
                            D_binaryEff, D_knudsenEff)
    
    # DGM_kl = dgm.calculate_D_DGM()
    x_tpb, P_tpb = dgm.solveDGM()
    
    
    w_zhou = x2w(x_tpb, M)
    # print(f'Zhou Compostion by weight @tpb \t[H2, H20] is {w_zhou}')
    # %%
    print(f'Compostion by weight @channel \t[H2, H20] is {w}')
    print(f'Compostion by weight @tpb     \t[H2, H20] is {w_zhou}')
    # print(f'The pressure difference is {dp_total:.5} Pa,')
    print(f'the total pressure is {(P_tpb)*1e-5:.5} bar')
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    
       
        
        
        
        
        
        
        
        
        
        
    
    
