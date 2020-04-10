# Project: ABE6035 Remote Sensing Project - Spring 2020
# Description: Code to calculate the apparent temperature (T_AP)

# This will produce radiative forcing (power density)

# Authors: Stephen Lantin, J. Barrett Carter
# Date Created: April 2, 2020

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt

# Ka absorption coefficient for methane and CO2 (Ulaby & Long, Chapter 5)
# Theoretical model for Ka_CO2 in the troposphere can be found in: Absorption coefficient of carbon dioxide
# across atmospheric troposphere layer (Wei et al., 2018)

# EXPONENTIAL WIDE-BAND MODEL PARAMETERS

h = 6.626*10**(-34) # J*s
##c = 299792458 # m/s
c = 2.998e10 #cm/s
kB = 1.381 * 10**(-23) #J/K
t_0 = 100 # K

# CO2
co2_m = [3,0,0,0,0,0] # deg. of freedom (?)
co2_eta1s = [1351,1351,1351,1351,1351,1351] # cm^(-1)
co2_eta2s = [667,667,667,667,667,667] # cm^(-1)
co2_eta3s = [2396,2396,2396,2396,2396,2396] # cm^(-1)
co2_g1s = [1,1,1,1,1,1]
co2_g2s = [2,2,2,2,2,2]
co2_g3s = [1,1,1,1,1,1]
co2_bands = [15,10.4,9.4,4.3,2.7,2.0]  # um
co2_band_centers = [667,960,1060,2410,3660,5200] # cm^(-1)
co2_d1s = [0,-1,-1,0,1,2]
co2_d2s = [1,0,0,0,0,0]
co2_d3s = [0,1,1,1,1,1] # see footnote b in Table 9.2
co2_bs = [1.3,1.3,1.3,1.3,1.3,1.3] # pressure parameter at t_0 = 100K
co2_ns = [0.7,0.8,0.8,0.8,0.65,0.65] # pressure parameter at t_0 = 100K
co2_alpha_0s = [19.0,2.47*10**9,2.48*10**9,110.0,4.0,0.066] # see footnote b in Table 9.2
co2_beta_0s = [0.06157,0.04017,0.11888,0.24723,0.13341,0.39305]
co2_omega_0s = [12.7,13.4,10.1,11.2,23.5,34.5]

# CH4
ch4_m = [4,0,0,0] # deg. of freedom (?)
ch4_etas = [2914,1526,3020,1306] # cm^(-1)
ch4_eta1s = [2914,2914,2914,2914]
ch4_eta2s = [1526,1526,1526,1526]
ch4_eta3s = [3020,3020,3020,3020]
ch4_eta4s = [1306,1306,1306,1306]
ch4_gs = [1,2,3,3]
ch4_g1s = [1,1,1,1]
ch4_g2s = [2,2,2,2]
ch4_g3s = [3,3,3,3]
ch4_g4s = [3,3,3,3]
ch4_bands = [7.66,3.31,2.37,1.71]  # um
ch4_band_centers = [1310,3020,4220,5861] # cm^(-1)
ch4_d1s = [0,0,1,1]
ch4_d2s = [0,0,0,1]
ch4_d3s = [0,1,0,0]
ch4_d4s = [1,0,1,1]
ch4_bs = [1.3,1.3,1.3,1.3] # pressure parameter at t_0 = 100K
ch4_ns = [0.8,0.8,0.8,0.8] # pressure parameter at t_0 = 100K
ch4_alpha_0s = [28.0,46.0,2.9,0.42]
ch4_beta_0s = [0.08698,0.06973,0.35429,0.68598]
ch4_omega_0s = [21.0,56.0,60.0,45.0]

# assign values to parameter matrix: band, band center, d1, d2, d3, b, n, alpha_0, beta_0, omega_0
co2 = np.zeros((6,17))
co2[:,0] = co2_m
co2[:,1] = co2_bands
co2[:,2] = co2_band_centers
co2[:,3] = co2_d1s
co2[:,4] = co2_d2s
co2[:,5] = co2_d3s
co2[:,6] = co2_bs
co2[:,7] = co2_ns
co2[:,8] = co2_alpha_0s
co2[:,9] = co2_beta_0s
co2[:,10] = co2_omega_0s
co2[:,11] = co2_eta1s
co2[:,12] = co2_eta2s
co2[:,13] = co2_eta3s
co2[:,14] = co2_g1s
co2[:,15] = co2_g2s
co2[:,16] = co2_g3s


# assign values to parameter matrix: band, band center, d1, d2, d3, d4, b, n, alpha_0, beta_0, omega_0
ch4 = np.zeros((4,20))
ch4[:,0] = ch4_m
ch4[:,1] = ch4_bands
ch4[:,2] = ch4_band_centers
ch4[:,3] = ch4_d1s
ch4[:,4] = ch4_d2s
ch4[:,5] = ch4_d3s
ch4[:,6] = ch4_d4s
ch4[:,7] = ch4_bs
ch4[:,8] = ch4_ns
ch4[:,9] = ch4_alpha_0s
ch4[:,10] = ch4_beta_0s
ch4[:,11] = ch4_omega_0s
ch4[:,12] = ch4_eta1s
ch4[:,13] = ch4_eta2s
ch4[:,14] = ch4_eta3s
ch4[:,15] = ch4_eta4s
ch4[:,16] = ch4_g1s
ch4[:,17] = ch4_g2s
ch4[:,18] = ch4_g3s
ch4[:,19] = ch4_g4s


# now that we have psi, we can use that to calculate alpha

def psi_fun(temperature,matrix): # temperature, degrees of freedom, matrix
    psi = np.zeros((np.size(matrix,0),1))
    iterations = int(100)
    dof = int(matrix[0,0])
    for p in range(0,np.size(matrix,0)):

        psi_num_prev = 1 #previous numerator of equation for psi
        psi_den_prev = 1 #previous denominator of equation for psi

        for k in range(0,dof):

            psi_num = 0 #numerator of equation for psi
            psi_den = 0 #denominator of equation for psi

            if matrix[p,k+3]<0:
                nu_0 = int(abs(matrix[p,k+3]))
            else:
                nu_0 = int(0)

            u = h*c*matrix[p,k+8+dof]/(kB*temperature)

            for nu in range((nu_0),int(nu_0+iterations)):

                psi_num = psi_num + math.factorial(nu+matrix[p,k+8+2*dof]+abs(matrix[p,k+3])-1)/\
                        (math.factorial(matrix[p,k+8+2*dof]-1)*math.factorial(nu))*np.exp(-1*u*nu)

            for nu in range(0,int(iterations)):
            
                psi_den = psi_den + math.factorial(nu+matrix[p,k+8+2*dof]-1)/\
                        (math.factorial(matrix[p,k+8+2*dof]-1)*math.factorial(nu))*np.exp(-1*u*nu)

            psi_num = psi_num*psi_num_prev
            psi_den = psi_den*psi_den_prev

            psi_num_prev = psi_num
            psi_den_prev = psi_den

                
        psi[p,0] = psi_num/psi_den

    
    return psi

def phi_fun(temperature,matrix): # temperature, degrees of freedom, matrix
    phi = np.zeros((np.size(matrix,0),1))
    iterations = int(100)
    dof = int(matrix[0,0])
    for p in range(0,np.size(matrix,0)):

        phi_num_prev = 1 #previous numerator of equation for phi
        phi_den_prev = 1 #previous denominator of equation for phi

        for k in range(0,dof):

            phi_num = 0 #numerator of equation for phi
            phi_den = 0 #denominator of equation for phi

            if matrix[p,k+3]<0:
                nu_0 = int(abs(matrix[p,k+3]))
            else:
                nu_0 = int(0)

            u = h*c*matrix[p,k+8+dof]/(kB*temperature)

            for nu in range((nu_0),int(nu_0+iterations)):

                phi_num = phi_num + (math.factorial(nu+matrix[p,k+8+2*dof]+abs(matrix[p,k+3])-1)/\
                        (math.factorial(matrix[p,k+8+2*dof]-1)*math.factorial(nu))*np.exp(-1*u*nu))**(0.5)

##            for nu in range(0,int(iterations)):
            
                phi_den = phi_den + math.factorial(nu+matrix[p,k+8+2*dof]+abs(matrix[p,k+3])-1)/\
                        (math.factorial(matrix[p,k+8+2*dof]-1)*math.factorial(nu))*np.exp(-1*u*nu)

            phi_num = phi_num*phi_num_prev
            phi_den = phi_den*phi_den_prev

            phi_num_prev = phi_num
            phi_den_prev = phi_den

                
        phi[p,0] = phi_num**2/phi_den

    
    return phi

def alpha(temperature,matrix): #returns alpha in m^2/(g*cm)
    global psi_fun

    alphas = np.zeros((np.size(matrix,0),1)) # array initialization
    psis = alphas
    psi_0s = alphas
    dof = int(matrix[0,0])
    psis = psi_fun(temperature,matrix)
    psi_0s = psi_fun(t_0,matrix) # psi function outside of the loop --> decreased computational time
#    print(psis)
#    print(psi_0s)
    for p in range(0,np.size(matrix,0)):

        sum_num = 0 # summation term in the numerator
        sum_den = 0 # summation term in the denominator

        for k in range(0,dof):

            u = h*c*matrix[p,k+8+dof]/(kB*temperature)
            u_0 = h*c*matrix[p,k+8+dof]/(kB*t_0)
            sum_num += u*matrix[p,k+3]
            sum_den += u_0*matrix[p,k+3]

        alphas[p] = matrix[p,dof+5]*(1-(np.exp(-sum_num)*psis[p]))/\
                (1-(np.exp(-sum_den)*psi_0s[p])) # this will take the last value of sum num from each k loop
    return alphas
#print(co2[:,8])
print(alpha(500,co2))

def beta(temperature,matrix):
    global phi_fun
    betas = np.zeros((np.size(matrix,0),1)) # array initialization
    dof = int(matrix[0,0])
    phis = phi_fun(temperature,matrix)
    phi_0s = phi_fun(t_0,matrix)
    for p in range(0,np.size(matrix,0)):

        betas[p] = matrix[p,dof+6]*(t_0/temperature)**(1/2)*phis[p]/phi_0s[p]
    
    return betas

#print(beta(500,co2))

def omega(temperature,matrix):
    omegas = np.zeros((np.size(matrix,0),1))
    dof = int(matrix[0,0])
    for p in range(0,np.size(matrix,0)):
        omegas[p] = matrix[p,dof+7]*(temperature/t_0)**(1/2)
    return omegas

#print(omega(500,co2))