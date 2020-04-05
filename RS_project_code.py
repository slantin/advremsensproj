# Project: ABE6035 Remote Sensing Project - Spring 2020
# Description: Code to calculate the apparent temperature (T_AP)

# This will produce radiative forcing (power density)

# Authors: Stephen Lantin, J. Barrett Carter
# Date Created: April 2, 2020

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.integrate as integrate

# Ka absorption coefficient for methane and CO2 (Ulaby & Long, Chapter 5)
# Theoretical model for Ka_CO2 in the troposphere can be found in: Absorption coefficient of carbon dioxide
# across atmospheric troposphere layer (Wei et al., 2018)

# EXPONENTIAL WIDE-BAND MODEL PARAMETERS

h = 6.626*10**(-34) # J*s
c = 299792458 # m/s
k = 1.381 * 10**(-23) #J/K
T0 = 100 # K

# CO2
co2_m = 3 # deg. of freedom (?)
co2_etas = [1351,667,2396] # cm^(-1)
co2_gs = [1,2,1]
co2_bands = [15,10.4,9.4,4.3,2.7,2.0]  # um
co2_band_centers = [667,960,1060,2410,3660,5200] # cm^(-1)
co2_d1s = [0,-1,-1,0,1,2]
co2_d2s = [1,0,0,0,0,0]
co2_d3s = [0,1,1,1,1,1] # see footnote b in Table 9.2
co2_bs = [1.3,1.3,1.3,1.3,1.3,1.3] # pressure parameter at T0 = 100K
co2_ns = [0.7,0.8,0.8,0.8,0.65,0.65] # pressure parameter at T0 = 100K
co2_alpha_0s = [19.0,2.47*10**9,2.47*10**9,110.0,4.0,0.066] # see footnote b in Table 9.2
co2_beta_0s = [0.06157,0.04017,0.04017,0.24723,0.13341,0.39305]
co2_omega_0s = [12.7,13.4,10.1,11.2,23.5,34.5]

# CH4
ch4_m = 4 # deg. of freedom (?)
ch4_etas = [2914,1526,3020,1306] # cm^(-1)
ch4_gs = [1,2,3,3]
ch4_bands = [7.66,3.31,2.37,1.71]  # um
ch4_band_centers = [1310,3020,4220,5861] # cm^(-1)
ch4_d1s = [0,0,1,1]
ch4_d2s = [0,0,0,1]
ch4_d3s = [0,1,0,0]
ch4_d4s = [1,0,1,1]
ch4_bs = [1.3,1.3,1.3,1.3,1.3,1.3] # pressure parameter at T0 = 100K
ch4_ns = [0.8,0.8,0.8,0.8] # pressure parameter at T0 = 100K
ch4_alpha_0s = [28.0,46.0,2.9,0.42]
ch4_beta_0s = [0.08698,0.06973,0.35429,0.13219]
ch4_omega_0s = [21.0,56.0,60.0,45.0]

# assign values to parameter matrix: band, band center, d1, d2, d3, b, n, alpha_0, beta_0, omega_0
co2 = np.zeros((6,10))
co2[:,0] = co2_bands
co2[:,1] = co2_band_centers
co2[:,2] = co2_d1s
co2[:,3] = co2_d2s
co2[:,4] = co2_d3s
co2[:,5] = co2_bs
co2[:,6] = co2_ns
co2[:,7] = co2_alpha_0s
co2[:,8] = co2_beta_0s
co2[:,9] = co2_omega_0s

# assign values to parameter matrix: band, band center, d1, d2, d3, d4, b, n, alpha_0, beta_0, omega_0
ch4 = np.zeros((4,11))
ch4[:,0] = ch4_bands
ch4[:,1] = ch4_band_centers
ch4[:,2] = ch4_d1s
ch4[:,3] = ch4_d2s
ch4[:,4] = ch4_d3s
ch4[:,5] = ch4_d4s
ch4[:,6] = ch4_bs
ch4[:,7] = ch4_ns
ch4[:,8] = ch4_alpha_0s
ch4[:,9] = ch4_beta_0s
ch4[:,10] = ch4_omega_0s

# def psi_co2(temperature,dof,matrix): # temperature, degrees of freedom, matrix
#     for k in range(0,dof):
#         for 



# alpha_0 = 
# def alpha(temperature,alpha_0):
#     alpha = alpha_0


#Ka_CH4 = ...
#Ka_CO2 = ...



# # U.S. Standard Atmosphere, 1976
# def temperature(z):

#     conds = [(z >= 0) & (z < 11000),
#              (z >= 11000) & (z < 20000),
#              (z >= 20000) & (z < 32000),
#              (z >= 32000) & (z < 47000),
#              (z >= 47000) & (z < 51000),
#              (z >= 51000) & (z < 71000),
#              (z >= 71000)]
    
#     # slopes are in K/m, b term is the temperature at each boundary
#     funcs = [lambda z: -0.0065*z+288.15,
#              lambda z: 216.65,
#              lambda z: 0.001*z+216.65,
#              lambda z: 0.0028*z+228.65,
#              lambda z: 270.65,
#              lambda z: -0.0028*z+270.65,
#              lambda z: -0.002*z+214.65]

#     return np.piecewise(np.float(z),conds,funcs)

# # # Calculation of upwelling and downwelling temperature


# # to the tropopause only (as per definition of radiative forcing)
# def t_up(theta,H,kappa):
#     # optical thickness (probably have to assume a constant kappa, unless you find new data)
#     # if the kappa we find is dependent on z, then we'll have to do symbolic integration
#     # because tau contains an integration over z that feeds into the overall integration of z'.
#     result = (1/np.cos(theta))*integrate.quad(lambda zp: kappa*temperature(zp)*np.exp(-(kappa*H-kappa*zp)*1/np.cos(theta)),0,H)
#     return result

# def t_dn(theta,H,kappa):
#     result = (1/np.cos(theta))*integrate.quad(lambda zp: kappa*temperature(zp)*np.exp(-(kappa*zp)*1/np.cos(theta)),0,np.inf)
# #     return result
# print(temperature(0))
# print(temperature(12000))


# integrate the spectral power density and integrate to get the overall power density and radiative forcing