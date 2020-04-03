# Project: ABE6035 Remote Sensing Project - Spring 2020
# Description: Code to calculate the apparent temperature (T_AP)
# Authors: Stephen Lantin, J. Barrett Carter
# Date Created: April 2, 2020

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.integrate as integrate

# Ka absorption coefficient for methane and CO2 (Ulaby & Long, Chapter 5)
# Theoretical model for Ka_CO2 in the troposphere can be found in: Absorption coefficient of carbon dioxide
# across atmospheric troposphere layer (Wei et al., 2018)

#Ka_CH4 = ...
#Ka_CO2 = ...

# U.S. Standard Atmosphere, 1976
def temperature(z):

    conds = [(z >= 0) & (z < 11000),
             (z >= 11000) & (z < 20000),
             (z >= 20000) & (z < 32000),
             (z >= 32000) & (z < 47000),
             (z >= 47000) & (z < 51000),
             (z >= 51000) & (z < 71000),
             (z >= 71000)]
    
    # slopes are in K/m, b term is the temperature at each boundary
    funcs = [lambda z: -0.0065*z+288.15,
             lambda z: 216.65,
             lambda z: 0.001*z+216.65,
             lambda z: 0.0028*z+228.65,
             lambda z: 270.65,
             lambda z: -0.0028*z+270.65,
             lambda z: -0.002*z+214.65]

    return np.piecewise(np.float(z),conds,funcs)

# # Calculation of upwelling and downwelling temperature

# def t_up(theta,H,kappa):
#     # optical thickness (probably have to assume a constant kappa, unless you find new data)
#     # if the kappa we find is dependent on z, then we'll have to do symbolic integration
#     # because tau contains an integration over z that feeds into the overall integration of z'.
#     result = (1/np.cos(theta))*integrate.quad(lambda zp: kappa*temperature(zp)*np.exp(-(kappa*H-kappa*zp)*1/np.cos(theta)),0,H)
#     return result

# def t_dn(theta,H,kappa):
#     result = (1/np.cos(theta))*integrate.quad(lambda zp: kappa*temperature(zp)*np.exp(-(kappa*zp)*1/np.cos(theta)),0,np.inf)
#     return result
print(temperature(0))
print(temperature(12000))
