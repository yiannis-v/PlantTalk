
# coding: utf-8

# # Photosynthesis parameters calculation from ACi curves

# This file will be using the calculations developped by Sharkey et al 2007.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline
from scipy.optimize import curve_fit

#define constants
gm = 0.874 #from excel file, check constant
T = 23 #temperature in celcius // eventually read this in
R = 0.008314 #ideal gas constant
O = 21 #oxygen conc. in air // assume 21
Kc = np.exp(35.9774 - (80.99/(R*(273 + T))))
Ko = np.exp(12.3772 - (23.72/(R*(273 + T))))
g_star = O*(np.exp(11.187 - (24.46/(R*(273 + T)))))/21

#read in A values, Ci values, Rd, Temp
data = pd.read_csv('/Users/JackieJ/Downloads/160607_soybean_16r.csv')
cleandata = data.loc[data.Obs!='Remark='] #now only rows with data points remain
A_light = cleandata.Photo[1:11] #light Photo
Ci_light = (101*0.001)*cleandata.Ci[1:11] #light Ci in Pa

#calculate other variables from data
Cc = Ci_light - A_light/gm
TPU = (A_light + Rd)/3
J = (A_light + Rd)*(4*Cc + 8*g_star)/(Cc - g_star)

#define function for A in terms of Cc

def func(x, a, g_star, Kc, O, Ko, Rd):
  return a*(x - g_star)/(x + Kc*(1 + O/Ko)) - Rd

popt, pcov = curve_fit(func, Cc, A_light)

plt.plot(func(Cc, *popt))
plt.scatter(Ci_light, A_light)
plt.ylabel('Photo')
plt.xlabel('Ci')
plt.show()








