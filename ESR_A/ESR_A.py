#%% init

import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from scipy.optimize import curve_fit as cfit

from uncertainties import unumpy
from uncertainties import ufloat
from uncertainties.umath import *
from uncertainties.unumpy import nominal_values as uval
from uncertainties.unumpy import std_devs as uerr


# helper function for plotting data and regression
def one4all(xdata,ydata,yerr=0,xerr=0,mode="none",f=None,xlabel="x",ylabel="y",show=True):
    fig = plt.figure(dpi=300)
    plt.errorbar(xdata,ydata,yerr,xerr,"o",label="Data")

    
    if mode == "none":
        fit= []
        
        
    
    elif mode == "linear":
        fit = linregress(xdata,ydata)
        f = lambda x,a,b: a*x+b
        #fit_label = "Regression: y=" + str(fit.slope) + str("x+") + str(fit.intercept)
        plt.plot(xdata,f(xdata,fit.slope,fit.intercept),"-.",label="Regression")
    
       
        
    elif mode == "0 intercept":
        f = lambda x,a: a*x
        fit = cfit(f,xdata,ydata)
        plt.plot(xdata,f(xdata,*fit[0]),"-.",label="Regression")
        
    elif mode == "general function":
        if f==None:
            raise TypeError

        fit = cfit(f,xdata,ydata)
        plt.plot(xdata,f(xdata,*fit[0]),"-.",label="Regression")
    
    
    else:
        print("mode='",mode,"' is not definted!")
        raise TypeError


    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.grid()
    if mode != "none":
        plt.legend()
        
    if show:
        plt.show()
    
    return (fig,fit)

def Rsqrue(x,y):
    RSS = np.sum((y-x)**2)
    TSS = np.sum(y**2)
    return 1 - RSS/TSS

def Reg_print(fit):
    m = ufloat(fit.slope,fit.stderr*2)
    b = ufloat(fit.intercept,fit.intercept_stderr*2)
    print("==> y =(",m,")x + (",b,") , R^2=",fit.rvalue**2)

#part 0
V_max = 565e-3 #volt
V_min = 395e-3 #volt
V = np.mean([V_max,V_min])
I = 0.5804 #A
r = V/I
print(f"r={r}")
#%% part 1
r = 0.82 #ohm
r_err = 0.05*r

min_C=0
min_f = 96.6e6 #Hz
min_f_err = 0

max_C = 60
max_f = 95.5e6 #Hz
max_f_err = 0

C = 55.5
C_err = 0.5

v_RF = (max_f - min_f)/(max_C - min_C) * (C - min_C) + min_f

h = 6.624e-34 # J*s
g = 2.0023
mu_B = 0.927e-23 # J/T
mu0 = 1.25663706212e-6
B_res = h * v_RF / (g * mu_B)

print(f"v_RF={v_RF/1e6}MHz\nB_res={B_res*1e3}mT")

#%% part 1
V_min = 460e-3 #volt (negativ)
V_max = 440e-3 #volt

V = np.mean([V_min,V_max])

V = 432.623e-3 # volt
I = V/r
H = B_res / mu0

k1 = H / I

print(f"k_1={k1}")


#%% part 2

I = 0.4766 #A
I_second = 0.5120#A
I_third = 0.5241#A

k2 = H / I

print(f"k_2={k2}")

R_data = pd.read_csv("part2 0.csv",header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
t = R_data["Time (s)"]
x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]

plt.plot(x,y)
plt.grid()


#%% part 3

I = 0.5513 #A

k3 = H / I

print(f"k_3={k3}")


R_data = pd.read_csv("part3 0.csv",header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
t = R_data["Time (s)"]
x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]

plt.plot(x,y)
plt.grid()

#%% k theory

