#%% init

import dis
import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from scipy.optimize import curve_fit as cfit
from scipy.constants import physical_constants

from uncertainties.unumpy import uarray
from uncertainties import ufloat
from uncertainties.umath import *
from uncertainties.unumpy import nominal_values as uval
from uncertainties.unumpy import std_devs as uerr


# helper function for plotting data and regression
def one4all(xdata,ydata,yerr=0,xerr=0,mode="none",f=None,xlabel="x",ylabel="y",show=True,label="Data"):
    fig = plt.figure(dpi=300)
    plt.errorbar(xdata,ydata,yerr,xerr,".",label=label)

    
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
        plt.plot(xdata,f(xdata,*fit[0]),"-.",label="fit")
    
    
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
    m = ufloat(fit.slope,fit.stderr)
    b = ufloat(fit.intercept,fit.intercept_stderr)
    print("==> y =(",m,")x + (",b,") , R^2=",fit.rvalue**2)

                            ######### Station Number: 3 #########
#%% part - Pleateau
parameters={'step_voltage':0,'preset_time':0} # [V],[sec]
source_chosen=" "

file_name = ""
df = pd.read_csv(file_name) # read the data.
counts = np.array(df["counts"]) # need to change to relevant name
V = np.array(df["voltage"]) # need to change to relevant names
one4all(V,counts,xlabel="V[V]",ylabel="counts per second")

'''
Nativ:
is counts actually the cps? or should we divide by the elapsed time in seconds
from the video explaining the experimental setup it seems like the counts is always rising.
Alternative-
time_elapsed =np.array(df["time_elapsed"])# need to change to relevant name
cps=countes/time_elapsed
one4all(V,cps,xlabel="V[V]",ylabel="counts per second")
'''
operatin_V = 0 # volt

#%% part - Statistics of counting

# 2. background meas
counts = 0
time = 100 # sec
time_err = 0
background_rate = counts/time 

# 3. get measurments
source_name = ''

# 4.
source_to_counter_distance = '' # at which slot the source was put

time = 1 #sec
rates = np.array([]) / time
n_bar = np.mean(rates)

#5.
#m_prime=150*n_bar #number of measurments

file_name = ""
df = pd.read_csv(file_name) # read the data.
counts = np.array(df["counts"]) # need to change to relevant names
rates = counts/time

m = len(rates)
n_bar = np.mean(rates)
n_std = np.std(rates)
n_bar_std = n_std/np.sqrt(m-1)
print(f"n_bar={n_bar}+-{n_bar_std}\nn_std={n_std}")

K3 = 1/(m-1) * np.sum((rates-n_bar)**3)
K3_std = np.sqrt(np.var((rates-n_bar)**3)/(m-1))
print(f"K3={K3}+-{K3_std}")

#%% Inverse square law


file_name = ""
df = pd.read_csv(file_name) # read the data.
counts = np.array(df["counts"]) # need to change to relevant names
R = counts/time

# 1.
x = np.array([]) #m
x_err = 0
m = len(R)

# 3.
beta_source_chosen= ""


#R_b = np.random.normal(loc=background_rate,sclae=?,size=m)
R_b = np.random.poisson(lam=n_bar,size=m)
Y = 1/np.sqrt(R-R_b)
fig,fit = one4all(x,Y,xlabel="x[m]",ylabel=r"$frac{1}{\sqrt(R-R_b)}$",mode="regression")
Reg_print(fit)

a = fit.intercept

fig,fit = one4all(x[x<0.02],Y[x<0.02],xlabel="x[m]",ylabel=r"$frac{1}{\sqrt(R-R_b)}$",mode="regression")
Reg_print(fit)

a = fit.intercept

#%% part - Range of alpha particles
time = 0 # sec
time_err = 0 
counts = np.array([]) 
count_err = np.array([])
R = counts/time
x = np.array([]) # meter

R = R - R_b
R = R * (x+a)**2
fig,fit = one4all(x+a,R,xlabel="range[m]",ylabel="rate[cps]",mode="linear")


#%% part - Absorption of Beta Particles and Beta Decay Energy
# background
time = 0
time_err =0
counts = 0
R_b = counts/time

thickness = np.array([])
counts = np.array([])

time = 0
time_err = 0
thick_err =0
R = counts / time

#not sure what to do
one4all(thickness,R-R_b,xlabel="thickness",ylabel="rate [cps]")