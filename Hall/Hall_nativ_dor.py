#%% init

import dis
import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from scipy.optimize import curve_fit as cfit

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
    m = ufloat(fit.slope,fit.stderr*2)
    b = ufloat(fit.intercept,fit.intercept_stderr*2)
    print("==> y =(",m,")x + (",b,") , R^2=",fit.rvalue**2)

#%% part 0
                            ######### Station Number: __ #########

d = 1e-3 #m
L = 16e-3 #m
W= 10e-3 #m
e = 1.6e-19 #C

Ip = np.array([    -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30]) * 1e-3 #A
I_err = np.array([ 0]) #A
Up = np.array([    0]) #V
U_err = np.array([ 0]) #V

fig,fit = one4all(xdata=Ip,ydata=Up,xerr=I_err,yerr=U_err,xlabel=r"$I_p[A]$",ylabel=r"$U_p[V]$",mode="linear")
Reg_print(fit)

R0 = ufloat(fit.slope,fit.stderr*2)
print("r=",R0,"ohm")

#%% part 1

B = 250e-3 #T
B_err =0
B = ufloat(B,B_err)


Ip = np.array([    -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30]) * 1e-3 #A
I_err = np.array([ 0]) #A
UH = np.array([    0]) #V
U_err = np.array([ 0]) #V

fig,fit = one4all(xdata=Ip,ydata=UH,xerr=I_err,yerr=U_err,xlabel=r"$I_p[A]$",ylabel=r"$U_H[V]$",mode="linear")
Reg_print(fit)

slope = ufloat(fit.slope,fit.stderr*2)
RH = slope * d / B
print("RH=",RH,"m^-3*C-1")

if RH>0:
    print("The s/c is of P-type!")
else:
    print("The s/c is of N-type")

