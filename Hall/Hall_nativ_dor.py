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

rho0 = R0*d*W/L

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
RH1 = slope * d / B
print("RH=",RH,"m^-3*C-1")

if RH1>0:
    print("The s/c is of P-type!")
else:
    print("The s/c is of N-type!")

n = 1/np.abs(RH1)
print("n=",n,"m^-3")

mu = np.abs(RH1)/rho0
print("mu=",mu,"m^-4*(ohm*C)-1")


#%% part 2
I = 30e-3 #A
I_err = 0
I = ufloat(I,I_err)

B = np.array([    -300, -280, -260, -240, -220, -200, -180, -160, -140, -120, 
                  -100,  -80,  -60,  -40,  -20,    0,   20,   40,   60,   80,
                   100,  120,  140,  160,  180,  200,  220,  240,  260,  280,
                   300]) * 1e-3 #T
B_err = np.array([ 0]) #T
UH = np.array([    0]) #V
U_err = np.array([ 0]) #V


fig,fit = one4all(xdata=B,ydata=UH,xerr=B_err,yerr=U_err,xlabel=r"$B[T]$",ylabel=r"$U_H[V]$",mode="linear")
Reg_print(fit)

slope = ufloat(fit.slope,fit.stderr*2)
RH2 = slope * d / I
print("RH=",RH2,"m^-3*C-1")

if RH2>0:
    print("The s/c is of P-type!")
else:
    print("The s/c is of N-type!")

n = 1/np.abs(RH2)
print("n=",n,"m^-3")

mu = np.abs(RH2)/rho0
print("mu=",mu,"m^-4*(ohm*C)-1")

rel_dif = np.abs(RH1-RH2)/max(np.abs(RH1.n),np.abs(RH2.n)) * 100
print("The differnce between RH1 to RH2 is:",rel_dif,"%")

#%% part 3
Ip = 30e-3 # A
Ip_err = 0
Ip = ufloat(Ip,Ip_err)

# Take only the higest precesion
B = np.array([     0]) * 1e-3 #T
B_err = np.array([ 0]) #T
Up = np.array([    0]) #V
Up_err = np.array([ 0]) #V

B = uarray(B,B_err)
Up = uarray(Up,Up_err)

R = Up/Ip

fig,fit = one4all(xdata=uval(B),ydata=uval(R),xerr=uerr(B_err),yerr=uerr(R),xlabel=r"$B[T]$",ylabel=r"$R[\Omega]$",mode="none")


#%% part 4
I = 30e-3 #A
I_err = 0
I = ufloat(I,I_err)

T = np.array([     0]) #C
T_err = np.array([ 0]) #C
Up = np.array([    0]) #V
Up_err = np.array([ 0]) #V

T = uarray(T,T_err)
Up = uarray(Up,Up_err)

fig,fit = one4all(xdata=uval(1/T),ydata=uval(1/Up),xerr=uerr(1/T_err),yerr=uerr(1/Up),xlabel=r"$T^{-1}[K^{-1}]$",ylabel=r"$U_p^{-1}[V^{-1}]$",mode="none")

fig,fit = one4all(xdata=uval(1/T),ydata=uval(np.log(Up)),xerr=uerr(1/T_err),yerr=uerr(np.log(Up),xlabel=r"$T^{-1}[K^{-1}]$",ylabel=r"$\ln(U_p)$",mode="linear")

#%% part 5
B= 250e-3 #T
B_err =0 
B = ufloat(B,B_err)