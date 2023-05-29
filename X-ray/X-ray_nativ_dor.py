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

####
from scipy.signal import find_peaks


# helper function for plotting data and regression
def one4all(xdata,ydata,yerr=0,xerr=0,mode="none",f=None,xlabel="x",ylabel="y",show=True,label="Data"):
    fig = plt.figure(dpi=300)
    if (np.sum(np.abs((uerr(xdata))))>0) or (np.sum(np.abs((uerr(ydata))))>0):
        plt.errorbar(uval(xdata),uval(ydata),uerr(ydata),uerr(xdata),".",label=label)
        xdata = uval(xdata)
        ydata = uval(ydata)
    else:
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


####
def smooth(y, box_pts):
     box = np.ones(box_pts)/box_pts #running average
     y_smooth = np.convolve(y, box, mode='same')
     return y_smooth

def plot_spec(file_name,f=None, height = 100, prominence=40):
    Mes1 = pd.read_csv('Mes1.txt',sep='\t',header=1) # read the data.
    Counts = np.array(Mes1['Impulses/#']) # Impulses
    Counts_smoothed=smooth(Counts, 10) # smooth the data over 10 channels
    Channels = np.array(Mes1['Channel/#']) # Channel
    
    peaks, _ = find_peaks(Counts_smoothed, height=height, prominence=prominence)

    xlabel = "Channels"
    if f is not None:
        Channels = f(np.array(Channels))
        xlabel = "Energy [eV]"

    plt.figure(dpi=300)
    plt.plot(Channels ,Counts_smoothed, label='Original Mo-tube spectrum')
    plt.plot(Channels[peaks] ,Counts_smoothed[peaks],"x",color='red',
    label='Lines')
    plt.ylabel('Impulses')
    plt.xlabel(xlabel)
    plt.legend()
#%% Information and configuration
                            ######### Station Number: _ #########

# congirusation we used to setup the experiment system:
# ____ : _____
# ____ : _____
# ____ : _____


#%%

channel = np.array([])
E = np.array([]) # eV

fig,fit = one4all(channel,E,xlabel="#Channel",ylabel="E[eV]",mode="linear")

f = lambda x: fit.slope*x+fit.intercept


plot_spec("some file.csv") # with channels 
plot_spec("some file.csv",f=f) # with energy
