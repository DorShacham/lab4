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

def plot_spec(file_name,f=None, height = 10, prominence=8):
    Mes1 = pd.read_csv(file_name,sep='\t',header=1) # read the data.
    Counts = np.array(Mes1['Impulses/#']) # Impulses
    Counts_smoothed=smooth(Counts, 10) # smooth the data over 10 channels
    Channels = np.array(Mes1['Channel/#']) # Channel
    
    peaks, _ = find_peaks(Counts_smoothed, height=height, prominence=prominence)

    xlabel = "Channels"
    if f is not None:
        Channels = f(np.array(Channels))/1e3
        xlabel = "Energy [KeV]"

    plt.figure(dpi=300)
    plt.plot(Channels ,Counts_smoothed, label='Original Mo-tube spectrum')
    plt.plot(Channels[peaks] ,Counts_smoothed[peaks],"x",color='red',
    label='Lines')
    plt.ylabel('Impulses')
    plt.xlabel(xlabel)
    plt.grid()
    plt.legend()
    
    if f is not None:
        return f(np.array(peaks))
    return peaks
    

#%% Information and configuration
                            ######### Station Number: 3 #########

# congirusation we used to setup the experiment system:
    
# The tume is Mo
# offset : 3
# gain : 2
# angle : 1.9 degree
# 

#%%
#42 Mo 17,479.34 17,374.3 19,608.3 
#28 Ni 7,478.15 7,460.89 8,264.66
#29 Cu 8,047.78 8,027.83 8,905.29
#30 Zn 8,638.86 8,615.78 9,572.0 
#26 Fe 6,403.84 6,390.84 7,057.98
#82 Pb 10,551.5 10,449.5 12,613.7 12,622.6 14,764.4 2,345.5 


channel = np.array([1973, 2240,763, 859,827,942,909,1019,635,715, 1136 ,1385 ,1646])
E = np.array([17400, 19608,7469.52, 8264.66,8037.805, 8905.29,8627.32,9572,6397.34,7057.98,10500.5,12618.15,14764.4]) # eV

fig,fit = one4all(channel,E,xlabel="#Channel",ylabel="E[eV]",mode="linear")

f = lambda x: fit.slope*x+fit.intercept



#plot_spec("pre_calibreation.txt", height = 10, prominence=9,f=f) # with energy



#%%
peak = plot_spec("pre_calibreation.txt", height = 10, prominence=9)
# I = 0.02 mA
peak = plot_spec("Ni.txt", height = 10, prominence=9)

peak = plot_spec("Cu.txt", height = 10, prominence=9)

peak = plot_spec("Zn.txt", height = 10, prominence=9)
# I = 0.03 mA
peak = plot_spec("Fe.txt", height = 10, prominence=10)

peak = plot_spec("Pb.txt", height = 10, prominence=11)

print(peak)


#%% alloys
# Number 11 - FeTiO3
peak = plot_spec("alloy11.txt", height = 10, prominence=11,f=f)

print(peak)

#%% alloys
# Number 23 - Sb2S3
peak = plot_spec("alloy23.txt", height = 10, prominence=11,f=f)

print(peak)

#%% a measure with only the rubber and the holding
# Number 23 - Sb2S3
peak = plot_spec("background.txt", height = 10, prominence=11,f=f)

print(peak)