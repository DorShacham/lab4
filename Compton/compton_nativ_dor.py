#%% init
import pandas as pd # read the data from files
import matplotlib.pyplot as plt # plot the data
import numpy as np # just for fun :)
from scipy.signal import find_peaks # find the first estimate to a line.
from scipy.stats import linregress # for calibration and compton fit
from scipy.optimize import curve_fit # for gausian fit
import scipy.constants as const # physical constants.


from uncertainties.unumpy import uarray
from uncertainties import ufloat
from uncertainties.umath import *
from uncertainties.unumpy import nominal_values as uval
from uncertainties.unumpy import std_devs as uerr

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


#Define the Gaussian function

def Gauss(x, H, A, x0, sigma):
 return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

# Mean over near chanels.
def smooth(y, box_pts):
 box = np.ones(box_pts)/box_pts
 y_smooth = np.convolve(y, box, mode='same')
 return y_smooth

# just linear function
def lin(x,a,b):
 return a*x+b

# convertion of energies in eV to wavelength in m
def e2lam(e):
 hc=1.24*10**-6; # eV*m
 return hc/e

# compton shift in whavelength
def comp(theo_angles):
 h=6.626*10**-34 # eV*s
 c=2.9*10**8; # m/s
 me=9.1*10**-31; # kg
 return h/(me*c)*(1-np.cos(theo_angles * np.pi/180))+e2lam(17479)

# The Klein-Nishina formula to energy 17479 eV
def klein_nish(A,theo_angles):
 E0=17479 # eV . The line energy without the Plexiglas
 lam0=e2lam(E0)
 return 1/2*A**2*(lam0/comp(theo_angles))**2*(lam0/comp(theo_angles)+comp(theo_angles)/lam0-np.sin(theo_angles * np.pi /180)**2)


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

def plot_spec(file_name,f=None, height = 10, prominence=8, name= "Mo-tube", background = None):
    Mes1 = pd.read_csv(file_name,sep='\t',header=1) # read the data.
    Counts = np.array(Mes1['Impulses/#']) # Impulses
    Counts_smoothed=smooth(Counts, 10) # smooth the data over 10 channels
    Channels = np.array(Mes1['Channel/#']) # Channel
    
    peaks, _ = find_peaks(Counts_smoothed, height=height, prominence=prominence)

    xlabel = "Channels"
    if f is not None:
        Channels = f(np.array(Channels))/1e3
        xlabel = "Energy [KeV]"

    if background is not None:
        new_peaks = []
        for p in peaks:
            if np.min(np.abs(f(p)-np.array(background))) > 50:
                new_peaks.append(p)
        final_peaks = np.array(new_peaks)
    else:
        final_peaks = peaks

    plt.figure(dpi=300)
    plt.plot(Channels ,Counts_smoothed, label=str(f'Original {name} spectrum'))
    plt.plot(Channels[final_peaks] ,Counts_smoothed[final_peaks],"x",color='red',
    label='Lines')
    plt.ylabel('Impulses')
    plt.xlabel(xlabel)
    plt.grid()
    plt.legend()
    
    if f is not None:
        return f(np.array(final_peaks))
    else:
        return final_peaks
    


#%% Preparation Questions
'''
Descloizite = (Pb,Zn)2VO4OH
#82 Pb 74,969.4 72,804.2 84,936 (10,551.5 10,449.5 12,613.7 12,622.6 14,764.4) 2,345.5 
#30 Zn (8,638.86 8,615.78 9,572.0) 1,011.7 1,011.7 1,034.7 
#23 V (4,952.20 4,944.64 5,427.29) 511.3 511.3 519.2 
#8 O 524.9 
#1 H None


#42 Mo (17,479.34 17,374.3 19,608.3) 2,293.16 2,289.85 2,394.81 2,518.3 2,623.5 


'''

#%% Information and configuration
                            ######### Station Number: 3 #########

# congirusation we used to setup the experiment system:
    
# The tube is __
# offset : _
# gain : _
# angle : ___ degree
# 

#%% finding the conversion between channel and energy
channel = np.array([]) ##
E = np.array([]) #ev

fig,fit = one4all(channel,E,xlabel="#Channel",ylabel="E[eV]",mode="linear")

conv = lambda x: fit.slope*x+fit.intercept
Reg_print(fit)
#%% The loop
# Define the line energy and amplitudes
comp_amp=[]
comp_eng=[]
# Define 14 colors
colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple',
'pink', 'brown', 'black', 'gray', 'cyan', 'magenta', 'lime', 'navy']
# The loop run on different angles
plt.figure(dpi=300)
for i in range(1,15):
 # import the data
 data=pd.read_csv('run{}'.format(i), sep='\t', header=1)
 chanals=np.array(data['Channel/#'])
 Impulses=np.array(data['Impulses/#'])
 Imp_smooth=smooth(Impulses,50)

 # cut relevant interval
 x=conv(chanals);
 y=Imp_smooth;
 indss = (x>16000) & (x<18500)
 x=x[indss];
 y=y[indss];

 # plot the relevant interval
 plt.plot(x,y,':',color=colors[i-1])

 # first estimate the line energy
 peaks, properties = find_peaks(y, prominence=5,width=20,distance=1000)
 plt.plot(x[peaks], y[peaks], "rx") # plot the estimation

 # Fit the line to gaussian. p0 is the initial guess
 parameters, covariance = curve_fit(Gauss, x, y,p0=[10,y[peaks].item(),x[peaks].item(), 500]);

 plt.plot(x,Gauss(x,parameters[0],parameters[1],parameters[2],parameters[3]),'--',color=colors[i-1])

 #acumulate the line energies and amplitudes
 comp_amp.append(Gauss(parameters[2],parameters[0],parameters[1],parameters[2],parameters[3]))
 comp_eng.append(parameters[2]) # eV

 #%% compton
comp_amp = np.array(comp_amp)
comp_eng = np.array(comp_eng)
comp_lam = e2lam(comp_eng)

E0=17479 # eV . The line energy without the Plexiglas
lam0=e2lam(E0)
delta_lam = comp_lam -  lam0

x = np.cos(comp_eng * np.pi/180)
fig,fit = one4all(x,delta_lam,mode="linear")
Reg_print(fit)
plt.figure(fig)
plt.plot(x,comp(comp_eng),label="Theory")
plt.legend()

m = ufloat(fit.slope,fit.stderr)
lam_e = - m

h=6.626*10**-34 # eV*s
c=2.9*10**8; # m/s
me=9.1*10**-31; # kg
lam_e_theory =  h/(me*c)

print(f"Lam_e = {lam_e}")
print(f"Lam_e theory = {lam_e_theory}")

#%% The Klein-Nishina formula

fig, fit = one4all(comp_eng,comp_amp)
plt.figure(fig)

re = 2.8179403262e-15 # m
plt.plot(comp_eng,klein_nish(re,comp_eng),label="Theory")
plt.legend()