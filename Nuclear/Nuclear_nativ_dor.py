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

                            ######### Station Number: 3 #########
#%% part - Pleateau
parameters={'step_voltage':15,'preset_time':2} # [V],[sec]
source_chosen="SR-90"

file_name = "plat_itai.tsv"
df = pd.read_csv(file_name,sep='\t+', header=9) # read the data.
counts = np.array(df["Counts"]) # need to change to relevant name
V = np.array(df["Voltage"]) # need to change to relevant names

counts = uarray(counts,np.sqrt(counts))

np.random.seed(1)
noise = np.random.normal(loc=0,scale=uerr(counts))
counts = counts + uarray(noise,0)
cps = counts/parameters["preset_time"] 


V = uarray(V,15)
noise = np.random.normal(loc=0,scale=15)
V = V + noise
print(type(V[0]),type(cps[0]))
fig,fit = one4all(V,cps,xlabel="V[V]",ylabel="counts per second",show=False)
plt.figure(fig)
plt.savefig("fig/part1_1.png")
'''
Nativ:
is counts actually the cps? or should we divide by the elapsed time in seconds
from the video explaining the experimental setup it seems like the counts is always rising.
Alternative-
time_elapsed =np.array(df["time_elapsed"])# need to change to relevant name
cps=countes/time_elapsed
one4all(V,cps,xlabel="V[V]",ylabel="counts per second")
'''
operatin_V = 1000 # volt

#%% part - Statistics of counting




# 2. background meas
counts = ufloat(29,np.sqrt(29))
time = 100 # sec
time_err = 0
background_rate = counts/time 

# 3. get measurments
source_chosen="Co-60"

# 4.
source_to_counter_distance = '' # at which slot the source was put

file_name = "stat1.tsv"
df = pd.read_csv(file_name,sep='\t+', header=9) # read the data.
counts = np.array(df["Counts"]) # need to change to relevant name
time = 1 #sec
R = counts / time
n_bar = np.mean(R)

#5.
#m_prime=150*n_bar #number of measurments

file_name = "stat2.tsv"
df = pd.read_csv(file_name,sep='\t+', header=9) # read the data.
counts = np.array(df["Counts"]) # need to change to relevant names
R = counts/time

m = len(R)
n_bar = np.mean(R)
n_std = np.std(R)
n_bar_std = n_std/np.sqrt(m-1)
print(f"n_bar={n_bar}+-{n_bar_std}\nn_std={n_std}")

K3 = 1/(m-1) * np.sum((R-n_bar)**3)
K3_std = np.sqrt(np.var((R-n_bar)**3)/(m-1))

#def K3(data): n=np.mean(data) K3=np.sum((data-n)**3) K3=K3/(len(data)-1) return K3

print(f"K3={K3}+-{K3_std}")


plt.figure(dpi=300)
counts, bins = np.histogram(R) 
plt.stairs(counts, bins)

x = np.array(range(15))
pois = lambda x,lam: lam**x/np.math.factorial(int(x))*np.exp(-lam)
y = np.array([pois(xi,n_bar) for xi in x]) * m
plt.plot(x,y,"-.",label="Theory")
plt.grid()
plt.legend()
plt.xlabel("Number of counts in 1 second")
plt.ylabel("Number of occurrences")
plt.savefig("fig/part2_1.png")

#%% Inverse square law


counts = np.array([1034,1241,1249,1035,1049,1084,1094,1091,1206,1344])
time = np.array([69,75,60,40,32,24,19,12,8,5])

counts = uarray(counts,np.sqrt(counts))
time = uarray(time,1)
R = counts/time

# 1.
total_length = 101.5 *1e-3 #m
total_err = 0.5e-3 #m
one_level = total_length / 10
x = np.arange(total_length,0,-one_level)
x_err = total_err / 10

x = uarray(x,x_err)
m = len(R)

# 3.
beta_source_chosen= "Sr-90"

R_b = background_rate
Y = np.array([1/sqrt(r-R_b) for r in R])
fig,fit = one4all(x,Y,xlabel="x[m]",ylabel=r"$\frac{1}{\sqrt{R-R_b}}$",mode="linear",show=False)
Reg_print(fit)
plt.figure(fig)
plt.savefig("fig/part3_1.png")

m = ufloat(fit.slope,fit.stderr)
b = ufloat(fit.intercept,fit.intercept_stderr)
a = b/m
print(f"a={a}")



#%% part - Range of alpha particles
time_err = 1 
counts = np.array([6,5,27,54,71,109,140,386,238,211,296,290]) 
time = np.array([20,20,21,21,21,21,22,60,22,21,20,21]) #sec
counts = uarray(counts,np.sqrt(counts))
time = uarray(time,1)

R = counts/time
x = np.arange(20,8,-1) *1e-3
x_err = 1e-4
x = uarray(x,x_err)
#x = np.array([x[-2],x[-2]-1e-3,x[-2]-2e-3,x[-2]-3e-3,x[-2]-4e-3,x[-2]-5e-3,x[-2]-6e-3,x[-2]-7e-3,x[-2]-8e-3,x[-2]-10e-3,x[-2]-11e-3,x[-2]-12e-3]) # meter

R = R - background_rate
R = R * (x+a)**2
fig,fit = one4all(x+a,R,xlabel="range[m]",ylabel="rate[cps]",mode="linear")
Reg_print(fit)

m = ufloat(fit.slope,fit.stderr)
b = ufloat(fit.intercept,fit.intercept_stderr)
zero_rate = -b/m
print(f"zero_rate={zero_rate}")

#from graph we got 
E = ufloat(5.0625,0.3125) #MeV

#%% part - Absorption of Beta Particles and Beta Decay Energy


beta_source_chosen= "Sr-90"
thickness = np.array([0,40,80,160,240,320,400])*1e-6 #m
counts = np.array([1305,1124,1145,1055,1067,1267,1055])
time = np.array([14,13,15,15,17,22,20]) #sec
R = counts / time
one4all(thickness,R-background_rate,xlabel="thickness",ylabel="rate [cps]")


beta_source_chosen= "Ti-204"
thickness = np.array([0,40,80,160,240,320,400])*1e-6 #m
counts = np.array([])
time = np.array([]) #sec

#not sure what to do
one4all(thickness,np.log10(R-R_b),xlabel="thickness",ylabel="rate [cps]")