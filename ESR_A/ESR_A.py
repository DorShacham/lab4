#%% init

import dis
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

#%%
#part 0
V_max = ufloat(565e-3,1e-3)#volt
V_min = ufloat(395e-3,1e-3) #volt
V = np.mean([V_max,V_min])
I = ufloat(0.5804,1e-4) #A
r = V/I
print(f"r={r}")
#%% part 1
# r = 0.82 #ohm
# r_err = 0.05*r
# r = ufloat(r,r_err)

min_C=ufloat(0,0.5)
min_f = ufloat(96.6e6,0.1e6) #Hz

max_C = ufloat(60,0.5)
max_f = ufloat(95.5e6,0.1e6) #Hz

C = ufloat(55.5,0.5)

v_RF = (max_f - min_f)/(max_C - min_C) * (C - min_C) + min_f


from scipy.constants import h, physical_constants
h = h # Planck constant in J*s
g = 2.0023
mu_B = physical_constants['Bohr magneton'][0] # Bohr magneton in J/T
mu0 = physical_constants['vacuum mag. permeability'][0] # vacuum magnetic permeability
B_res = h * v_RF / (g * mu_B)

# h = 6.624e-34 # J*s
# g = 2.0023
# mu_B = 0.927e-23 # J/T
# mu0 = 1.25663706212e-6
# B_res = h * v_RF / (g * mu_B)


print(f"v_RF={v_RF/1e6}MHz\nB_res={B_res*1e3}mT")

#%% part 1
V_min = ufloat(460e-3,1e-3) #volt (negativ)
V_max = ufloat(440e-3,1e-3) #volt

V = np.mean([V_min,V_max])

V = ufloat(432.623e-3,0.001e-3) # volt
I = V/r
H = B_res / mu0

k1 = H / I

print(f"k_1={k1}")


#%% part 2

I_first = 0.4766 #A
I_second = 0.5120#A
I_third = 0.5241#A
I = np.array([I_first,I_second,I_third])
I_err = 0.0001
I = unumpy.uarray(I,I_err)
k2 = H / I

print(f"k_2={k2}")

R_data = pd.read_csv("part2 0.csv",header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
t = R_data["Time (s)"]
x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]

plt.plot(x,y)
plt.grid()


#%% part 3

I = ufloat(0.5513,0.0001) #A

k3 = H / I

print(f"k_3={k3}")


R_data = pd.read_csv("part3 0.csv",header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
t = R_data["Time (s)"]
x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]

plt.plot(x,y)
plt.grid()

#%% k theory


#%% ESR_B

# 5)
amplitude_of_wave_gen_sin= 0 #mV

# 7)

amp_X= 0 #mV
amp_Y= 0 #mV
phase_diff_Y_X= 0  # or pi?


avg_X=0 #mV
I0=avg_X/r #mA
'''
## alternative way to calculate avg_X (?)
max_X=0 #mV
min_X=0 #mV
alt_avg_X=(max_X+min_X)/2
'''


# 8) measure amplitude of absorption and phase difference for different I0 values (calculated from Avg_X)

avg_Vx_arr = np.array([]) #
amp_Vy_arr=np.array([])
phase_diff_arr=np.array([])

I0_arr=avg_Vx_arr/r

diff_I0= np.diff(I0_arr)
diff_Y=np.diff(amp_Vy_arr)

discrete_derivative_absorption_by_I0= np.divide(diff_Y,diff_I0)
plt.figure(dpi=300)
plt.plot(discrete_derivative_absorption_by_I0,I0_arr[:-1])
plt.grid(True)
plt.title('derivative of absorption by I0')
plt.xlabel('I0')
plt.ylabel('dA/dI0')


# 9
absorption_reconstructed=scipy.integrate.cumtrapz(discrete_derivative_absorption_by_I0, I0, dx=1.0, axis=-1, initial=None)

#10
plt.figure(dpi=300)
plt.plot(absorption_reconstructed,I0_arr[:-1])
plt.grid(True)
plt.title('absorption reconstructed by I0')
plt.xlabel('I0')
plt.ylabel('Abs_reconstructed')


#11



'''
R_data = pd.read_csv("ESR_B_0.csv",header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
t = R_data["Time (s)"]
x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]
'''


def calculate_theoretical_k(d,h,L,N):
    cosa=(h/2)/np.sqrt((h/2)**2+(d/2)**2)
    k=cosa*N/h
    return k

d_inner=38 #mm
d_outer=43 #mm
d_linspace=np.linspace(d_inner,d_outer,4)
L=70 #mm
h=L
N=440

k_arr=[calculate_theoretical_k(d,h,L,N/d_linspace.size) for d in d_linspace]

print(f'k for each layer {k_arr}')
print(f'k sum {np.sum(k_arr)}')

k_inner=calculate_theoretical_k(d_inner,h,L,N)
print(f'k inner = {k_inner}')
k_outer=calculate_theoretical_k(d_outer,h,L,N)
print(f'k inner = {k_outer}')