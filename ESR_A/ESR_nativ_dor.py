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
min_C=ufloat(0,0.5)
min_f = ufloat(96.6e6,0.1e6) #Hz

max_C = ufloat(60,0.5)
max_f = ufloat(95.5e6,0.1e6) #Hz

C = ufloat(42.0,0.5)

v_RF = (max_f - min_f)/(max_C - min_C) * (C - min_C) + min_f
print(f'V_rf is {v_RF}')

#%%
# 5)
amplitude_of_wave_gen_sin= 655# mVpp

# where we search for the absorbption
Vmax= 14.46
Vmin=9
V_change_of_phase=13.13

# we chose an averaging for Y to be
avg=64

# 7)

amp_X= 2.87 #mV
amp_Y= 0 #mV
phase_diff_Y_X= 0  # or pi?


avg_X=490 #mV
I0=avg_X/r #mA
'''
## alternative way to calculate avg_X (?)
max_X=0 #mV
min_X=0 #mV
alt_avg_X=(max_X+min_X)/2
'''
print(f'I0 = {I0}')
#%%
# 8) measure amplitude of absorption and phase difference for different I0 values (calculated from Avg_X)

pp_Vx_arr_old = np.array([1.8,1.9,1.6,1.6,1.4, 1.5, 1.4, 1.5, 1.5, 1.4, 1.4,
                      1.4,1.4, 1.4, 1.4,1.5]) #mV
amp_Vy_arr_old=np.array([23, 18,26, 30, 58, 66, 64, 65, 47, 8,28,
                     65, 62, 53, 42,29 ]) #mV
V_avg_old=np.array([487,488,476,474,475, 473,468, 471,467, 462,460
                ,456,452,448,445,441]) #mV
phase_diff_arr_old=(-1)**(np.array([0,0,0,0,0,0,0,0,0, 1,1,
                                1,1,1,1,1]))



pp_Vx_arr = np.array([2,1.5 ,1.7, 1.8, 1.9, 2, 2, 2,
                     1.4,1.4, 1.4, 1.4,1.5]) #mV
pp_Vx_arr_err=np.array([ 0.1, 0.1, 0.1, 0.1, 0.1,0.1,
                        0.1,0.1,0.1,0.1,0.1]) #mv




amp_Vy_arr=np.array([ 17,30, 51, 60, 62, 44,29,
                    
                     65, 62, 53, 42,29 ,11]) #mV
amp_Vy_arr_err=np.array([ 2, 2, 2, 2, 1,1,
                         1,2,2,2,1]) #mV

phase_diff_arr=(-1)**(np.array([0,0,0,0,0, 0,1,
                                1,1,1,1,1,1]))




V_avg=np.array([479, 477 , 475, 472 ,468, 465,462

                ,456,452,448,445,441, 435]) #mV
V_avg_err=np.array([1,1,1,1,1,1,1,1,1 ,1,1,1,1,1,1,1,1]) #mV

#%%
r=0.82
I0_arr=V_avg/r

'''
diff_I0= np.diff(I0_arr)
diff_Y=np.diff(amp_Vy_arr)
'''

#discrete_derivative_absorption_by_I0= np.divide(amp_Vy_arr,pp_Vx_arr)
discrete_derivative_absorption_by_I0=amp_Vy_arr/pp_Vx_arr*phase_diff_arr

plt.figure(dpi=300)
plt.plot(I0_arr,discrete_derivative_absorption_by_I0,".")
plt.grid(True)
plt.title('derivative of absorption by I0')
plt.xlabel('I0')
plt.ylabel('dA/dI0')

#%%
# 9
absorption_reconstructed=scipy.integrate.cumtrapz(discrete_derivative_absorption_by_I0, I0_arr, dx=1.0, axis=-1, initial=None)

#10
plt.figure(dpi=300)
plt.plot(I0_arr[:-1],absorption_reconstructed,".")
plt.grid(True)
plt.title('absorption reconstructed by I0')
plt.xlabel('I0')
plt.ylabel('Abs_reconstructed')


#11
# u as a function of B as a function of I
# The reason there weren't more points is because we couldn't controll the sensitivity of the X measurment in the scope without changing the measurment of AVG_1

'''
R_data = pd.read_csv("ESR_B_0.csv",header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
t = R_data["Time (s)"]
x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]
'''