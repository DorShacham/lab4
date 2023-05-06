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

#%% part 1 graphs
files = ["part1_1mod6.5 0","part1_2mod7 0","part1_3mod7.1 0","part1_4mod8 0","part1_5 0"]
for file in files:
    R_data = pd.read_csv(str(file+".csv"),header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
    t = R_data["Time (s)"]
    x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]
    Verr = 1e-4
    Terr = 1e-6
    fig,fit = one4all(t,y,Verr,Terr,xlabel="Time[sec]",ylabel="V[V]",label="~A",show=False)
    plt.figure(fig)
    plt.errorbar(t,x,Verr,Terr,".",label=r"$V_r$[A]")
    plt.legend()
    fig.savefig(str("figs/A_"+file+".png"))
    if file == "part1_3mod7.1 0":
        x3 = x
    if file == "part1_1mod6.5 0":
        x1 = x
    if file == "part1_5 0":
        x4 = x




#%% part 1
# V_min = ufloat(460e-3,1e-3) #volt (negativ)
# V_max = ufloat(440e-3,1e-3) #volt
for x in [x3,x1,x4]:
    V_min = ufloat(abs(min(x)),1e-4)
    V_max = ufloat(abs(max(x)),1e-4)

    V = np.mean([V_min,V_max])

    # V = ufloat(432.623e-3,0.001e-3) # volt
    I = V/r
    H = B_res / mu0

    k1 = H / I

    print(f"k_1={k1}")


V_min = ufloat(abs(min(x3)),abs(min(x3)-min(x4)))
V_max = ufloat(abs(max(x3)),abs(max(x3)-max(x4)))

V = np.mean([V_min,V_max])

# V = ufloat(432.623e-3,0.001e-3) # volt
I = V/r
H = B_res / mu0

k1 = H / I

print(f"k_1={k1}")



#%% part 2

I_first = 0.4766 #A
I_second = 0.5120#A
I_third = 0.5241#A
I = I_third
I_err = I_third-I_first
I = uarray(I,I_err)
k2 = H / I

print(f"k_2={k2}")

R_data = pd.read_csv("part2 0.csv",header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
t = R_data["Time (s)"]
x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]

fig,fit = one4all(t,y,Verr,Terr,xlabel="Time[sec]",ylabel="V[V]",label="~A",show=False)
plt.figure(fig)
plt.errorbar(t,x,Verr,Terr,".",label=r"$V_r$[A]")
plt.legend()
fig.savefig("figs/A_part2_ty")

plt.figure(dpi=300)
plt.plot(x,y,".")
plt.xlabel(r"$V_r$ [Volt]", fontsize=14)
plt.ylabel("~A [V]",fontsize=14)
plt.plot(2*[np.mean([min(x),max(x)])],[min(y),max(y)],"-.",color="orange")
plt.grid()
plt.savefig("figs/A_part2_xy")

#%% part 3

I = ufloat(0.5513,0.0001) #A

k3 = H / I

print(f"k_3={k3}")


R_data = pd.read_csv("part3 0.csv",header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
t = R_data["Time (s)"]
x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]

plt.figure(dpi=300)
plt.plot(x,y,".")
plt.xlabel(r"$V_r$ [Volt]", fontsize=14)
plt.ylabel("~A [V]",fontsize=14)
plt.plot(2*[np.mean([min(x),max(x)])],[min(y),max(y)],"-.",color="orange")
plt.grid()
plt.savefig("figs/A_part3")


#%% k theory
h = L = 70e-3
d1 = 43e-3
N = 440
cosa=(h/2)/np.sqrt((h/2)**2+(d1/2)**2)

k_theory = N/L*cosa
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
amplitude_of_wave_gen_sin= ufloat(655,1)# mVpp

# where we search for the absorbption
Vmax= ufloat(14.46,0.01)
Vmin=ufloat(9,0.01)
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
                      1.4,1.4, 1.4, 1.4,1.5])*1e-3 #V
amp_Vy_arr_old=np.array([23, 18,26, 30, 58, 66, 64, 65, 47, 8,28,
                     65, 62, 53, 42,29 ])*1e-3 #V
V_avg_old=np.array([487,488,476,474,475, 473,468, 471,467, 462,460
                ,456,452,448,445,441])*1e-3 #V
phase_diff_arr_old=(-1)**(np.array([0,0,0,0,0,0,0,0,0, 1,1,
                                1,1,1,1,1]))



pp_Vx_arr = np.array([2,1.5 ,1.7, 1.8, 1.9, 2, 2, 2,
                     1.4,1.4, 1.4, 1.4,1.5])*1e-3 #V
pp_Vx_arr_err=np.array([ 0.1, 0.1, 0.1, 0.1, 0.1,0.1,
                        0.1,0.1,0.1,0.1,0.1,0.1,0.1])*1e-3 #v




amp_Vy_arr=np.array([ 17,30, 51, 60, 62, 44,29,
                    
                     65, 62, 53, 42,29 ,11])*1e-3 #V
amp_Vy_arr_err=np.array([ 2, 2, 2, 2, 1,1,
                         1,2,2,2,1,1,1])*1e-3 #V

phase_diff_arr=-(-1)**(np.array([0,0,0,0,0, 0,1,
                                1,1,1,1,1,1]))




V_avg=np.array([479, 477 , 475, 472 ,468, 465,462

                ,456,452,448,445,441, 435])*1e-3 #V
V_avg_err=np.array([1,1,1,1,1,1,1,1,1 ,1,1,1,1])*1e-3 #V


V_avg = uarray(V_avg,V_avg_err)[::-1]
amp_Vy_arr = uarray(amp_Vy_arr,amp_Vy_arr_err)[::-1]
pp_Vx_arr = uarray(pp_Vx_arr,pp_Vx_arr_err)[::-1]

#%%

# r=0.82
I0_arr=V_avg/r
discrete_derivative_absorption_by_I0=-amp_Vy_arr/pp_Vx_arr*phase_diff_arr
# I0_arr,discrete_derivative_absorption_by_I0 = zip(*sorted(zip(I0_arr,discrete_derivative_absorption_by_I0)))
# I0_arr = np.array(I0_arr)
# discrete_derivative_absorption_by_I0 = np.array(discrete_derivative_absorption_by_I0)

# C0 =k3.n * mu0 * g * mu_B / (h / (2*np.pi))
# f = lambda I0, w, c, T2: c*C0*w*T2*((w*T2)**2-(C0*T2*I0)**2+1)/(1+T2**2*(w-C0*I0)**2)**2
fig,fit = one4all(uval(I0_arr),uval(discrete_derivative_absorption_by_I0),xerr=uerr(I0_arr),yerr=uerr(discrete_derivative_absorption_by_I0),xlabel=r"$I_0[A]$",ylabel=r"~$\frac{dA}{dI_0}$",mode="none")
fig.savefig("figs/B_deriv")
# plt.figure(dpi=300)
# plt.plot(I0_arr,discrete_derivative_absorption_by_I0,"o",label="measured")
# plt.xlabel(r"$I_0[A]$",fontsize=14)
# plt.ylabel(r"$~\frac{dA}{dI_0}$[mV]",fontsize=14)
# plt.grid()
# bounds=([0.56*C0,0,0],[0.57*C0,np.inf,1e-7])
# fit = cfit(f=f,xdata=I0_arr,ydata=discrete_derivative_absorption_by_I0,bounds=bounds)
# I0 = np.linspace(min(I0_arr),max(I0_arr),100)
# plt.plot(I0,f(I0,*fit[0]),"-.")
# plt.plot(I0,f(I0,*(6.13350773e+11,2.03368596e+01, 6.52595646e-11)),"-.")
# w,c,T2 = (6.13350773e+11,2.03368596e+01, 6.52595646e-11)
# plt.ylim(min(discrete_derivative_absorption_by_I0),max(discrete_derivative_absorption_by_I0))
# print(fit[0])

#%%
# 9
I0_arr = uval(I0_arr)
discrete_derivative_absorption_by_I0 = uval(discrete_derivative_absorption_by_I0)

absorption_reconstructed=scipy.integrate.cumtrapz(discrete_derivative_absorption_by_I0, I0_arr, dx=1.0, axis=-1, initial=0)

#10
fig,fit = one4all(I0_arr,absorption_reconstructed,xlabel=r"$I_0[A]$",ylabel=r"~$A[A]$")
fig.savefig("figs/B_A(I)")

#11
# u as a function of B as a function of I
# The reason there weren't more points is because we couldn't controll the sensitivity of the X measurment in the scope without changing the measurment of AVG_1

'''
R_data = pd.read_csv("ESR_B_0.csv",header=1, usecols=["Time (s)", "1 (VOLT)", "2 (VOLT)"]) #header = 1 means that the line 1 (the second line) is the header
t = R_data["Time (s)"]
x, y = R_data["1 (VOLT)"], R_data["2 (VOLT)"]
'''

#%% transforming to w,w0 and calculating T2
w_measured = v_RF * 2*np.pi
B = k3.n * mu0 * I0_arr 
w0 = B * g * mu_B / (h / (2*np.pi))

f = lambda w0,w,T2,C: w0*w*C*T2/(1+(w0-w)**2*T2**2)
plt.figure(dpi=300)
plt.plot(w0,absorption_reconstructed,"o",label="measured")
plt.xlabel(r"$\omega_0[rad/sec]$",fontsize=14)
plt.ylabel(r"$~A[A]$",fontsize=14)
plt.grid()
fit = cfit(f=f,xdata=w0,ydata=absorption_reconstructed,bounds=([6e8,0,0],[6.15e8,1e-7,1e-9]))

x = np.linspace(min(w0),max(w0),100)
plt.plot(x,f(x,*fit[0]),"-.",label="fit")
plt.legend()
plt.savefig("figs/B_A(w0)")

errs  = np.sqrt(np.diag(fit[1]))
w_fit = ufloat(fit[0][0],2*errs[0])
T2_fit = ufloat(fit[0][1],2*errs[1])

print(f"w form fit ={w_fit}\nw from measure = {w_measured}")
print(f"T2 form fit ={T2_fit}")
print(f"R^2={Rsqrue(f(w0,*fit[0]),absorption_reconstructed)}")
#%%
C0 =k3.n * mu0 * g * mu_B / (h / (2*np.pi))
f = lambda I0, w, c0, c, T2: 2*w*c0*c*T2**3*(w-c0*I0) / (1+(w-c0*I0)**2*T2**2)**2

#%%
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