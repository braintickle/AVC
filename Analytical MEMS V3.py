# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 21:32:04 2016

@author: rahishah
"""

import numpy as np
import matplotlib.pyplot as plt

pts = 50
pts2 = 5    
Angle = np.linspace(0,8,num = pts)
Ang_rad = np.deg2rad(Angle)

#Dimensions of fingers
T_f = 25e-6        #Thickness of fingers
L_f = 130e-6       #Length of fingers
L_fo = 50e-6       #Offset from tip of stationary fingers to axis of rotation
N_f = 60           #Number of fingers
eta_o = 8.854e-12  #Permittivity of air
G_s = 5e-6   #Lateral spacing between moving and stationary fingers
k_sp = 4.02e-07
T_t = 30e-6        #Thickness of moving fingers
T_s = 20e-6        #Thickness of fixed fingers
a = 30e-6
beta =  np.linspace(12,16,num = pts2)*np.pi/180 #np.arctan(T_f/L_f)
ho = np.linspace(1,10,num = pts2)*1e-06#np.linspace(0,5,num = pts2)*1e-06        #half offset 

def Cap_AVC(theta, beta):
    if theta < beta - (T_f/L_f):
        
        Cap = (eta_o*N_f*0.5/G_s)*np.square(T_f - L_fo*(beta - theta))/(beta - theta)
        
    else:
        
        Cap = (eta_o*N_f*0.5/G_s)*2*(T_f*(L_f - L_fo) - 0.5*(L_f**2 - L_fo**2)*(beta - theta))    
    return Cap

Cap_SVC = lambda theta: (eta_o*N_f*0.5/G_s)*(L_f**2 - L_fo**2)*theta
    
##############################################Our Model###########################################

def Cap_Mod(theta, beta, ho):
    
    Ha = (0.5*T_t - ho)
    Hb = T_s - ((T_t*0.5) + ho) + a*theta
    Hl = L_fo + a*(theta**2 - 1)
    Hb2 = T_t + 2*ho - T_s    
    t_b = T_t - (L_f - a)*beta
    
    if theta < ((T_t + 2*ho - t_b - T_s)/(L_f)):
        
        Cap2 = (eta_o*N_f*0.25/G_s)*np.square((Ha*(1 - theta*(beta - theta))) - (Hl*(beta - theta)) + Hb)/(beta - theta)
    else:
        Cap2 = (eta_o*N_f*0.25/G_s)*(((Ha*(1 - theta*(beta - theta))) - (Hl*(beta - theta)) + Hb) + (L_f*theta + t_b - Hb2))*(L_f - L_fo)
    
    return Cap2
    

##############################Del C#############################################################    

def DelC_Mod(theta, beta, ho):
    
    Ha = (0.5*T_t - ho)
    Hb = T_s - ((T_t*0.5) + ho) + a*theta
    Hl = L_fo + a*(theta**2 - 1)   
    t_b = T_t - (L_f - a)*beta
    
    if theta < ((T_t + 2*ho - t_b - T_s)/(L_f)):
        DelC = (eta_o*N_f*0.25/G_s)*(np.square(((Ha*(1 - theta*(beta - theta))) - (Hl*(beta - theta)) + Hb)/(beta - theta)) + (2*(Ha*(2*theta - beta)+ 2*a*theta*(beta - theta) - a*(theta**2) + a + L_fo )*((Ha*(1 - theta*(beta - theta))) - (Hl*(beta - theta)) + Hb)*(1/(beta - theta))))    
    else:
        DelC = (eta_o*N_f*0.25/G_s)*(L_f - L_fo)*((a*(3*(theta**2) - 1 - 2*beta*theta)) - beta*Ha + 2*Ha*theta + L_fo + L_f + 3*a) 
    
    return DelC


def DelC_AVC(theta, beta):
    if theta < beta - (T_f/L_f):    
        DelC = (eta_o*N_f*0.5/G_s)*(T_f + (L_fo*(theta - beta)))*(T_f + (L_fo*(beta - theta)))/np.square(beta - theta)
    else:
        DelC = (eta_o*N_f*0.5/G_s)*(L_f**2 - L_fo**2)
    return DelC
    
DelC_SVC = lambda theta: (eta_o*N_f*0.5/G_s)*(L_f**2 - L_fo**2)

################################Voltage#########################################################
def Vol_Mod(theta, beta, ho):
    return np.sqrt((2*k_sp*theta)/(DelC_Mod(theta, beta, ho)))
def Vol_AVC(theta,beta):
    return  np.sqrt((2*k_sp*theta)/(DelC_AVC(theta, beta)))

Vol_SVC = lambda theta:  np.sqrt((2*k_sp*theta)/(DelC_SVC(theta))) 
#################################Pull-in##########################################################

#PI = lambda theta: DelC_Mod(theta) - (eta_o*N_f*0.5/G_s)*0.5*theta*2*(np.square(Ha*beta + Ha + Hb))/((beta - theta)**3)

##############################################Plot time!########################################################    


#
#PullIn = []
#for P_COUNT in range (0,pts):
#   PC= PI(Ang_rad[P_COUNT])
#   PullIn.append(PC)

Capacitance_AVC = np.zeros((pts2,pts))
Capacitance_Mod = np.zeros((pts2,pts))
Capacitance_SVC = []
DelCap_Mod = np.zeros((pts2,pts))
DelCap_AVC = np.zeros((pts2,pts))
DelCap_SVC = []
V_Mod = np.zeros((pts2,pts))
V_AVC = np.zeros((pts2,pts))
V_SVC = []
V_Mod1 = np.zeros((pts2,pts))
V_AVC1 = np.zeros((pts2,pts))
V_SVC1 = []


for ii in range (0,pts2):
    
    for i in range (0,pts):
        z= Cap_Mod(Ang_rad[i], beta[ii], 2.5e-06)
        Capacitance_Mod[ii,i] = (z*1e12)
        z1= Cap_AVC(Ang_rad[i], beta[ii])
        Capacitance_AVC[ii,i] = (z1*1e12)
    
for ic in range (0,pts):
   z2= Cap_SVC(Ang_rad[ic])
   Capacitance_SVC.append(z2*1e12)    

for jj in range (0,pts2):
    
    for j in range (0,pts):
        z3= DelC_Mod(Ang_rad[j], beta[jj], 2.5e-06)
        DelCap_Mod[jj,j] = (z3)
        z4= DelC_AVC(Ang_rad[j], beta[jj])
        DelCap_AVC[jj,j] = (z4)
        
for jc in range (0,pts):        
    z5 = DelC_SVC(Ang_rad[jc])
    DelCap_SVC.append((z5))
    
for kk in range (0,pts2):
    for k in range (0,pts):
        z7= Vol_Mod(Ang_rad[k], beta[0], ho[kk])
        V_Mod[kk,k] = (z7)
        z8= Vol_AVC(Ang_rad[k], beta[0])
        V_AVC[kk,k] = (z8)
    
for kc in range (0,pts):
    z9= Vol_SVC(Ang_rad[kc])
    V_SVC.append(z9)    
    
for kk1 in range (0,pts2):
    for k1 in range (0,pts):
        z10= Vol_Mod(Ang_rad[k1], beta[kk1], 2.5e-06)
        V_Mod1[kk1,k1] = (z10)
        z11= Vol_AVC(Ang_rad[k1], beta[kk1])
        V_AVC1[kk1,k1] = (z11)
    
for kc1 in range (0,pts):
    z12= Vol_SVC(Ang_rad[kc1])
    V_SVC1.append(z12)  
 
#for count in range (0,pts):
#    z9= Force_Mod(Ang_rad[count])
#    F_Mod.append(z9)
#    z10= Force_AVC(Ang_rad[count])
#    F_Mod.append(z10) 
#    z11= Force_SVC(Ang_rad[count])
#    F_Mod.append(z11)
plt.figure(figsize=(16, 10))
plt.subplot(2, 2, 1)
plt.plot(Angle,Capacitance_Mod[0,:],'b-^', label=r'$\theta_{i} = 12$')
plt.plot(Angle,Capacitance_Mod[1,:],'b-*', label=r'$\theta_{i} = 13$')
plt.plot(Angle,Capacitance_Mod[2,:],'b-.', label=r'$\theta_{i} = 14$')
plt.plot(Angle,Capacitance_AVC[1,:],'r-',lw = 5, label = r'AVC $(\theta_{i} = 12)$')
plt.plot(Angle,Capacitance_SVC,'k-', lw = 2,label = 'SVC')
plt.xlabel(r'Angle $[\theta]$')
plt.ylabel("Capacitace[pF]")
plt.legend(loc = 2, fontsize=14)



plt.subplot(2, 2, 2)
plt.plot(Angle,DelCap_Mod[0,:],'b-*', label=r'$\theta_{i} = 12$')
plt.plot(Angle,DelCap_Mod[1,:],'b-.', label=r'$\theta_{i} = 13$')
plt.plot(Angle,DelCap_Mod[2,:],'b-^', label=r'$\theta_{i} = 14$')
plt.plot(Angle,DelCap_AVC[1,:],'r-',lw = 5, label = r'AVC $(\theta_{i} = 12)$')
plt.plot(Angle,DelCap_SVC,'k-',lw = 2, label = 'SVC')
plt.xlabel(r'Angle $[\theta]$')
plt.ylabel( r'$\frac{\partial C}{\partial \theta}$', fontsize=25,rotation=0)
plt.legend(loc = 4, fontsize=14)

x1 = np.argmax(V_Mod[0,:])
x2 = np.argmax(V_AVC[1,:])
x3 = np.argmax(V_AVC[2,:])
#V_max = V_Mod[x]
#Ang_max = Ang_rad[x]
#Force_Max = k_sp*Ang_max
#Force_Max_FEA = 4.4e-08
#Force_Discrepency = Force_Max - Force_Max_FEA
plt.subplot(2, 2, 3)
plt.plot(V_Mod1[0,0:x1], Angle[0:x1],'b-*', label=r'$\theta_{i} = 12$')
plt.plot(V_Mod1[1,0:x2], Angle[0:x2],'b-.', label=r'$\theta_{i} = 13$')
plt.plot(V_Mod1[2,0:x3], Angle[0:x3],'b-^', label=r'$\theta_{i} = 14$')
plt.plot(V_AVC1[1,:], Angle[:],'r-',lw = 5, label = r'AVC $\theta_{i} = 12$')
plt.plot(V_SVC1[:], Angle[:],'k-',lw = 2, label = 'SVC')
plt.xlabel(r'Voltage [V]')
plt.ylabel(r'Angle $[\theta]$')
plt.legend()
plt.legend(loc = 2, fontsize=14) 

plt.subplot(2, 2, 4)
plt.plot(V_Mod[0,0:x1], Angle[0:x1],'b-*', label=r'$h{o} = 1\mu m$')
plt.plot(V_Mod[1,0:x2], Angle[0:x2],'b-.', label=r'$h{o} = 3\mu m$')
plt.plot(V_Mod[2,0:x3], Angle[0:x3],'b-^', label=r'$h{o} = 5\mu m$')
plt.plot(V_AVC[1,:], Angle[:],'r-',lw = 5, label = r'AVC $h{o} = 1\mu m$')
plt.plot(V_SVC[:], Angle[:],'k-',lw = 2, label = 'SVC')
plt.xlabel(r'Voltage [V]')
plt.ylabel(r'Angle $[\theta]$')
plt.legend()
plt.legend(loc = 2, fontsize=14)    


plt.tight_layout()
plt.show()