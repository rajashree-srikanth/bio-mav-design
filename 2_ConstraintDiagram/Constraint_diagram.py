#!/usr/bin/env python
# coding: utf-8

# In[70]:


import numpy as np
import math
import matplotlib.pyplot as plt


####################SOME BASIC DATA - NOT IMP###############################
a = np.zeros(4)
g = np.zeros(4)
lo = 11.0168 #latitude data
h0 = 0.5 #kilometres
for i in range(0,4):
    a[i] = (lo + 9.04371713*(10**(-3)))*2 #for finding g value, accurate to altitude
print(a)
coslo = np.cos(a)
for i in range(0,4):
    g[i] = 9.80616 - (0.025928*coslo[i]) + (0.000069*((coslo[i])**2)) - [h0*3.086*(10**(-6))]
print(g) #final g, precise

Ta = 15 - (0.0065*h0) #temp of air
p = 1013*((1-(h0*2.26*(10**(-5))))**5.265) #pressure
rho = 1.226*(p/1013)*(288/(Ta+273)) #air density
nu = (1.466 + 0.09507*h0 + 0.010470*(h0**2))*(10**(-5)) #kinematic viscosity
print("Temperature,pressure,air density and kinematic viscosity are: ", Ta, p, rho, nu)
alp = ((Ta + 273.16)/(Ta + 273.16) - (0.001981*h0))*((1 - (((0.001981*h0)/288.16)))**5.256) #correction factor for SL thrust
print("Alpha : ",alp)

#tadarida australis:
# Enter take-off weight (in kg) 0.0353
# Enter wingspan (in m) 0.425
# Enter cruise velocity (in m/s) 12


# In[66]:


##########################MAIN DATA for TADARIA AUSTRALIS OBTAINED FROM SOURCES#################
wto = 0.0353
#wto = float(input('Enter take-off weight (in kg)')) #take-off wt, kg
e = 0.8 #oswald eff fact
b = 0.425
#b = float(input('Enter wingspan (in m)')) #wingspan, m
########################INPUT VELOCITIES###################################
U = 11.9
#U = float(input('Enter cruise velocity (in m/s)')) #cruise velocity, m/s
ustall = U/1.5 #stall velocity
U_climb = ustall*1.2 #climb vel, m/s
turn_speed = 8.4 #m/s 
#turn_speed = int(input('Enter turn velocity (in m/s)'))#turn vel, m/s
hor_acc = 1 #horizontal acc, m/s2
climb_acc = 0.6 #climb acc, m/s2
climb_rate = 1.2 #ROC

####################DYNAMIC VELOCITIES#############################
qturn = 0.5*rho*(turn_speed**2)
qcruise = 0.5*rho*(U**2) #cruise dyn pressure
qclimb = 0.5*rho*(U_climb**2) #climb dyn pressure


Re = 0.3*b*U/nu #reynolds number
k = 0.115
Clmax = 1.2
psi = 3.2 #ratio of the parasite drag coeff of flap wing (CDP) to the friction drag coeff (Cf) for a flat sheet
print("Reynolds number: ", Re)
n1 = 1.0 #load factor cruise
#n2 = 1.154 #load factor climb ###HOW?
n2 = math.cos(math.asin(climb_rate/U_climb))
print('Climb load factor:', n2)

##########drag-induced############

dind1 = ((wto*n1)**2)/(e*(3.1415)*(b**2)*qcruise) #induced drag at cruise
dind2 = ((wto*n2)**2)/(e*(3.1415)*(b**2)*qclimb)  #induced drag at climb
print("Induced Drag for load factor 1 : ",dind1)
print("Induced Drag for load factor",n2,'=', dind2)

cf = 0.455*((np.log10(Re))**(-2.58))
print("Cf : ",cf)

######################WING LOADING CALC###############################
dzbydt = ((climb_rate+(U*climb_acc/9.832))) #dZ/dt for acc climb

wbys = list(range(1,600))
TbyW1 = np.zeros(len(wbys))
TbyW2 = np.zeros(len(wbys))
TbyW3 = np.zeros(len(wbys))
TbyW4 = np.zeros(len(wbys))
TbyW5 = np.zeros(len(wbys))
wbys_stall = 0.5*rho*(ustall**2)*Clmax #stall condition
print("Stall loading = ",wbys_stall, "N/m^2")
#print("Stall loading = ",wbys_stall)
for i in range(0,599):
    TbyW1[i] = ((qcruise/wbys[i])*((k*(((n1*wbys[i])/qcruise)**2)) + (2*psi*cf) ))/alp #cruise
    TbyW2[i] = (((qclimb/wbys[i])*((k*(((n2*wbys[i])/qclimb)**2)) + (2*psi*cf))) + (climb_rate/U_climb))/alp #constant climb
    TbyW3[i] = (((qcruise/wbys[i])*((k*(((n1*wbys[i])/qcruise)**2)) + (2*psi*cf))) + (hor_acc*U/(U*9.81) ))/alp #horizontal acc
    TbyW4[i] = ((qclimb/wbys[i])*((k*(((n2*wbys[i])/qclimb)**2)) + (2*psi*cf)) + dzbydt/U_climb )/alp #acc climb
    TbyW5[i] = ((qturn/wbys[i])*((k*(((n1*wbys[i])/qturn)**2)) + (2*psi*cf)) )/alp #constant turn
     
#TbyW4_constrained = ((qclimb/wbys_stall)*((k*(((n2*wbys_stall)/qclimb)**2)) + (2*psi*cf)) + dzbydt/U )/alp
TbyW4_constrained = ((qclimb/wbys_stall)*((k*(((n2*wbys_stall)/qclimb)**2)) + (2*psi*cf)) + dzbydt/U_climb )/alp #obtaining the constraint pt
print("Minimum T/W = ", TbyW4_constrained)
print("load factor", math.cos(math.asin(climb_rate/U_climb)))
print("climb angle", 180/math.pi*(math.asin(climb_rate/U_climb)))
print("Mass of aircraft = ", 13*0.02584/9.81)


# In[67]:


###PLOTTING###
limits = [0, 100, 0, 1]                
plt.axis(limits)
plt.xlabel("W/S(N/m2)")
plt.ylabel("T/W") 
plt.axvline(x = wbys_stall, label = 'Hand launch - stall speed', color = 'red')    
plt.plot(wbys,TbyW1, label = 'Cruise')
plt.plot(wbys,TbyW2, label = 'Climb')
plt.plot(wbys,TbyW3, label = 'Horizontal acceleration')
plt.plot(wbys,TbyW4, label = 'Climb accelerated' , color = '#23395d')
plt.plot(wbys,TbyW5, label = 'Constant altitude and turn speed' , color = 'purple')
plt.plot(wbys_stall, TbyW4_constrained, "o",color = '#666666',markeredgecolor='black', label = 'minimum T/W')
plt.plot(13,0.425, 'ko', label='Design Point') #this point is obtained from section below
plt.legend()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) #loc of legend box
plt.show()


print("Mass of aircraft = ", wbys_stall*0.02584/9.81)


# In[68]:


print(180/math.pi*(math.asin(1)))
print(climb_rate)
print(U_climb)
print("load factor", math.cos(math.asin(climb_rate/U_climb)))
print("climb angle", 180/math.pi*(math.asin(climb_rate/U_climb)))


# In[69]:


#finding point of intersection - run code from beginning, since the values of wbys and TbyW4 changes due to code below
wbys = np.round(wbys,3)
#print(wbys)
TbyW3 = np.round(TbyW3,3)
TbyW4 = np.round(TbyW4,3)
#print(TbyW3)
for i in range(0,len(wbys)):
    if TbyW3[i]==TbyW4[i]:
        print ('ha')
        tbyw_dp = TbyW3[i]
        wbys_dp = wbys[i]
print(tbyw_dp)
print(wbys_dp)
        


# In[64]:


S = 0.02584 #FOR TADARIA AUSTRALIS
#W_dp = wbys_dp*S
power = tbyw_dp*U*wbys_dp*S
print("Power required = ", power, 'W')


# In[ ]:





# In[ ]:




