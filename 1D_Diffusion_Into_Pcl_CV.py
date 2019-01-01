#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 11:23:31 2018

Simulation of the equations from Zhang2000


   y=0     _____________________________ y=1            electrolyte
pcl center                           pcl surface       


For 1-D linear diffusion equation (dm/dt = -D d2m/dx2) the stability criteria is dt < dx^2/D

@author: damon
"""
import numpy
import matplotlib.pyplot as plt

###INSERT ONE OF THESE EQNS INTO LINE 43 ############
#Li-ion intercalation plateaus from Zhang2000
#The equation from Dolye1996 and Zhang2000 is typo-fixed error-fixed:  
#Ueq = 4.19829 + 0.0565661*numpy.tanh(-14.5546*y[-1] + 8.60942) - 0.0275479 * ((0.998432 - y[-1])**-0.492465 - 1.90111) - 0.157123*numpy.exp(-0.04738*y[-1]**8) + 0.810239 * numpy.exp(-40*y[-1] + 5.355)
#One intercalation plateau (see Excel file)
#Ueq=25*numpy.exp(-((1-y1)+0.02)/0.04)+2-35*numpy.exp(((1-y1)-1.03)/0.04)-(1-y1)*0.05
#Three intercalation plateaus (see Excel file)
#Ueq=25*numpy.exp(-((1-y1)+0.04)/0.02)-0.025*(1+numpy.erf(((1-y1)-0.66)/(0.03*numpy.sqrt(2))))-0.025*(1+numpy.erf(((1-y1)-0.34)/(0.03*numpy.sqrt(2))))+2-25*numpy.exp(((1-y1)-1.03)/0.02)-(1-y1)*0.03
#Sloping Three Peaks
#Ueq=25*numpy.exp(-((1-y1)+0.04)/0.02)-0.015*(1+numpy.erf(((1-y1)-0.66)/(0.03*numpy.sqrt(2))))-0.015*(1+numpy.erf(((1-y1)-0.34)/(0.03*numpy.sqrt(2))))+2.3-25*numpy.exp(((1-y1)-1.03)/0.02)-(1-y1)*0.5
#Sloping
#Ueq=25*numpy.exp(-((1-y1)+0.04)/0.02)+2.3-25*numpy.exp(((1-y1)-1.03)/0.02)-(1-y1)*0.5
#Experimental curves of Ueq vs y[-1]
#Ueq=numpy.interp(y[-1],imported_y,imported_Ueq)
###INSERT ONE OF THESE EQNS INTO LINE 43 ############
T =    numpy.float64(298) # Kelvin  
F =    numpy.float64(96500) # Faradays Constant
R =    numpy.float64(8.3) #(J/mole/Celcius)
beta = numpy.float64(0.5) #Butler-Volmer transfer coefficient
cl =   numpy.float64(1.0) #(mol / L)  Li+ concentration in liquid phase  (constant!)
r0 =   numpy.float64(15E-4) #(cm)  radius of spherical particle
c0 =   numpy.float64(2.28) #(mol / L) initial concentration of intercalated-Li throughout pcl
ct =   numpy.float64(24) #(mol / L) total concentration of intercalation sites (occupied or unoccpied)
Utop=  numpy.float64(2.3) #(mV /s)  CV scan rate
Ubottom=numpy.float64(1.5) #(mV /s)  CV scan rate
Ueq=lambda y1: 25*numpy.exp(-((1-y1)+0.04)/0.02)+2.3-25*numpy.exp(((1-y1)-1.03)/0.02)-(1-y1)*0.5
z=     numpy.float64(-1.0) #charge of the intercalating ion
kb =   numpy.float64(0.05) ###(cm^5/2 s^-1 mol^1/2)
D =    numpy.float64(1.0E-8) ###(cm^2 /s)
dtau=  numpy.float64(0.0001) #dtau=dt*D/r0/r0 ## 
tau_nodes=numpy.arange(0.,4000.0)*dtau####
saved_time_spacing = numpy.float64(30)####
vreal= numpy.float64(0.01)###(mV /s) scan rate
v=     vreal*(r0*r0/D)

#Create arrays to hold data
#dx_real = numpy.float64(1E-4)      # in units of cm.  dx is 100 nm
dx=       numpy.float64(0.05) #  non-dimensional x      
x_nodes=  numpy.arange(0,21)*dx  # 0 to 0.1 mm
#dt=       numpy.float64(0.001)
dtau=     dtau      #It was defined above
tau_nodes=tau_nodes #It was defined above
print('dt is ', dtau*r0*r0/D, ' seconds')
print('D dt/dx^2 is ',dtau/dx/dx, ' and must be less than 0.5')
y=        numpy.ones((x_nodes.size))*c0/ct  # array of y , y=c/ct
Ustart=Ueq(y[-1])# the initial applied potential
dydx=     numpy.zeros((x_nodes.size))  # array of y , y=c/ct
d2ydx2=   numpy.zeros((x_nodes.size))  # array of y , y=c/ct
dydtau=     numpy.zeros((x_nodes.size))
i_array=numpy.zeros(tau_nodes.size)


#STABILITY dt < dx^2/D  , in my case here dt = 0.001  and  dx^2=0.00001 and D=1


j=0
cs=y[-1]*ct
Ucollector=Ustart
eta = Ucollector-Ueq(y[-1])
dydx[-1]=r0*z*kb*(ct-cs)/cl/ct/D*(numpy.exp(-beta*R*T/F*eta) - numpy.exp((1-beta)*R*T/F*eta) )
dydx[0]=0
d2ydx2[0]=(2*y[1] - 2*y[0]) / dx / dx
d2ydx2[-1]=(dydx[-1]-dydx[-2])/dx
dydtau = d2ydx2 #+ 2/x_nodes*dydx
#print(j, Ucollector,eta,cs,dydx[-1],d2ydx2[-1])

#Create the arrays that will save the chronological results to disk, and viewed later for insights
saved_time_spacing = saved_time_spacing  #It was defined above
tau_array_saved = numpy.concatenate([[0.,1.,2.,4.],numpy.arange(6,tau_nodes.size,saved_time_spacing)])
y_saved=numpy.zeros((tau_array_saved.size,x_nodes.size))
i_saved=numpy.zeros((tau_array_saved.size))
Ucollector_saved=numpy.zeros((tau_array_saved.size))
saved_row_index=0
y_saved[saved_row_index,:]=y
i_saved[saved_row_index]=-F*z*D*dydx[-1]*ct/r0
Ucollector_saved[saved_row_index]=Ucollector

for j in range(1,tau_nodes.size):   #loop over time
    Ucollector=Ucollector+v*dtau
    if Ucollector >= Utop and numpy.sign(v)==1 : 
        v=-v
        Ucollector=Ucollector+v*dtau*2
    if Ucollector <= Ubottom and numpy.sign(v)==-1: 
        v=-v
        Ucollector=Ucollector+v*dtau*2
    y[:]=y[:]+dydtau*dtau                #update the concentration field inside the pcl 
    cs=y[-1]*ct                          #calculate the surface concentraction of intercalated ions
    eta = Ucollector-Ueq(y[-1])
    for i in range(1,x_nodes.size-1): dydx[i] = (y[i+1] - y[i-1]) / 2 / dx
    #I used this following equation for deriving the dydx[-1] equation seen below: Ueq=4.0 + R*T/F*numpy.log((ct/cs-1)/cl)
    dydx[-1]=r0*z*kb*(ct-cs)/cl/ct/D*(numpy.exp(-beta*R*T/F*eta) - numpy.exp((1-beta)*R*T/F*eta) )
    #print(Ucollector,dydx[-1])
    dydx[0]=0
    for i in range(1,x_nodes.size-1): d2ydx2[i] = (y[i-1] - 2*y[i] + y[i+1]) / dx / dx
    d2ydx2[-1]=(dydx[-1]-dydx[-2])/dx
    d2ydx2[0]=(2*y[1] - 2*y[0]) / dx / dx
    dydtau = d2ydx2 #+ 2/x_nodes*dydx
    
    
    if numpy.any(j==tau_array_saved):
        saved_row_index=saved_row_index+1
        print(Ucollector,eta,y[0],y[-3],y[-1])
        #print(j, Ucollector,eta,cs,dydx[-1],d2ydx2[-1])
        y_saved[saved_row_index,:]=y
        i_saved[saved_row_index]=-F*z*D*dydx[-1]*ct/r0
        Ucollector_saved[saved_row_index]=Ucollector
        
        
    
print('elapsed time is ',r0*r0/D*tau_nodes[-1], ' seconds')
plt.plot(Ucollector_saved,i_saved)

