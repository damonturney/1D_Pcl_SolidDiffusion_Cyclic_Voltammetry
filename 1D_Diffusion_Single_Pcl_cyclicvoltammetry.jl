"""
Created on Fri Dec 14 11:23:31 2018

Simulation of the equations from Zhang2000

  current collector | ------------------- y=0 ------------------------ y=1|  electrolyte
                                       pcl center

For 1-D linear diffusion equation (dm/dt = -D d2m/dx2) the stability criteria is dt < dx^2/D

@author: damon
"""

function main()
    ####### EQUATIONS FOR FREE ENERGY OF THE ELECTRON IN THE ELECTROCHEMICAL REACTION ############
    #Li-ion intercalation plateaus from Zhang2000
    #The equation from Dolye1996 and Zhang2000 is typo-fixed error-fixed:
    #Ueq = 4.19829 + 0.0565661*tanh(-14.5546*y[-1] + 8.60942) - 0.0275479 * ((0.998432 - y[-1])**-0.492465 - 1.90111) - 0.157123*exp(-0.04738*y[-1]**8) + 0.810239 * exp(-40*y[-1] + 5.355)
    #One intercalation plateau (see Excel file)
    #Ueq=25*exp(-((1-y1)+0.02)/0.04)+2-35*exp(((1-y1)-1.03)/0.04)-(1-y1)*0.05
    #Three intercalation plateaus (see Excel file)
    #Ueq=25*exp(-((1-y1)+0.04)/0.02)-0.025*(1+erf(((1-y1)-0.66)/(0.03*sqrt(2))))-0.025*(1+erf(((1-y1)-0.34)/(0.03*sqrt(2))))+2-25*exp(((1-y1)-1.03)/0.02)-(1-y1)*0.03
    #Sloping Three Peaks
    #Ueq=25*exp(-((1-y1)+0.04)/0.02)-0.015*(1+erf(((1-y1)-0.66)/(0.03*sqrt(2))))-0.015*(1+erf(((1-y1)-0.34)/(0.03*sqrt(2))))+2.3-25*exp(((1-y1)-1.03)/0.02)-(1-y1)*0.5
    #Sloping
    #Ueq=25*exp(-((1-y1)+0.04)/0.02)+2.3-25*exp(((1-y1)-1.03)/0.02)-(1-y1)*0.5
    #Experimental curves of Ueq vs y[-1]
    #Ueq=interp(y[-1],imported_y,imported_Ueq)
    ######## INSERT ONE OF THESE EQNS INTO THE NEXT LINE ############
    Ueq(y1) = 25*exp(-((1-y1)+0.04)/0.02)+2.3-25*exp(((1-y1)-1.03)/0.02)-(1-y1)*0.5
    T =    298.0    # Kelvin
    F =    96500.0  # Faradays Constant
    R =    8.3      #(J/mole/Celcius)
    beta = 0.5      #Butler-Volmer transfer coefficient
    cl =   1.0      #(mol / L)  Li+ concentration in liquid phase  (constant!)
    r0 =   15E-4    #(cm)  radius of spherical particle
    c0 =   2.28     #(mol / L) initial concentration of intercalated-Li throughout pcl
    ct =   24       #(mol / L) total concentration of intercalation sites (occupied or unoccpied)
    Utop=  2.3      #(mV /s)  CV scan rate
    Ubottom=1.5     #(mV /s)  CV scan rate
    z=     -1.0     #charge of the intercalating ion
    kb =   0.05     ###(cm^5/2 s^-1 mol^1/2)
    D =    1.0E-8   ###(cm^2 /s)
    dtau=  0.0001   #dtau=dt*D/r0/r0 ##
    tau_nodes=Array(0:4000000.0)*dtau####
    saved_time_spacing = 30000####
    vreal= 0.00001###(mV /s) scan rate
    v=     vreal*(r0*r0/D)

    #Create arrays to hold data
    #dx_real = float64(1E-4)      # in units of cm.  dx is 100 nm
    dx=       0.05 #  non-dimensional x
    x_nodes=  Array(0:20)*dx  # 0 to 0.1 mm
    #dt=       float64(0.001)
    dtau=     dtau      #It was defined above
    tau_nodes=tau_nodes #It was defined above
    print("dt is ", dtau*r0*r0/D, " seconds\n")
    print("D dt/dx^2 is ",dtau/dx/dx, " and must be less than 0.5\n")
    y=        ones(size(x_nodes)[1])*c0/ct  # array of y , y=c/ct
    Ustart=Ueq(y[end])# the initial applied potential
    dydx=     zeros(size(x_nodes)[1])  # array of y , y=c/ct
    d2ydx2=   zeros(size(x_nodes)[1])  # array of y , y=c/ct
    dydtau=     zeros(size(x_nodes)[1])

    #STABILITY dt < dx^2/D  , in my case here dt = 0.001  and  dx^2=0.00001 and D=1
    j=1
    cs=y[end]*ct
    Ucollector=Ustart
    eta = Ucollector-Ueq(y[end]) #eta is overpotential
    dydx[end]=r0*z*kb*(ct-cs)/cl/ct/D*(exp(-beta*R*T/F*eta) - exp((1-beta)*R*T/F*eta) )
    dydx[1]=0
    d2ydx2[1]=(2*y[2] - 2*y[1]) / dx / dx
    d2ydx2[end]=(dydx[end]-dydx[end-1])/dx
    dydtau = d2ydx2 #+ 2/x_nodes*dydx
    #print(j, Ucollector,eta,cs,dydx[-1],d2ydx2[-1])

    #Create the arrays that will save the chronological results to disk, and viewed later for insights
    saved_time_spacing = saved_time_spacing  #It was defined above
    tau_array_saved = cat([1.,2.,4.],6:saved_time_spacing:size(tau_nodes)[1],dims=1)
    y_saved=zeros(size(tau_array_saved)[1],size(x_nodes)[1])
    i_saved=zeros((size(tau_array_saved)[1]))
    Ucollector_saved=zeros((size(tau_array_saved)[1]))
    y_saved[1,:]=y
    i_saved[1]=-F*z*D*dydx[end]*ct/r0
    Ucollector_saved[1]=Ucollector
    k=1

    for j in 2:size(tau_nodes)[1]   #loop over time
        Ucollector=Ucollector+v*dtau
        if Ucollector >= Utop && sign(v)==1
            v=-v
            Ucollector=Ucollector+v*dtau*2
        end
        if Ucollector <= Ubottom && sign(v)==-1
            v=-v
            Ucollector=Ucollector+v*dtau*2
        end
        y[:]=y[:]+dydtau*dtau                #update the concentration field inside the pcl
        cs=y[end]*ct                          #calculate the surface concentraction of intercalated ions
        eta = Ucollector-Ueq(y[end])
        for i in 2:size(x_nodes)[1]-1
            dydx[i] = (y[i+1] - y[i-1]) / 2 / dx
        end
        #I used this following equation for deriving the dydx[-1] equation seen below: Ueq=4.0 + R*T/F*log((ct/cs-1)/cl)
        dydx[end]=r0*z*kb*(ct-cs)/cl/ct/D*(exp(-beta*R*T/F*eta) - exp((1-beta)*R*T/F*eta) )
        #print(Ucollector,dydx[-1])
        dydx[1]=0
        for i in 2:size(x_nodes)[1]-1
            d2ydx2[i] = (y[i-1] - 2*y[i] + y[i+1]) / dx / dx
        end
        d2ydx2[end]=(dydx[end]-dydx[end-1])/dx
        d2ydx2[1]=(2*y[2] - 2*y[1]) / dx / dx
        dydtau = d2ydx2 #+ 2/x_nodes*dydx

        if any(j.==tau_array_saved)
            k=k+1
            #println(k)
            #print(Ucollector,eta,y[1],y[end-2],y[end])
            #println(j, Ucollector,eta,cs,dydx[end],d2ydx2[end])
            y_saved[k,:]=y
            i_saved[k]=-F*z*D*dydx[end]*ct/r0
            Ucollector_saved[k]=Ucollector
        end
    end
    print("elapsed time is ",r0*r0/D*tau_nodes[end], " seconds\n")
end

using Plots
main()
