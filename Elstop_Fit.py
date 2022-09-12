#This file performs conversions from digitised Antiproton data sets into a constant electronic stopping power profile
#Allows for many data points to be used
#Uses multiple fitting routines, starting from polynomial fitting with straight line extrapolaions, and later using a smoothing spline to remove any experimental noise

import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy import interpolate

Include_ZBL = 0
ZBL_Vel_Cutoff = 0.4e8
ZBL_File ="My_Elstop.in"
#from scipy.optimize import curve_fit
#smoothness = 17
#transition=1.2e7

foil = ['Ti','Ag','Au','Cu','Si','Al']
smoothnesses = [10,17,60,20,1,1]
transitions = [1.2e7, 1.2e7,1.2e7,1.2e7,1.2e7,1.2e7]


def exp(x, a, b):
    return a*np.exp(b*x)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def poly_lin_fit(x, x_data, y_data, poly_n, transition):
    """Fit a compound curve to the given data, return values over range
x.
    The curve consists of a straight line through the origin,
    a polynomial fit of degree poly_n,
    and an exponential.
    User selected transition marks the point where polynomial ends
    and exponential begins."""
   
    i_transition = find_nearest(x_data, transition)
    z_lin1 = np.polyfit((0.0, x_data[1]), (0.0, y_data[1]), 1) # from origin to first point
    p_lin1 = np.poly1d(z_lin1)
    y_lin1 = p_lin1(x[x <= x_data[1]])
   
    z_poly = np.polyfit(x_data[1:], y_data[1:], 7)
    p_poly = np.poly1d(z_poly)   
    y_poly = p_poly(x[(x >= x_data[1]) & (x <= x_data[i_transition])])

    popt, pcov = curve_fit(exp, x_data[i_transition:], y_data[i_transition:], p0=(4, -0.4e-7)) # initial guesses for this example p0=(4, -0.4e-7)) 
    y_exp = exp(x[x >= x_data[i_transition]], *popt)
   
    y_tot = np.concatenate((y_lin1, y_poly, y_exp))
    return y_tot

def splinefit(xnew, xdata, ydata,smoothness=8):
    xdata,ydata = zip(*sorted(zip(xdata,ydata)))
    #print(xdata)
    #print(ydata)
    tck = interpolate.splrep(xdata,ydata,s=smoothness)
    spline =  interpolate.splev(xnew,tck)
    return(spline)

def do_elstop(foil,transition,smoothness):
    shared = "/media/sf_MDrange_Shared/"
    El_digi = np.genfromtxt(shared+"Pbar_"+foil+"_Dataset.txt",delimiter=",")   #Digitised data set in its native units, please see Thesis entitled
    #"COMPLETE END TO END ANTIPROTON SIMULATIONS OF TRANSPORT,DEGRADATION AND EARLY TRAPPING" for detailed explanation of where the data comes from.

    El_digi[:,1] = El_digi[:,1] / 10 
    proton_mass = 1.6726219e-27
    kev_2_j = 1.60218e-16
    El_digi[:,0] = np.sqrt((2*( El_digi[:,0] * kev_2_j))/proton_mass)

    ZBL_Vel_Cutoff = max(El_digi[:,0])
    #print(El_digi)
    if Include_ZBL == 1:
        zbl = np.genfromtxt(ZBL_File,unpack=True)
        for i in range(len(zbl[0,:])):
            if zbl[0,i] >= ZBL_Vel_Cutoff:
                El_digi = np.vstack((El_digi,zbl[:,i]))

                #print(El_digi)
        #np.savetxt(shared+sys.argv[1],El_digi,delimiter=" ")#
    df = pd.DataFrame(El_digi,columns=["x","y"])
    print(df)


    df_array = df.to_numpy()

    x = np.linspace(0.0, np.max(df_array[:,0]), num=100)
    #x_long = np.linspace(0.0, np.max(df_array[:,0])+1e7, num=3000)
    x_long = np.linspace(0.0, 1e8,num=3000)


    y = poly_lin_fit(x_long, df_array[:,0], df_array[:,1], 7, transition)
    #y = poly_lin_fit(x_long, df_array[:,0], df_array[:,1], 3, 0.65e7)
    """
    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel('Velocity (m/s)')
    ax.set_ylabel(r'$\mathrm{Electronic\ stopping\ (eV/\AA)}$')
    ax.scatter(df_array[:,0], df_array[:,1], color='C0',marker="+")
    ax.plot(x_long, y, linestyle='--', color='C3',label="Poly+Exp")
    #fig.savefig('./el_stop_velocity_Nordlund_lin-poly-exp_fit.jpeg',dpi=150, bbox_inches='tight')
    plt.show()
    """


    #Fid = open(shared+"Python_Fit_SiPbar_No_ZBL.in","w");
    #Fid = open(shared+"Python_Fit_AlPbar_No_ZBL.in","w");
    #Fid = open(shared+"Python_Fit_AuPbar_No_ZBL.in","w");
    #Fid = open(shared+"Python_Fit_CuPbar_No_ZBL.in","w");
    #Fid = open(shared+"Python_Fit_TaPbar_No_ZBL.in","w");
    #Fid = open(shared+"Python_Fit_AgPbar_No_ZBL.in","w");
    #Fid = open(shared+"Python_Fit_TiPbar_No_ZBL.in","w");

    
    Fid = open(shared+"Python_Fit_"+foil+"Pbar_No_ZBL.in","w");
    
    splined=  splinefit(x_long,x_long,y,smoothness=smoothness)
    splined[0] = 0.00
    print(splined)
    
    fig1,ax2 = plt.subplots(figsize=(6.25,4))
    ax2.set_xlabel('Velocity (m/s)')
    ax2.set_ylabel(r'$\mathrm{Electronic\ stopping\ (eV/\AA)}$')
    ax2.scatter(df_array[:,0], df_array[:,1], color='black',marker="+",label="Experimentally measured datapoints")
    ax2.plot(x_long,splined,linestyle="-",color="royalblue",label="Secondary spline fit",alpha=0.6)
    ax2.plot(x_long,y,linestyle="--",color="red",label="Polynomial + exponential fit",alpha=0.6)
    ax2.legend(frameon=False)
    ax2.set_title("$\overline{p}$-"+foil+" electronic stopping power")
    ax2.set_xlim([0,1e8])
    ax2.set_ylim([0,max(y)*1.1])
    fig1.tight_layout()
    plt.show()


    for i in range(len(x_long)):
        # buf = "   %6.11f       %1.16f     \r\n"%(x_long[i],y[i])
        buf = "   %6.11f       %1.16f     \r\n"%(x_long[i],splined[i])
        #print(buf)
        #print(i)
        Fid.write(buf)
    Fid.close()


for i,j,k in zip(foil,transitions,smoothnesses):
    do_elstop(i,j,k)
