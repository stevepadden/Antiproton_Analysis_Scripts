#This file initially plots the depth profile of 100 keV antiprotons into infinite depth silicon.
#Following this is then plots a 5keV antiproton beam depth profile into silicon, the suggested thickness for 100 keV -> 5 keV degredation is thus
#calculated and plotted nicely for visualisation. This is not the true optimised depth, but used as a starter for later gaussian analysis.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gaussians import *
from lmfit.models import SkewedGaussianModel

pd.set_option('display.max_columns', None)  #Setting pandas to display as many columns as are present

folder = r"/home/mdrange/Pickle_Files/" #Location of the MDRange returned datasets
prefix = "Silicon/Silicon_End_" #File prefix
suffix = ".pkl" #File suffix

datasets = ["Range"]    #Which data sets we wish to look at within the pickle structures
names = [r"100 keV $\bar{p}$"]
model = SkewedGaussianModel()

#Binning and creating a histogram object of antiprotons
def Depth(fname,name):
    fname = folder+prefix+fname+suffix
    data = pd.read_pickle(fname)
    print(data)
    h1 =ax.hist(data["sz"],bins=100,edgecolor=(0,0,0,1),color="royalblue",label=name)
    print(h1)
    return h1   
#Plotting the histogram
fig,ax = plt.subplots()
xx = []
for i,j in zip(datasets,names):
    h1=Depth(i,j)
    xx.append(h1)

ax.set_xlabel("Range ($\AA$)")
ax.set_ylabel("Counts")
ax.set_title(r"Range profile of 100 keV $\bar{p}$ into silicon")

ax.legend(frameon=False)
#First guesses for the skew gaussian
Amp=200
Cen=15000
Wid=2400
g = -2.4
#Fitting the skew gaussian
model = SkewedGaussianModel()
params = model.make_params(amplitude=Amp,center=Cen,sigma=Wid,gamma=g)



h1 = xx[0]
y =np.array(h1[0],dtype=int)
x = np.array(h1[1],dtype=int)
x = x[:-1]

result = model.fit(y,params,x=x)

xnew = np.arange(min(x),max(x),500)
#Returning fit params
print(result.best_fit)
fit = result.best_fit
R = 1-result.residual.var() / np.var(y)
maxindex = fit.argmax()
mode = x[maxindex]
ax.plot(x,result.best_fit,label="Skewed Gaussian Fit E=100keV:\n$\mu$=%.2f\n$\sigma$=%2.f\n$A$=%.2f\n$\gamma$=%.2f\n$R^2$=%.3f"%(result.params['center'],result.params['sigma'],result.params['amplitude'],result.params['gamma'],R),color="black",linestyle="--")
cen = result.params['center']
ax.axvline(mode,color="red",label="Modal center 100keV=%.2f $\AA$"%mode)
#ax.plot(xnew,result.best_fit)
ax.set_xlim([0,max(x)*1.02])
ax.legend(frameon=False)
plt.show()

#Doing the same as above but now using the 5 keV dataset, if you dont have a 5 keV data set, set this to False.
kev_5 = True
#5kev stuff
fig,ax = plt.subplots()
if kev_5==True:
    data = pd.read_pickle(folder+"pbar_si_5kev.pkl")
    sz = data["sz"][data["sz"]>10]
    h2 =ax.hist(sz,bins=30,edgecolor=(0,0,0,1),color="indianred",label=r"5keV $\bar{p}$")
    Amp=200
    Cen=100
    Wid=600
    g = 0
    params2 = model.make_params(amplitude=Amp,center=Cen,sigma=Wid,gamma=g)
    y1 =np.array(h2[0],dtype=int)
    x1 = np.array(h2[1],dtype=int)
    x1 = x1[:-1]
    print(x1)
    print(y1)
    result1 = model.fit(y1,params2,x=x1)
    R2 = 1-result1.residual.var() / np.var(y1)

    ax.plot(x1,result1.best_fit,label="Skewed Gaussian Fit E=5keV:\n$\mu$=%.2f\n$\sigma$=%2.f\n$A$=%.2f\n$\gamma$=%.2f\n$R^2$=%.3f"%(result1.params['center'],result1.params['sigma'],result1.params['amplitude'],result1.params['gamma'],R2),color="black",linestyle="-.")
    
    print(result1.fit_report())
    fit1 = result1.best_fit
    maxindex = fit1.argmax()
    mode1 = x1[maxindex]

    diff = mode-mode1   #This is the suggested thickness!
    ax.axvline(mode1,color="gold",label="Modal center 5keV=%.2f $\AA$"%mode1)

    ax.annotate(
    '', xy=(mode1, 300), xycoords='data',
    xytext=(mode, 300), textcoords='data',
    arrowprops={'arrowstyle': '<->'})

    ax.annotate("Suggested thickness : %i $\AA$"%diff,xy=(max(x)/2, 301),ha='center')

    from matplotlib.legend import Legend
    #leg = Legend(ax,lines=["-."],loc="best")
    #ax.add_artist(leg)
    for i,j in zip(datasets,names):
        h1=Depth(i,j)
    ax.plot(x,result.best_fit,label="Skewed Gaussian Fit E=100keV:\n$\mu$=%.2f\n$\sigma$=%2.f\n$A$=%.2f\n$\gamma$=%.2f \n$R^2$=%.3f"%(result.params['center'],result.params['sigma'],result.params['amplitude'],result.params['gamma'],R),color="black",linestyle="--")
    cen = result.params['center']
    ax.axvline(mode,color="red",label="Modal center 100keV=%.2f $\AA$"%mode)
    #ax.plot(xnew,result.best_fit)
    ax.set_xlim([0,max(x)*1.02])
    ax.set_xlabel("Range ($\AA$)")
    ax.set_ylabel("Counts")
    ax.set_title(r"Range profile of 5 and 100 keV $\bar{p}$ into silicon")
    
    ax.legend(ncol=4,loc="upper center",frameon=False)
ax.set_xlim([10,max(x)*1.02])


plt.show()

