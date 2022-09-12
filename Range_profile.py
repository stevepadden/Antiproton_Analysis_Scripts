import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gaussians import *
from lmfit.models import SkewedGaussianModel

pd.set_option('display.max_columns', None)

folder = r"/home/mdrange/Pickle_Files/"
prefix = "Silicon/Silicon_End_"
suffix = ".pkl"

datasets = ["Range"]
names = [r"100 keV $\bar{p}$"]
model = SkewedGaussianModel()


def Depth(fname,name):
    fname = folder+prefix+fname+suffix
    data = pd.read_pickle(fname)
    print(data)
    h1 =ax.hist(data["sz"],bins=100,edgecolor=(0,0,0,1),color="royalblue",label=name)
    print(h1)
    return h1

fig,ax = plt.subplots()
xx = []
for i,j in zip(datasets,names):
    h1=Depth(i,j)
    xx.append(h1)

ax.set_xlabel("Range ($\AA$)")
ax.set_ylabel("Counts")
ax.set_title(r"Range profile of 100 keV $\bar{p}$ into silicon")
#nordlunddat = pd.read_pickle(folder+"endrecoil.pkl")
#print(nordlunddat)
#ax.hist(nordlunddat["sz"],bins=100,alpha=0.5,edgecolor=(0,0,0,1))
ax.legend(frameon=False)

#plt.show()


#model = Model(gaussian)
#model = Model(bimodal)
#a2 = 200
#c2 = 1000
#w2 = 400
Amp=200
Cen=15000
Wid=2400
g = -2.4

model = SkewedGaussianModel()
params = model.make_params(amplitude=Amp,center=Cen,sigma=Wid,gamma=g)



h1 = xx[0]
y =np.array(h1[0],dtype=int)
x = np.array(h1[1],dtype=int)
x = x[:-1]
print(y)
print(x)
#result = model.fit(y,x=x, amp1=Amp, cen1=Cen, wid1=Wid,amp2=a2,wid2=w2,cen2=c2)
#result = model.fit(y,x=x, amp=Amp, cen=Cen, wid=Wid)
result = model.fit(y,params,x=x)

#tamp = result.params['amp'].value
#twid = result.params['wid'].value
#tcen = result.params['cen'].value
#a1 = result.params['amp1'].value
#a2 = result.params['amp2'].value
#c1 = result.params['cen1'].value
#c2 = result.params['cen2'].value
#w1 = result.params['wid1'].value
#w2 = result.params['wid2'].value

#xnew=np.arange(min(x),max(x),50)
xnew = np.arange(min(x),max(x),500)
#const_fit = gaussian(xnew, tamp, tcen, twid)  # y values from model
#const_fit = bimodal(xnew, a1,c1,w1,a2,c2,w2)
print(result.fit_report())
#hp = ax.plot(xnew,const_fit,color="red")
#ax.plot(x,result.init_fit)
#const_fit = SkewedGaussianModel(x=xnew,amplitude=result.params['amplitude'],center=result.params['center'],sigma=result.params['sigma'],gamma=result.params['gamma'])
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

    diff = mode-mode1
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

