#This file does all analysis on an output foil, it produces multiple graphs such as;
#Energy distributions, angular distributions, Range profiles and spatial heatmaps
#This file is used extensivley in checking the antiproton distributions - it is strongly suggested to understand what is happening here before using it.
#Bootstrap resampling is used to find errors related to given quantities, this is found at the bottom end of the code. It is not explained how to bootstrap resample
#within this code, however the infomation here in conjunction with some prior knowledge should be sufficient to follow and adjust the code
#as required.


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit.models import SkewedGaussianModel
pd.set_option('display.max_columns', None)

Foil_type = "Si"
Foil_fullname = "silicon"
optimised_thick = 15044 #Thickness used
srim_thick = 8571   #Thickness from SRIM
capture_energy = 5000   #Max trappable energy
#importing data
suffix1 = "/home/mdrange/Pickle_Files/" #Pkl file location (ignore it being called suffix, this is related to another code)
#Range Data
Range_suffix = "Range_Profiles/"
range_path = suffix1 + Range_suffix+Foil_type+".pkl"
range_data = pd.read_pickle(range_path)
figsize = [6.25,4]
figsize = np.array(figsize)
figsize = figsize * 1.4
cmap="viridis"
savestring = "/media/sf_MDrange_Shared/Si_Gauss_Images/" #Where the images save - make sure this is adjusted 
#Optimised Data
opt_suffix = "Si/" 
opt_path = suffix1+opt_suffix+Foil_type+"_%s_withgauss.pkl"%optimised_thick
opt_data = pd.read_pickle(opt_path)

nogauss = pd.read_pickle(suffix1+opt_suffix+Foil_type+"_%s_nogauss.pkl"%optimised_thick)
#transit_data = opt_data[opt_data["sz"]>optimised_thick]
nogauss_transit = nogauss[nogauss["sz"]>optimised_thick]

#histogram the range profile
fig,ax = plt.subplots(figsize=figsize)
ax.set_xlabel("Range ($\AA$)")
ax.set_ylabel("Counts")
ax.set_title("Range profile of 100 keV $\\bar{p}$ into %s \n with G4Beamline input parameters"%Foil_type)
Num = range_data["sz"].count()
num_past = range_data["sz"][range_data["sz"]>optimised_thick].count()
percent = num_past/Num * 100
past_srim =range_data["sz"][range_data["sz"]>srim_thick].count()
psrim = past_srim/Num * 100
h1=ax.hist(range_data["sz"],bins=100,edgecolor=(0,0,0,1),color="royalblue",label="100 keV $\overline{p}$ into %s"%Foil_type)
ax.axvline(optimised_thick,color="red",lw=2, label=("Optimised $\overline{p}$-%s thickness: %s $\AA$"% (Foil_type,optimised_thick)))
tstr="Optimised proton thickness: %i $\AA$\n%.2f%% of $\overline{p}$ travel past\nproton thickness"%(srim_thick,psrim)
ax.axvline(srim_thick,color="grey",ls="--",lw=2,alpha=0.6,label=tstr)
#ax.plot([],[],' ', label=("%.2f %% travel further than \noptimised thickness"% percent))
#ax.title("%.2f %% travel further than \noptimised thickness"% percent)
Amp=1734367
Cen=16000
Wid=2480
g =-3
model = SkewedGaussianModel()
params = model.make_params(amplitude=Amp,center=Cen,sigma=Wid,gamma=g)
y=np.array(h1[0],dtype=int)
x = np.array(h1[1],dtype=int)
x = x[:-1]
result = model.fit(y,params,x=x)
print(result.fit_report())
xnew = np.arange(min(x),max(x),500)
R = 1-result.residual.var() / np.var(y)
ax.plot(x,result.best_fit,label="Skewed Gaussian Fit:\n$\mu$=%.2f\n$\sigma$=%2.f\n$A$=%.2f\n$\gamma$=%.2f\n$R^2$=%.3f"%(result.params['center'],result.params['sigma'],result.params['amplitude'],result.params['gamma'],R),color="black",linestyle="--")


ax.legend(frameon=False,title="%.2f %% travel further than \noptimised $\overline{p}$-%s thickness"% (percent,Foil_type))
ax.set_xlim([0,1.02*max(range_data["sz"])])

fig7,ax7 = plt.subplots()
ax7.hist(opt_data["sz"],bins=50,label="Range profile with input distribution",edgecolor=(0,0,0,1),alpha=0.7,color="royalblue")
ax7.hist(nogauss["sz"],bins=50,label="Range profile without input distribution",edgecolor=(0,0,0,1),alpha=0.7)
ax7.set_xlabel("Range ($\AA$)")
ax7.set_ylabel("Counts")
ax7.set_title("Range profile comparisson based on input distribution \n for 100 keV $\overline{p}$ in %s"%Foil_fullname)
ax7.legend(frameon=False)



#extract the successful transit dataframe
fig2,ax2 = plt.subplots(figsize=figsize)
ax2.set_xlabel("Energy of transmitted $\overline{p}$ (eV)")
ax2.set_ylabel("Counts")
ax2.set_title("MDrange 100 keV $\overline{p}$ into %i $\AA$ %s \n with G4Beamline input parameters"%(optimised_thick,Foil_fullname))
transit_data = opt_data[opt_data["sz"]>optimised_thick]
num_sim = opt_data["sz"].count()
print(transit_data)
num_trappable = transit_data["energy"][transit_data["energy"]<capture_energy].count()
num_trans = transit_data["energy"].count()
per_trappable = num_trappable/num_trans *100
tot = num_sim*(per_trappable/100) * (percent/100)
h1=ax2.hist(transit_data["energy"],bins=50,edgecolor=(0,0,0,1),color="royalblue",label="Energy of transited $\overline{p}$")
ax2.axvline(capture_energy,color="red",lw=2,label="Maximum trappable energy")
ax2.set_xlim([0,max(transit_data["energy"])*1.1])
Amp=2500
Cen=130
Wid=2000
g =1
model = SkewedGaussianModel()
params = model.make_params(amplitude=Amp,center=Cen,sigma=Wid,gamma=g)
y=np.array(h1[0],dtype=int)
x = np.array(h1[1],dtype=int)
x = x[:-1]
result = model.fit(y,params,x=x)
print(result.fit_report())
xnew = np.arange(min(x),max(x),500)
R = 1-result.residual.var() / np.var(y)
ax2.plot(x,result.best_fit,label="Skewed Gaussian Fit:\n$\mu$=%.2f\n$\sigma$=%2.f\n$A$=%.2f\n$\gamma$=%.2f\n$R^2$=%.3f"%(result.params['center'],result.params['sigma'],result.params['amplitude'],result.params['gamma'],R),color="black",linestyle="--")
ax2.legend(frameon=False,title=("%.2f %% of transmitted $\overline{p}$ trappable\n%i $\overline{p}$ in total from %i simulated"%(per_trappable,tot,num_sim)))



#plotting 2d heatmap spatial

def findmax(x, y):
    x1 = max(x)
    x2 = abs(min(x))
    y1 = max(y)
    y2 = abs(min(y))
    if x2 > x1:
        maxx = x2
        minx = -x2
    else:
        maxx = x1
        minx = -x1

    if y2 > y1:
        maxy = y2
        miny = -y2
    else:
        maxy = y1
        miny = -y1
    endlim = (max(maxx, maxy))
    return minx, maxx, miny, maxy, endlim

fig3,ax3 = plt.subplots()
heatmap,xedges,yedges = np.histogram2d(transit_data["sx"]*1e-7,transit_data["sy"]*1e-7,bins=50)
rmsx = np.sqrt(np.sum(np.power(transit_data["sx"]*1e-7,2) )/len(transit_data["sx"]) )
rmsy = np.sqrt(np.sum(np.power(transit_data["sy"]*1e-7,2) )/len(transit_data["sy"]) )
ax3.set_xlabel("x (mm)")
ax3.set_ylabel("y (mm)")
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
minx, max3x, miny, max3y, endlim = findmax(transit_data["sy"]*1e-7, transit_data["sy"]*1e-7)
ax3.set_xlim((-endlim, endlim))
ax3.set_ylim((-endlim, endlim))
tt = ax3.imshow(heatmap.T, extent=extent, origin='lower', cmap=cmap, aspect="auto")
#tt = ax3.imshow(heatmap.T, extent=extent, origin='lower', cmap=cmap, aspect="equal")
ax3.set_facecolor('#440154')
str = "$\sigma_{RMSX} =$ %.2f (mm) \n$\sigma_{RMSY} =$ %.2f (mm)" % (rmsx, rmsy)
ax3.plot([], [], ' ', label=str)
ax3.legend(frameon=False, facecolor='#440154', labelcolor="w")
cbar = plt.colorbar(tt, ax=ax3)
ax3.set_title("Spatial heatmap of $\overline{p}$ leaving MDRange simulated %s foil \n of %s $\AA$ thickness"% (Foil_fullname,optimised_thick))
fig3.set_tight_layout(tight=True)


print(transit_data)
#Calculate theta and phi
vx = transit_data['vx']
vy = transit_data['vy']
vz = transit_data['vz']
theta = np.arccos(vz/np.sqrt((vx**2) + (vy**2) + (vz**2)))
theta = np.degrees(theta)
phi = np.arctan2(vy,vx)# * 180 / np.pi
phi = np.degrees(phi)
#plot theta
fig4,ax4 = plt.subplots(figsize=figsize)
h1=ax4.hist(theta,bins=60,edgecolor=(0,0,0,1),color="royalblue",label="Theta")
#h1 = np.histogram(theta,bins=200)
ax4.set_title("$\\theta$ distribution resulting from %i $\AA$ of %s \n including G4Beamline start parameters"%(optimised_thick,Foil_fullname))
ax4.set_xlabel("$\\theta$ ($^{\circ}$)")
ax4.set_ylabel("Counts")
ax4.set_xlim([0,90])

Amp=22579
Cen=10.3
Wid=32
g =2.86

model = SkewedGaussianModel()
params = model.make_params(amplitude=Amp,center=Cen,sigma=Wid,gamma=g)
y=np.array(h1[0],dtype=int)
x = np.array(h1[1],dtype=int)
x = x[:-1]
result = model.fit(y,params,x=x,method="least_squares")
print(result.fit_report())
xnew = np.arange(0,90,500)

R = 1-result.residual.var() / np.var(y)
ax4.plot(x,result.best_fit,label="Skewed Gaussian Fit:\n$\mu$=%.2f\n$\sigma$=%2.f\n$A$=%.2f\n$\gamma$=%.2f\n$R^2$=%.3f"%(result.params['center'],result.params['sigma'],result.params['amplitude'],result.params['gamma'],R),color="black",linestyle="--")
ax4.legend(frameon=False)

print(R)
print(result.residual.var())
print(np.var(y))

fig5,ax5 = plt.subplots(figsize=figsize)
ax5.hist(phi,bins=50,edgecolor=(0,0,0,1),color="royalblue",label="Phi")
ax5.set_title("$\phi$ distribution resulting from %i $\AA$ of %s \n including G4Beamline start parameters"%(optimised_thick,Foil_fullname))
ax5.set_xlabel("$\phi$ ($^{\circ}$)")
ax5.set_ylabel("Counts")
ax5.set_xlim([-180,180])

fig6,ax6 = plt.subplots()

no_inp_mu = nogauss_transit["energy"].mean() 
no_inp_std = nogauss_transit["energy"].std()
no_inp_err = no_inp_mu/np.sqrt(no_inp_std)

inp_mu = transit_data["energy"].mean()
inp_std = transit_data["energy"].std()
inp_err = inp_mu/np.sqrt(inp_std)

if abs(inp_mu-no_inp_mu)< 3*np.sqrt(inp_err**2 + no_inp_err**2):
    print("RESULTS ARE CONSISTENT")

str1= "Energy profile with\ninput distribution:\n$\mu$: %i $\pm$ %i eV" % (inp_mu,inp_err)
str2 = "Energy profile without\ninput distribution:\n$\mu$: %i $\pm$ %i eV" % (no_inp_mu,no_inp_err)
ax6.hist(transit_data["energy"],bins=50,color="royalblue",alpha=0.7,edgecolor=(0,0,0,1),label = str1 ,zorder=2)
ax6.hist(nogauss_transit["energy"],bins=50,label=str2,color="red",alpha=0.7,edgecolor=(0,0,0,1))
ax6.set_title("Energy distributions of 100 kev $\overline{p}$ emerging from %i $\AA$ thick %s foil"%(optimised_thick,Foil_fullname))
ax6.set_xlabel("Energy (eV)")
ax6.set_ylabel("Counts")
ax6.legend(frameon=False)


#plot ranges of theta and no theta
rangesuffix = "/home/mdrange/Pickle_Files/Range_Profiles/"
ng = pd.read_pickle(rangesuffix+"Si.pkl")

g = pd.read_pickle(rangesuffix+"Range_gauss.pkl")
g = g.sample(20000)

no_inp_mu = ng["sz"].mean() 
no_inp_std = ng["sz"].std()
no_inp_err = no_inp_mu/np.sqrt(no_inp_std)

inp_mu = g["sz"].mean()
inp_std = g["sz"].std()
inp_err = inp_mu/np.sqrt(inp_std)

if abs(inp_mu-no_inp_mu)< 3*np.sqrt(inp_err**2 + no_inp_err**2):
    print("RANGE RESULTS ARE CONSISTENT")

str1= "Range profile without input distribution:\n$\mu$: %i $\pm$ %i $\AA$" % (no_inp_mu,no_inp_err)
str2= "Range profile with input distribution:\n$\mu$: %i $\pm$ %i $\AA$" % (inp_mu,inp_err)
    
fig8,ax8 = plt.subplots()
ax8.hist(ng["sz"],bins=50,color="red",alpha=0.7,edgecolor=(0,0,0,1),label=str1)
ax8.hist(g["sz"],bins=50,color="royalblue",alpha=0.7,edgecolor=(0,0,0,1),label=str2)
ax8.set_xlabel("Range ($\AA$)")
ax8.set_ylabel("Counts")
ax8.set_title("Range profile of $\overline{p}$ with and without input distributions")
ax8.set_xlim([0,max(g["sz"])*1.01])
ax8.legend(frameon=False)




fig.savefig(savestring+Foil_type+"_Range_GAUSS.png")
fig2.savefig(savestring+Foil_type+"_Energy_GAUSS.png")
fig3.savefig(savestring+Foil_type+"_Spatial_Heatmap_GAUSS.png")
fig4.savefig(savestring+Foil_type+"_Theta_GAUSS.png")
fig5.savefig(savestring+Foil_type+"_Phi_GAUSS.png")
fig6.savefig(savestring+Foil_type+"_Energy_gauss_diff.png")
#fig7.savefig(savestring+Foil_type+"_Range_Compare.png")
fig8.savefig(savestring+Foil_type+"_Range_Compare.png")

#plt.show()

#Bootstrap resampling for errors
n_replicas = 50
bsfrac = 0.05
print(range_data["energy"].count())
opt_data = pd.read_pickle(opt_path)
range_data = pd.read_pickle(range_path)
print(opt_data)
earray =[]
frac_array = []
per_trappable = []
rmsx_a = []
rmsy_a = []
for i in range(n_replicas):
    r =opt_data.sample(frac=bsfrac,replace=True) #Sample the data
    x = r[r["sz"]>optimised_thick] #Taking only those that pass the thickness

    pp = range_data.sample(frac=bsfrac,replace=True)
    Nparttrans = pp["sz"][pp["sz"]>optimised_thick].count()
    fractrans = Nparttrans/(pp["sz"].count())
    count = x["energy"][x["energy"]<capture_energy].count() #Count those that have an energy less than the capture energy
    rmsx = np.sqrt(np.sum(np.power(x["sx"]*1e-7,2) )/len(x["sx"]) )
    rmsy = np.sqrt(np.sum(np.power(x["sy"]*1e-7,2) )/len(x["sy"]) )
    #tot = r["energy"].count()
    tot = x["energy"].count()
    per_trap = count/tot * 100
    per_trappable.append(per_trap)
    per = count/tot * 100 * fractrans
    earray.append(per)
    frac_array.append(fractrans * 100)
    rmsx_a.append(rmsx)
    rmsy_a.append(rmsy)
earray = np.array(earray)
sd = earray.std()
mean=(earray.mean())
per_trappable= np.array(per_trappable)
rmsx_a = np.array(rmsx_a)
rmsy_a = np.array(rmsy_a)
                            

N = len(earray)
#err = sd/np.sqrt(mean)
err1 = sd/np.sqrt(N)
frac_array = np.array(frac_array)
err2 = frac_array.mean()/np.sqrt(N)
err = ((err1/mean)**2 + (err2/frac_array.mean())**2)**(1/2)


print(r"Average  %.2f $\pm$ %.2f" %( mean, err))

print(r"Percent trappable %.2f \pm %.2f" %(per_trappable.mean(), (per_trappable.std()/np.sqrt(N))))


print(r"frac trann trappable %.2f \pm %.2f" %(frac_array.mean(), (frac_array.std()/np.sqrt(N))))
rxa = rmsx_a.mean()
rxe = rmsx_a.std()/np.sqrt(N)
rya = rmsy_a.mean()
rye = rmsy_a.std()/np.sqrt(N)
print(r"RMS_x : %.2f $\pm$ %.4f : RMS_y : %.2f $\pm$ %.4f" % (rxa,rxe,rya,rye))

