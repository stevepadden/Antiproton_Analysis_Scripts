#Creates a heatmap of proton/antiproton spatial positions upon leaving a foil of a given thickness
#Overlays both onto one spatial heatmap with edge distributions for clarity,
#Useful to show the differences in either 2 foil types, or too show the difference between 2 particle species traversing foils

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter, MaxNLocator

from numpy import linspace
import seaborn as sns

mpl.rcParams["legend.frameon"] = False
mpl.rcParams["legend.loc"] = "upper left"
srimpath = "/media/sf_MDrange_Shared/TRANSMIT_100_5_8571_SI.txt"    #Srim transmitted dataset
srim_ini_num = 20000
srimrange = 8571
md_range = 15044
cmap_pbar = "Reds"
cmap_proton="Blues"
elcol="blue"
facecol = "#faf6f5"
#facecol=""
n_replicas = 50 # Number of samples used in boostrap resampling
bsfrac = 0.05 #Fraction used in each resample



#Cleaning up the SRIM data
with open(srimpath,'r') as file:
    text = file.read()
clean = text.replace("T","")

with open("cleaned.txt","w+") as file:
    file.write(clean)



transmitted = np.genfromtxt("cleaned.txt",skip_header=12,invalid_raise=False,dtype="float")
transcols = ["Ion Numb","Atom Numb","Energy","Depth_X","Lateral_Y","Lateral_Z","Cos_X","Cos_Y","Cos_Z"]
Srim = pd.DataFrame(transmitted,columns=transcols)
print(Srim)

mdpath = "/home/mdrange/Pickle_Files/Si/Si_15044_nogauss.pkl"   #Pickled dataset resulting from Molecular Dynamics simulations



md = pd.read_pickle(mdpath)
md_transmitted = md[md["sz"]>md_range]
pbarrmsxarr = []
pbarrmsyarr = []
protonrmsxarr = []
protonrmsyarr = []
for i in range(n_replicas):
    md_transmitteds = md_transmitted.sample(frac=bsfrac,replace=True)
    Srims = Srim.sample(frac=bsfrac,replace=True)
    heatmap,xedges,yedges = np.histogram2d(md_transmitteds["sx"],md_transmitteds["sy"],bins=50)
    print(heatmap)
    print(np.max(heatmap))
    heatmap = np.divide(heatmap,np.max(heatmap))
    pbar_rmsxs = np.sqrt(np.sum(np.power(md_transmitteds["sx"],2) )/len(md_transmitteds["sx"]) )
    pbar_rmsys = np.sqrt(np.sum(np.power(md_transmitteds["sy"],2) )/len(md_transmitteds["sy"]) )
    heatmap,xedges,yedges = np.histogram2d(Srims["Lateral_Y"],Srims["Lateral_Z"],bins=50)
    heatmap = np.divide(heatmap,np.max(heatmap))
    proton_rmsxs = np.sqrt(np.sum(np.power(Srims["Lateral_Y"],2) )/len(Srims["Lateral_Y"]) )
    proton_rmsys = np.sqrt(np.sum(np.power(Srims["Lateral_Z"],2) )/len(Srims["Lateral_Z"]) )
    pbarrmsxarr.append(pbar_rmsxs)
    pbarrmsyarr.append(pbar_rmsys)
    protonrmsxarr.append(proton_rmsxs)
    protonrmsyarr.append(proton_rmsys)

pbarrmsxarr = np.array(pbarrmsxarr)
pbarrmsyarr = np.array(pbarrmsyarr)
protonrmsxarr = np.array(protonrmsxarr)
protonrmsyarr = np.array(protonrmsyarr)

pbarrmsx = pbarrmsxarr.mean()
pbarrmsy = pbarrmsyarr.mean()
protonrmsx = protonrmsxarr.mean()
protonrmsy = protonrmsyarr.mean()

sdrmsxpb = pbarrmsxarr.std()
sdrmsypb = pbarrmsyarr.std()
sdrmsxpr = protonrmsxarr.std()
sdrmsypr = protonrmsyarr.std()

N = len(pbarrmsxarr)
errpbrmsx = sdrmsxpb/np.sqrt(N)
errpbrsmy = sdrmsypb/np.sqrt(N)
errprrmsx = sdrmsxpr/np.sqrt(N)
errprrmsy = sdrmsypr/np.sqrt(N)

print(pbarrmsx)
print(errpbrmsx)


def ellipse(ra,rb,ang,x0,y0,Nb=100):
    xpos,ypos=x0,y0
    radm,radn=ra,rb
    an=ang
    co,si=np.cos(an),np.sin(an)
    the=np.linspace(0,2*np.pi,Nb)
    X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
    Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
    return X,Y

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

print(md)
print(md_transmitted)


#step1 load datasets

mdpath = "/home/mdrange/Pickle_Files/Si/Si_15044_nogauss.pkl"   #Reloading the dataset for antiprotons
md = pd.read_pickle(mdpath)
pbar=md[md["sz"]>md_range]
pbar= pbar.sample(10000)

proton = pd.DataFrame(transmitted,columns=transcols)
proton = proton.sample(10000)
data1 = {"Particle species":"proton","sx":proton["Lateral_Z"],"sy":proton["Lateral_Y"]}
data2 = {"Particle species":"antiproton","sx":pbar["sx"],"sy":pbar["sy"]}
df1  = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)

print(df1)
print(df2)
df = pd.concat([df2,df1],ignore_index=True)
print(df)



pbar_cmap = "Blues"
pbar_contours = "blue"
proton_cmap= "Reds"
proton_contours="red"

pbar_rmsx = pbarrmsx
pbar_rmsy = pbarrmsy

proton_rmsy = protonrmsy
proton_rmsx = protonrmsx


pbar_string = "$\overline{p}$ : $\sigma_{rmsx}$ : %.2f$\pm$%.2f ($\AA$) : $\sigma_{rmsy}$ %.2f$\pm$%.2f ($\AA$)" % (pbar_rmsx,errpbrmsx,pbar_rmsy,errpbrsmy)
proton_string = "protons : $\sigma_{rmsx}$ : %.2f$\pm$%.2f ($\AA$) : $\sigma_{rmsy}$ %.2f$\pm$%.2f ($\AA$)" % (proton_rmsx,errprrmsx,proton_rmsy,errprrmsy)


palette = {'proton' : "tab:red", "antiproton":"tab:blue"}
g = sns.jointplot(data=df,x="sx",y="sy",hue="Particle species",kind="hist",palette=palette,joint_kws={"alpha":0.7},marginal_kws={"palette":palette})

ax=g.ax_joint
g.fig.suptitle("Spatial heatmaps of 10000 simulated $\overline{p}$ and protons\n leaving an optimised silicon degrader foil")
g.fig.tight_layout()
g.ax_joint.set_xlabel("x ($\AA$)")
g.ax_joint.set_ylabel("y ($\AA$)")
lim = 8500
g.ax_marg_x.set_xlim(-lim,lim)
g.ax_marg_y.set_ylim(-lim,lim)

labels=[pbar_string,proton_string]
texts = g.ax_joint.legend_.texts
for t,label in zip(texts,labels):
    t.set_text(label)

g.ax_joint.legend_.set_title(None)


pbar_sdx = np.std(pbar["sx"])
pbar_sdy = np.std(pbar["sy"])
pbar_xcen = np.mean(pbar["sx"])
pbar_ycen = np.mean(pbar["sy"])

proton_sdx = np.std(proton["Lateral_Y"])
proton_sdy = np.std(proton["Lateral_Z"])
proton_xcen = np.mean(proton["Lateral_Y"])
proton_ycen = np.mean(proton["Lateral_Z"])

rotation = 0
pbar_X1,pbar_Y1 = ellipse(pbar_sdx,pbar_sdy,rotation,pbar_xcen,pbar_ycen)

ax.plot(pbar_X1,pbar_Y1,"k:",ms=1,linewidth=2.0,zorder=7,color=pbar_contours)
xco = np.argmin(pbar_X1)
yco = np.argmin(pbar_Y1)


pbar_X2,pbar_Y2 = ellipse(pbar_sdx*2,pbar_sdy*2,rotation,pbar_xcen,pbar_ycen)
ax.plot(pbar_X2,pbar_Y2,"k:",ms=1,linewidth=2.0,zorder=8,color=pbar_contours)

pbar_X3,pbar_Y3 = ellipse(pbar_sdx*3,pbar_sdy*3,rotation,pbar_xcen,pbar_ycen)
ax.plot(pbar_X3,pbar_Y3,"k:",ms=1,linewidth=2.0,zorder=9,color=pbar_contours)


proton_X1,proton_Y1 = ellipse(proton_sdx,proton_sdy,rotation,proton_xcen,proton_ycen)

ax.plot(proton_X1,proton_Y1,"k:",ms=1,linewidth=2.0,zorder=2,color=proton_contours)

proton_X2,proton_Y2 = ellipse(proton_sdx*2,proton_sdy*2,rotation,proton_xcen,proton_ycen)
ax.plot(proton_X2,proton_Y2,"k:",ms=1,linewidth=2.0,zorder=3,color=proton_contours)

proton_X3,proton_Y3 = ellipse(proton_sdx*3,proton_sdy*3,rotation,proton_xcen,proton_ycen)
ax.plot(proton_X3,proton_Y3,"k:",ms=1,linewidth=2.0,zorder=4,color=proton_contours)


plt.show()

