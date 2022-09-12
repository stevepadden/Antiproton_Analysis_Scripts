#This file performs the gaussian analysis to find the optimum thickness for 5 keV transmission
# A number of datasets are used, found by iterating around the optimum thickness from Range_Find.py
#This produces a gaussian distribution, the modal centre of this gaussian is our new optimised thickness
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter,NullFormatter, MaxNLocator
from numpy import linspace

def pickle_read(path):
    Sample=True
    global Sample_Size
    df = pd.read_pickle(path)
    if df.shape[0] <= Sample_Size and Sample == True:
        Sample_Size = df.shape[0]
    return df
        
def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

def Quick_Hist(dataframe,name,marker='+',**kwargs):
    Sample=False
    if Sample==False:
        dataframe = dataframe.sample(Sample_Size)
    else:
        dataframe = dataframe
        
    try:
        counts, bins = np.histogram(dataframe['Depth(x)'],bins=300)
        average = np.average(dataframe['Depth(x)'])
    except:
        counts, bins = np.histogram(dataframe['range'],bins=300)
        average = np.average(dataframe['range'])
    counts = zero_to_nan(counts)
    buf = "%s , Average:%.2e $\AA$" % (name,average)
    ax.scatter(bins[:-1],counts,label=buf,marker=marker)


Sample_Size = 100000
pickle_folder = "/home/mdrange/Pickle_Files/"
Pbar_Si_New_Elstop = pickle_read(pickle_folder+"Pbar_Si_New_Elstop.pkl")
Pbar_Si_5kev = pickle_read(pickle_folder+"pbar_si_5kev.pkl")
fig,ax = plt.subplots()
Quick_Hist(Pbar_Si_New_Elstop,"New Elstop")
Quick_Hist(Pbar_Si_5kev,"5 Kev")
plt.legend(loc='upper left',frameon='False',fontsize=9)
plt.show()



average5 = np.average(Pbar_Si_5kev['range'])
print(Pbar_Si_5kev['energy'])
average100 = np.average(Pbar_Si_New_Elstop['range'])
Suggested_width = average100-average5
print("Average width @ 100 kev: %s \n Average width @ 5 kev : %s \n Suggested Width : %s "%(average100,average5,Suggested_width))
