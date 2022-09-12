import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))

def mgauss(x,a1,b1,c1):
    return (a1*np.exp(-((x-b1)/c1)**2))


def bimodal(x,amp1,cen1,wid1,amp2,cen2,wid2):
    return gaussian(x,amp1,cen1,wid1)+gaussian(x,amp2,cen2,wid2)

def histogram (data,name,numbins=200,marker='+',xlim=False,ylim=False,**kwargs):
    counts,bins = np.histogram(data,bins=numbins)
    average=np.average(data)
    buf = "%s , Average:%.2e $\AA$" % (name,average)
    fig,ax=plt.subplots()
    ax.set_ylabel("Counts")
    ax.set_xlabel(name)
    counts = zero_to_nan(counts)
    ax.scatter(bins[:-1],counts,label=buf,marker=marker)
    plt.legend(loc='upper left',frameon='False',fontsize=9)
    if xlim!=False:
        ax.xlim(xlim)
    if ylim!=False:
        ax.ylim(ylim)
    plt.show()
    plt.close()

def histogram_with_gaussian (data,name,Amp=1,Cen=1,Wid=1,numbins=200,marker='+',**kwargs):
    counts,bins = np.histogram(data,bins=numbins)
    average=np.average(data)
    #Gaussian fitting
    gmodel = Model(gaussian) #define the function returned model
    result = gmodel.fit(counts, x=bins[:-1], amp=Amp, cen=Cen, wid=Wid) #Fit the data (y, x to the model gmodel.fit)
    tamp = result.params['amp'].value
    twid = result.params['wid'].value
    tcen = result.params['cen'].value
    const_fit = gaussian(bins[:-1],tamp,tcen,twid) #y values from model
    buf = " Amplitude : %s \n Central : %s , Sigma : %s" %(tamp,tcen,twid)
    fig,ax=plt.subplots()
    plt.plot(bins[:-1],const_fit,label=buf)
    counts = zero_to_nan(counts)
    buf = "%s , Average:%.2e $\AA$" % (name,average)
    ax.set_ylabel("Counts")
    ax.set_xlabel(name)
    ax.scatter(bins[:-1],counts,label=buf,marker=marker)
    plt.legend(loc='upper left',frameon='False',fontsize=9)
    #plt.savefig(Shared_Path+"g4bl_"+name+".png")
    plt.show()
    
    plt.close()

    return result,counts,bins

def histogram_with_bimodal(data,name,Amp1=1,Cen1=1,Wid1=1,Amp2=1,Cen2=-1,Wid2=1,numbins=200,marker='+',**kwargs):
    counts,bins = np.histogram(data,bins=numbins)
    average=np.average(data)
    #Gaussian fitting
    gmodel = Model(bimodal) #define the function returned model
    result = gmodel.fit(counts, x=bins[:-1], amp1=Amp1, cen1=Cen1, wid1=Wid1, amp2=Amp2, cen2=Cen2, wid2=Wid2) #Fit the data (y, x to the model gmodel.fit)
    tamp1 = result.params['amp1'].value
    twid1 = result.params['wid1'].value
    tcen1 = result.params['cen1'].value
    tamp2 = result.params['amp2'].value
    twid2 = result.params['wid2'].value
    tcen2 = result.params['cen2'].value
    const_fit = bimodal(bins[:-1],tamp1,tcen1,twid1,tamp2,tcen2,twid2) #y values from model
    buf = " Amplitude1 : %s \n Central1 : %s \n Sigma1 : %s \n Amp2 : %s \n Cen2 : %s \n Sigma2 : %s" %(tamp1,tcen1,twid1,tamp2,tcen2,twid2)
    fig,ax=plt.subplots()
    plt.plot(bins[:-1],const_fit,label=buf)
    counts = zero_to_nan(counts)
    buf = "%s , Average:%.2e $\AA$" % (name,average)
    ax.set_ylabel("Counts")
    ax.set_xlabel(name)
    ax.scatter(bins[:-1],counts,label=buf,marker=marker)
    x= np.linspace(13.784,13.788,100)
    #y=bimodal(x,cen1=0.99,wid1=0.2,amp1=226,cen2=2.0,wid2=0.19,amp2=117)
    y=bimodal(x,cen1=13.7863,cen2=13.787,wid1=0.00025,wid2=0.00025,amp1=0.35,amp2=0.35)
    plt.plot(x,y)
    plt.legend(loc='upper left',frameon='False',fontsize=9)
    plt.show()
    plt.close()
    print(result.fit_report())

    return result,counts,bins

def TwoD_Hist(dataframe,name,arg1='px',arg2='py'):
    #dataframe = dataframe.sample(Sample_Size)
    sx = dataframe[arg1]
    sy = dataframe[arg2]
    xlims = [min(sx),max(sx)]
    ylims = [min(sy),max(sy)]
    xmin = min(xlims)
    xmax = max(xlims)
    ymin = min(ylims)
    ymax = max(ylims)
    sx = zero_to_nan(sx)
    sy = zero_to_nan(sy)
    nxbins = 50
    nybins = 50
    nbins = 100
    maxi = max(ymax,xmax)
    mini = min(ymin,xmin)
    xbins = np.linspace(start = xmin, stop = xmax, num = nxbins)
    ybins = np.linspace(start = ymin, stop = ymax, num = nybins)
    H, xedges,yedges = np.histogram2d(sy,sx,bins=(ybins,xbins))
    #xcenter = (xbins[0:-1]+xbins[1:])/2.0
    #ycenter = (ybins[0:-1]+ybins[1:])/2.0
    #aspectratio = 1.0*(xmax - 0)/(1.0*ymax - 0)
    #X = xcenter
    #Y = ycenter
    Z = H
    print(H.argmax())
    print(np.argwhere(H == H.max()))
    #axXY = plt.axes() # xy plot
    #cax = (axXY.imshow(H,interpolation='spline16',extent=[mini,maxi,mini,maxi],cmap="jet"))
    #plt.draw()
    #axXY.set_xlabel=arg1
    #axXY.set_ylabel=arg2
    #axXY.set_title="%s vs %s" %(arg1,arg2)
    extents = [xmin,xmax,ymin,ymax]
    #axXY.colorbar()
    #plt.show()
    plt.imshow(H,cmap="jet",interpolation='spline16',extent=[mini,maxi,mini,maxi])
    plt.colorbar()
    plt.xlabel(arg1)
    plt.ylabel(arg2)
    plt.title("%s vs %s" %(arg1,arg2))
    plt.show()
    return H,extents

def twohist_plot(data1,data2,name1,name2,xaxis="Momentum",marker1="+",marker2="+",numbins=200,**kwargs):
    counts1,bins1 = np.histogram(data1,bins=numbins)
    average1=np.average(data1)
    buf1 = "%s , Average:%.2e $\AA$" % (name1,average1)
    fig,ax=plt.subplots()
    ax.set_ylabel("Counts")
    ax.set_xlabel(xaxis)
    counts = zero_to_nan(counts1)
    ax.scatter(bins1[:-1],counts1,label=buf1,marker=marker1)
    counts2,bins2=np.histogram(data2,bins=numbins)
    average2=np.average(data2)
    buf2 = "%s , Average:%.2e $\AA$" % (name2,average2)
    ax.scatter(bins2[:-1],counts2,label=buf2,marker=marker2)
    plt.legend(loc='upper left',frameon='False',fontsize=9)
    plt.axvline(x=0.00)
    plt.show()

def Many_Hist(Dataarray,datatype,names,xaxis="energy",numbins=200,marker="+",**kwargs):
    #open figure
    fig,ax=plt.subplots()
    ax.set_ylabel("Counts")
    ax.set_xlabel(xaxis)
    i = 0
    for data in Dataarray:
        name = names[i]
        i = i+1
        counts,bins=np.histogram(data[datatype],bins=numbins)
        average=np.average(data[datatype])
        counts = zero_to_nan(counts)
        buf = "%s, Average:%.2e" % (name,average)
        ax.scatter(bins[:-1],counts,label=buf,marker=marker)
    plt.legend(loc='best', frameon="False", fontsize=9)
    plt.show()
                   



def histnoclose (data,name,numbins=200,marker='+',**kwargs):
    counts,bins = np.histogram(data,bins=numbins)
    average=np.average(data)
    buf = "%s , Average:%.2e $\AA$" % (name,average)
    fig,ax=plt.subplots()
    ax.set_ylabel("Counts")
    ax.set_xlabel(name)
    counts = zero_to_nan(counts)
    ax.scatter(bins[:-1],counts,label=buf,marker=marker)
    plt.legend(loc='upper left',frameon='False',fontsize=9)
