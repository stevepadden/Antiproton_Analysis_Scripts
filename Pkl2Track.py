#Converts a pickle file of particles that successfully traverse a foil into a beam readable within G4Beamline,
#Used to assess the ability of each foil within the particle trap itself
#To use this file call $python3 Pkl2Track.py <Input Data from MDRange (.pkl)> <Output data file name (.g4bl)> <Foil thickness (ANGSTROMS)>
import sys
import pandas as pd
from numpy import sqrt
"""
From MDrange document on http://beam.helsinki.fi/~knordlun/mdh/mdh_program.html
Internal quantity       Multiply with to obtain SI 
-----------------       --------------------------
       E                1.6022e-19 J
       r                1e-10 m
       F                1.6022e-9 N
       m                1.6606e-27*m(u) kg
       t                10.1805e-15*sqrt(m(u)) s
       v                9822.66/sqrt(m(u)) m/s
       a                9.6485e17/m(u) m/s^2

m(u) is atom mass in atomic mass units, 
"""

IncludeZ= False

data = pd.read_pickle(sys.argv[1])



data["vz"] = data["vz"].abs()
print(data)

pd.set_option("display.max_columns",None)
mu = 1
c = 299792458
def convert_dist(inp):
    return inp*1e-7

def convert_vel(inp):
    return (inp*9822.66)/sqrt(mu)

def convert_energy(inp):
    return inp*1.6022e-19

def convert_vel_2_P(inp):
    ms = inp*9822.66/sqrt(mu)
    mass = 938.272081358
    fracc = ms/c
    return mass*fracc

def transmit_energy(dataframe,Foil_Width):
    df = pd.DataFrame()
    for i,j in dataframe.iterrows():
        #print(i)
        
        range = dataframe.iloc[i].loc["range"]
        if range < float(Foil_Width):
            continue
        else:
            df = df.append(dataframe.iloc[i])
    return df

def use_range(dataframe,width):
    return dataframe[dataframe["range"]>width]



try :
    if sys.argv[3]:
        print("USING WIDTH %s"%(sys.argv[3]))
        x = int(sys.argv[3])
        data = data[data["sz"]>x]
except:
    print("No Width Input in argument 3")


data.reset_index(drop=True,inplace=True)


sx = convert_dist(data["sx"])

sy = convert_dist(data["sy"])
if IncludeZ == True:
    sz = convert_dist(data["sz"])
else:
    n = len(sx)
    sz= [0.0]* n
    
ranges = convert_dist(data["range"])
px = convert_vel_2_P(data["vx"])
py = convert_vel_2_P(data["vy"])
pz = convert_vel_2_P(data["vz"])
print(pz)
energy = convert_energy(data["energy"])

t = 0.0
pgid = -2212
eventid = 1
trackid = 1
parentid = 0
weight = 1



with open(sys.argv[2]+".txt",'wb') as f:
    f.write(b'#BLTrackFile inputbeam\r\n')
    f.write(b'#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\r\n')
    f.write(b'#mm mm mm MeV/c MeV/c MeV/c ns\r\n')
    for i in range(len(sx)):
        #print(i)
        f.write(b'%2.6f %2.6f %2.6f %2.6f %2.6f %2.6f %4.3f %i %i %i %i %i\r\n'% (sx[i],sy[i],sz[i],px[i],py[i],pz[i],t,pgid,eventid,trackid,parentid,weight))
        #f.write(b'%2.6f'%(sx[i]))
f.close()

fivekevdata = data[data["energy"]<5000]
fivekevdata.reset_index(inplace=True)

    
ranges = convert_dist(fivekevdata["range"])

px = convert_vel_2_P(fivekevdata["vx"])

py = convert_vel_2_P(fivekevdata["vy"])

pz = convert_vel_2_P(fivekevdata["vz"])

energy = convert_energy(fivekevdata["energy"])


sx = convert_dist(fivekevdata["sx"])

sy = convert_dist(fivekevdata["sy"])

if IncludeZ == True:
    sz = convert_dist(fivekevdata["sz"])
else:
    n = len(sx)
    sz= [0.0]* n

t = 0.0
pgid = -2212
eventid = 1
trackid = 1
parentid = 0
weight = 1


with open(sys.argv[2]+"_5kev.txt",'wb') as f:
    f.write(b'#BLTrackFile inputbeam\r\n')
    f.write(b'#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\r\n')
    f.write(b'#mm mm mm MeV/c MeV/c MeV/c ns\r\n')
    for i in range(len(sx)):
        #print(i)
        f.write(b'%2.6f %2.6f %2.6f %2.6f %2.6f %2.6f %4.3f %i %i %i %i %i\r\n'% (sx[i],sy[i],sz[i],px[i],py[i],pz[i],t,pgid,eventid,trackid,parentid,weight))
        #f.write(b'%2.6f'%(sx[i]))
f.close()
try:
    if sys.argv[3]:
        x=""
except:
        print("NO RANGE INPUT TO SYSTEM ARGUMENT 3")
