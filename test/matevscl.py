from os.path import getsize
#import pylab
#import matplotlib as mpl
import numpy as NP
from math import exp,log

#dC = NP.dtype([('re',float),('im',float)])

def load(fname,nb):
        sz=getsize(fname)
        l=int(sz/8)
	assert(nb*nb==l)
        fp=open(fname,'r')
        x=NP.fromfile(fname,dtype=float)
        x.shape=nb,nb
        return x

def evscl(x,escl,vscl,kbT):
	evs=0.0
	#print x.shape
	nb,nb2=x.shape
	assert (nb==nb2)
	for i in range(0,nb):
		ei=i*0.02-20.0
		for j in range(0,nb):
			vi=j*0.02-20.0
			if (x[i][j]>0.5):
				print ei,vi,ei*escl+vi*vscl,x[i][j]
				evs+=exp(-(ei*escl+vi*vscl)/kbT)*x[i][j]
	return evs

def getkBT(Temp):
	N_A = 6.022045000e+23;
	#//SI units, 2010 CODATA value, J/K = m2.kg/(s2.K) in SI base units
	#//const double k_B = 1.380648813e-23;
 	#//cal/K         1 steam table calorie = 4.1868 J
        k_B = 3.297623030e-24;
        return N_A*k_B*Temp/1000.0;

if __name__=="__main__":
	import sys
	fname=sys.argv[1]
	nb=int(sys.argv[2])
	m=load(fname,nb)
	escl=float(sys.argv[3])
	vscl=float(sys.argv[4])
	t=float(sys.argv[5])
	kbt=getkBT(t)
	print "#",escl,vscl,t,kbt
	bz=evscl(m,escl,vscl,kbt)
	print "# %16e"%bz	
