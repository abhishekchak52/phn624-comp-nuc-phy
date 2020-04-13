import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
# rc('text',usetex=True)
d = np.loadtxt('hot.out')
pl.figure(1)
ni=pl.subplot(221)
si=pl.subplot(223)
for t in np.arange(1.0,5.1,1.0):
  print(t)
  a=d[abs(d[:,0]-t)<1e-5,:]
  ni.plot(a[:,2],a[:,3],label='T = '+str(t)+' MeV')
  si.plot(a[:,2],a[:,4])
a=d[abs(d[:,0]-22222)<1e-5,:]
mu=pl.subplot(222)
mu.plot(a[:,1],a[:,2])
fe=pl.subplot(224)
fe.plot(a[:,1],a[:,3])

#Decorations start here
pl.suptitle(r'$^{84}$Zr (neutrons only at $\delta=0.2$)')
pl.subplots_adjust(hspace=0)

ni.legend(loc="upper right")
ni.set_xlim(15,70)
ni.set_ylim(0,1.1)
ni.set_ylabel(r'$n_i$')
ni.set_xticks([])
ni.axhline(y=0.5, color='r', linestyle='-')

si.set_xlim(15,70)
si.set_ylim(0,0.75)
si.set_ylabel(r'$s_i$')
si.set_xlabel(r'$e_i$')

mu.set_ylabel(r'$\mu$ (MeV)')
mu.set_xlim(0.5,5.0)
mu.set_ylim(42.7,44.7)
mu.set_xticks([])

fe.set_xlabel('$T$ (MeV)')
fe.set_ylabel(r'$F=E-TS$ (MeV)')
fe.set_xlim(0.5,5.0)
fe.set_ylim(695,785)
pl.savefig('nil_out.png')
