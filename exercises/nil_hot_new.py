import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from cycler import cycler
# rc('text',usetex=True)
plt.rc('axes', prop_cycle=(
    cycler(color=['k','g','b','r','c'])+
    cycler(linestyle=['-', '--', ':', '-.',(0, (3, 5, 1, 5, 1, 5))])))
fig, ax = plt.subplots(2,2,figsize=(9.7, 5.3))

#Decorations start here
plt.suptitle(r'$^{84}$Zr (neutrons only, at $\delta=0.2$)',fontsize=20)
plt.subplots_adjust(hspace=0,wspace=0.3)
                 
ax[0,0].set_xlim(15,80)
ax[0,0].set_ylim(0,1.1)
ax[0,0].set_ylabel(r'$n_i$')
ax[0,0].set_xticklabels([])
ax[0,0].axhline(y=0.5, color='k', linestyle='-',lw=1)

ax[1,0].set_xlim(15,80)
ax[1,0].set_ylim(0,0.75)
ax[1,0].set_ylabel(r'$s_i$')
ax[1,0].set_xlabel(r'$e_i$ (MeV)')

ax[0,1].set_ylabel(r'$\mu$ (MeV)')
ax[0,1].set_xlim(0.5,5.0)
ax[0,1].set_ylim(42.7,44.7)
ax[0,1].set_xticklabels([])

ax[1,1].set_xlabel('$T$ (MeV)')
ax[1,1].set_ylabel(r'$F=E-TS$ (MeV)')
ax[1,1].set_xlim(0.5,5.0)
ax[1,1].set_ylim(695,785)

for x in ax.flat:
  x.tick_params(direction='in',top='true',right='true')
  x.tick_params(which='both', width=1,labelsize=12)
  x.tick_params(which='major', length=7)
  x.minorticks_on()
  x.tick_params(which='minor', length=4, color='k', 
    direction='in',top='true',right='true')
  for axis in ['top','bottom','left','right']:
    x.spines[axis].set_linewidth(1)
  x.xaxis.label.set_fontsize(14)
  x.yaxis.label.set_fontsize(14)

#plotting starts here only
d = np.loadtxt('hot.out')
for t in np.arange(1.0,5.1,1.0):
  print(t)
  a=d[abs(d[:,0]-t)<1e-5,:]
  ax[0,0].plot(a[:,2],a[:,3],label=str(t),lw=3,zorder=1)
  ax[1,0].plot(a[:,2],a[:,4],lw=3,zorder=1)
a=d[abs(d[:,0]-22222)<1e-5,:]
ax[0,1].plot(a[:,1],a[:,2],lw=3,zorder=1)
ax[1,1].plot(a[:,1],a[:,3],lw=3,zorder=1)

#some more decorations
ax[0,0].legend(loc="upper right",frameon=False,fontsize=12,
    title=r'$T$ (MeV)',title_fontsize=14,handlelength=2.5)

plt.savefig('hot.pdf', bbox_inches = 'tight',
    pad_inches = 0)
plt.show()
