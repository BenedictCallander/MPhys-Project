import numpy as np 
import pandas as pd 
import matplotlib

import matplotlib.pyplot as plt 
import matplotlib.pylab as pyb
from matplotlib.gridspec import GridSpec 
from scipy.signal import medfilt, savgol_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as pat

dens = pd.read_csv("dens.csv")
met = pd.read_csv("met.csv")
sfr = pd.read_csv("sfr.csv")
hist = pd.read_csv("historgram.csv")
sfrprof = pd.read_csv("sfr2.csv")
star = pd.read_csv("starm.csv")
def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)
        


sprofile = savgol_filter(sfrprof['sfr'],window_length=5,polyorder=3)

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.style.use('dark_background')
ticks6 = np.arange(0,10.1,0.1)
fig = plt.figure(figsize=(30,15),constrained_layout = True,dpi=500)
gs = GridSpec(2,4, figure=fig)

ax1 = fig.add_subplot(gs[0,0]);ax2 = fig.add_subplot(gs[0,1]);ax3 = fig.add_subplot(gs[0,2])
ax4 = fig.add_subplot(gs[-1,:3]);ax5 = fig.add_subplot(gs[0,3]);ax6 = fig.add_subplot(gs[1,3])
divider1 = make_axes_locatable(ax1)  ;cax1 = divider1.append_axes('right', size='5%', pad=0.05)
divider2 = make_axes_locatable(ax2)  ;cax2 = divider2.append_axes('right', size='5%', pad=0.05)
divider3 = make_axes_locatable(ax3)  ;cax3 = divider3.append_axes('right', size='5%', pad=0.05)
divider5 = make_axes_locatable(ax5)  ;cax5 = divider5.append_axes('right', size='5%', pad=0.05)



ax1.set_title("$log_{10}$ Gas Column Density",fontsize=20); ax2.set_title("Subhalo Stellar Mass",fontsize=20); ax3.set_title("Subhalo Surface Metallicity",fontsize=20)
ax4.set_title("Subhalo Metallicity Profile",fontsize=20); ax5.set_title("Subhalo Star Formation Rate",fontsize=20);ax6.set_title("Subhalo Star Formation Rate Profile",fontsize=20)

ax1.set_xlabel('$\Delta x$ [kpc/h]',fontsize=20);ax1.set_ylabel('$\Delta y$ [kpc/h]',fontsize=20)
ax2.set_xlabel('$\Delta x$ [kpc/h]',fontsize=20);ax2.set_ylabel('$\Delta y$ [kpc/h]',fontsize=20)
ax3.set_xlabel('$\Delta x$ [kpc/h]',fontsize=20);ax3.set_ylabel('$\Delta y$ [kpc/h]',fontsize=20)
ax5.set_xlabel('$\Delta x$ [kpc/h]',fontsize=20);ax5.set_ylabel('$\Delta y$ [kpc/h]',fontsize=20)

#star minmax vmin=4.5, vmax = 7)
ax4.set_xlabel("Radial Distance from Galactic Centre",fontsize=20); ax4.set_ylabel(r"12+$log_{10}$ $\frac{M_Z}{M_{total}}$",fontsize=20)
ax6.set_xlabel("Radial Distance from Galactic Centre",fontsize=20);ax6.set_ylabel("Star Formation Rate [$M_\odot / yr$]",fontsize=20)


im1=ax1.scatter(dens['x'],dens['y'],c=(np.log10(dens['m'])),cmap='inferno', vmin=(min(np.log10(dens['m']))),vmax =(0.7*max(np.log10(dens['m']))))
im2=ax2.scatter(star['x'],star['y'],c =np.log10(star['m']), cmap = 'inferno', vmin=4.5, vmax=7 )
im3=ax3.scatter(met['x'],met['y'],c=(np.log10(met['met'])),cmap='inferno', vmin=(min(np.log10(met['met']))),vmax =(0.7*max(np.log10(met['met']))))
im4=ax4.hist2d(hist['rad'],12+np.log10(hist['met']),bins=[200,200], weights=1/hist['sfr'],cmap='PuOr')
im5=ax5.scatter(sfr['x'],sfr['y'],c=(np.log10(sfr['sfr'])),cmap='inferno')

im6=ax6.plot(ticks6,sprofile, 'r-' )

ax1.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
ax2.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
ax3.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
ax4.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
ax5.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)
ax6.tick_params(axis='both', which = 'both', direction='inout', length = 8, width =1)

fig.colorbar(im1, cax=cax1, orientation='vertical')
fig.colorbar(im2, cax=cax2, orientation='vertical')
fig.colorbar(im3, cax=cax3, orientation='vertical')
fig.colorbar(im5, cax=cax5, orientation='vertical')

an1 = pat.Annulus((0,0),13.55, width=1); dan1 = pat.Annulus((0,0),2*13.55, width=1); ax1.add_patch(an1);ax1.add_patch(dan1)
an2 = pat.Annulus((0,0),13.55, width=1); dan2 = pat.Annulus((0,0),2*13.55, width=1); ax2.add_patch(an2);ax2.add_patch(dan2)
an3 = pat.Annulus((0,0),13.55, width=1); dan3 = pat.Annulus((0,0),2*13.55, width=1); ax3.add_patch(an3);ax3.add_patch(dan3)
an5 = pat.Annulus((0,0),13.55, width=1); dan5 = pat.Annulus((0,0),2*13.55, width=1); ax5.add_patch(an5);ax5.add_patch(dan5)
fig.tight_layout()
fig.subplots_adjust(right=0.93)





fig.savefig("subtest1.png")



