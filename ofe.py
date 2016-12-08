import matplotlib as mpl
mpl.use('Agg') 
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(6, 6), dpi=100)
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel('[O/Fe]')
ax.set_xlabel('[Fe/H]')

ax.tick_params('both', which='major', length=7, width=1)
ax.tick_params('both', which='minor', length=3, width=1)
ax.tick_params('x', which='major', pad=6)

colors = [(0.0, 0.24705882352941178, 1.0),
          (1.0, 0.7686274509803922, 0.0),
          (0.9098039215686274, 0.0, 0.043137254901960784),
          (0.5411764705882353, 0.16862745098039217, 0.8862745098039215)]

with open('/data/shens/Eris_star_prop/feh_ofe_halo_agecut_binned.dat','r') as f:
    feh, ofe, ofe_16, ofe_84 = np.loadtxt(f, unpack=True)

ax.fill_between(feh[ofe>0.0], ofe_16[ofe>0.0], ofe_84[ofe>0.0], color=colors[0], alpha=0.3, edgecolor='none')
ax.plot(feh[ofe>0.0], ofe[ofe>0.0], c=colors[0], lw=1.5, label='halo (age-selected)')

with open('/data/shens/Eris_star_prop/feh_ofe_halo_zcut_binned.dat','r') as f:
    feh, ofe, ofe_16, ofe_84 = np.loadtxt(f, unpack=True)

ax.fill_between(feh[ofe>0.0], ofe_16[ofe>0.0], ofe_84[ofe>0.0], color=colors[1], alpha=0.3, edgecolor='none')
ax.plot(feh[ofe>0.0], ofe[ofe>0.0], c=colors[1], lw=1.5, label='halo (metallicity-selected)')

with open('/data/shens/Eris_star_prop/feh_ofe_bulge_agecut_binned.dat','r') as f:
    feh, ofe, ofe_16, ofe_84 = np.loadtxt(f, unpack=True)

ax.fill_between(feh[ofe>0.0], ofe_16[ofe>0.0], ofe_84[ofe>0.0], color=colors[3], alpha=0.3, edgecolor='none')
ax.plot(feh[ofe>0.0], ofe[ofe>0.0], c=colors[3], lw=1.5, label='bulge (age-selected)')

with open('/data/shens/Eris_star_prop/feh_ofe_bulge_zcut_binned.dat','r') as f:
    feh, ofe, ofe_16, ofe_84 = np.loadtxt(f, unpack=True)

ax.fill_between(feh[ofe>0.0], ofe_16[ofe>0.0], ofe_84[ofe>0.0], color=colors[2], alpha=0.3, edgecolor='none')
ax.plot(feh[ofe>0.0], ofe[ofe>0.0], c=colors[2], lw=1.5, label='bulge (metallicity-selected)')

with open('stars.dat','r') as f:
    feh, ofe = np.loadtxt(f, usecols=(2,3), unpack=True)
ax.scatter(feh, ofe, s=24, edgecolors='k', c='#ffffff', zorder=5)

with open('MgFe_only.dat','r') as f:
    feh, ofe = np.loadtxt(f, usecols=(2,3), unpack=True)
ax.scatter(feh, ofe, s=24, edgecolors='k', c='#ffffff', zorder=5)

ax.set_ylabel('[O/Fe]')
ax.set_xlabel('[Fe/H]')

plt.xlim(-10, 1)
plt.ylim(-0.2, 0.71)

plt.legend(loc='lower left',fontsize=12,handlelength=3,frameon=False,framealpha=0.0,
           labelspacing=.1,handletextpad=0.4,borderpad=0.2,scatterpoints=1)

plt.savefig("ofe.pdf",bbox_inches='tight')


