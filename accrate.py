import numpy as np
import matplotlib as mpl
mpl.use('Agg') 
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import z_at_value 
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a'] 
nplots_x = 3
nplots_y = 2
nplots = 5
plot_number = 0 

zlims=(13.0, 1.0)
zmin, zmax = zlims
z = np.linspace(zmin, zmax, num=50)

colors = [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
          (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
          (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
          (0.5058823529411764, 0.4470588235294118, 0.6980392156862745),
          (0.8, 0.7254901960784313, 0.4549019607843137),
          (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]

colors = [(0.0, 0.24705882352941178, 1.0),
          (1.0, 0.7686274509803922, 0.0),
          (0.9098039215686274, 0.0, 0.043137254901960784),
          (0.5411764705882353, 0.16862745098039217, 0.8862745098039215)]

fig = plt.figure(figsize=(6, 6), dpi=100)

ax = fig.add_subplot(1, 1, 1)
ax.set_xscale('log')
ax.set_xlim(zmin, zmax)

ax.tick_params('both', which='major', length=5, width=1)

with open('accretion_rate_halo_zcut.dat','r') as f:
    z, dmfedt_median, dmfedt_16, dmfedt_84  = np.loadtxt(f, usecols=(0,2,3,4), unpack=True)
ax.fill_between(1.+z, np.log10(dmfedt_16), np.log10(dmfedt_84), color=colors[1], alpha=0.3, edgecolor='none') 
ax.plot(1.+z, np.log10(0.64*dmfedt_median), color=colors[1], lw=1.5, label='halo (metallicity-selected)')

with open('accretion_rate_halo_agecut.dat','r') as f:
    z, dmfedt_median, dmfedt_16, dmfedt_84  = np.loadtxt(f, usecols=(0,2,3,4), unpack=True)
print z.max(), z.min()
ax.fill_between(1.+z, np.log10(dmfedt_16), np.log10(dmfedt_84), color=colors[0], alpha=0.3, edgecolor='none') 
ax.plot(1.+z, np.log10(0.64*dmfedt_median), color=colors[0], lw=1.5, label='halo (age-selected)')


with open('accretion_rate_bulge_zcut.dat','r') as f:
    z, dmfedt_median, dmfedt_16, dmfedt_84  = np.loadtxt(f, usecols=(0,2,3,4), unpack=True)
ax.fill_between(1.+z, np.log10(dmfedt_16), np.log10(dmfedt_84), color=colors[2], alpha=0.3, edgecolor='none') 
ax.plot(1.+z, np.log10(0.64*dmfedt_median), color=colors[2], lw=1.5, label='bulge (metallicity-selected)')

with open('accretion_rate_bulge_agecut.dat','r') as f:
    z, dmfedt_median, dmfedt_16, dmfedt_84  = np.loadtxt(f, usecols=(0,2,3,4), unpack=True)
ax.fill_between(1.+z, np.log10(dmfedt_16), np.log10(dmfedt_84), color=colors[3], alpha=0.3, edgecolor='none') 
ax.plot(1.+z, np.log10(0.64*dmfedt_median), color=colors[3], lw=1.5, label='bulge (age-selected)')

ax.set_ylim(-26, -17)

ax.set_ylabel(r'$\log_{10}\left(\dot m_\mathrm{Fe}/\mathrm{M}_\odot\mathrm{yr}^{-1}\right)$')

ax.set_xticklabels('')

plt.yticks(np.arange(-26,-16,2))

plt.legend(loc='upper right',fontsize=12,handlelength=3,frameon=False,framealpha=0.0,
           labelspacing=.1,handletextpad=0.4,borderpad=0.2,scatterpoints=1)

ages = np.array([13, 12, 10, 5, 1])*u.Gyr
ageticks = np.array([z_at_value(cosmo.lookback_time, age) for age in ages])
ageticks = (1.0+ageticks)
print ageticks, ages
ax2 = ax.twiny()
ax2.set_xscale('log')
plt.minorticks_off()
ax2.set_xlim(zmin, zmax)
ax2.tick_params('both', which='major', length=5, width=1)
ax2.set_xticks(ageticks)
ax2.set_xticklabels(['{:g}'.format(age) for age in ages.value])
ax2.set_xlabel(r'lookback time [Gyr]')

zs = np.array([1,2,3,4,5,6,7,8,9,10,12])
ax.set_xticklabels(['{:g}'.format(z-1) for z in zs])
ax.set_xticks(zs)
plt.minorticks_off()
ax.set_xlabel('redshift')

plt.savefig('accrate.pdf', bbox_inches='tight')


