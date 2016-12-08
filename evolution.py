import numpy as np
import matplotlib as mpl
mpl.use('Agg') 
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '14'
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
# zlims = (1.0,0.0)
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

def plot_gas_density(fig):

    mpl.rcParams['font.size'] = '14'

    ax = fig.add_subplot(nplots_x, nplots_y, plot_number+1)
    ax.set_xlim(zmin, zmax)
    ax.set_xscale('log')
    
    ax.tick_params('both', which='major', length=5, width=1)    

    with open('accretion_rate_halo_agecut.dat','r') as f:
        z, rhog_median, rhog_16, rhog_84 = np.loadtxt(f, usecols=(0,5,6,7), unpack=True)
    ax.fill_between(1.+z, np.log10(rhog_16), np.log10(rhog_84), color=colors[0], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(rhog_median), color=colors[0], lw=1.5)

    with open('accretion_rate_halo_zcut.dat','r') as f:
        z, rhog_median, rhog_16, rhog_84 = np.loadtxt(f, usecols=(0,5,6,7), unpack=True)
    ax.fill_between(1.+z, np.log10(rhog_16), np.log10(rhog_84), color=colors[1], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(rhog_median), color=colors[1], lw=1.5)

    with open('accretion_rate_bulge_zcut.dat','r') as f:
        z, rhog_median, rhog_16, rhog_84 = np.loadtxt(f, usecols=(0,5,6,7), unpack=True)
    ax.fill_between(1.+z, np.log10(rhog_16), np.log10(rhog_84), color=colors[2], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(rhog_median), color=colors[2], lw=1.5)

    with open('accretion_rate_bulge_agecut.dat','r') as f:
        z, rhog_median, rhog_16, rhog_84 = np.loadtxt(f, usecols=(0,5,6,7), unpack=True)
    ax.fill_between(1.+z, np.log10(rhog_16), np.log10(rhog_84), color=colors[3], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(rhog_median), color=colors[3], lw=1.5)

    redshift, rvir, nH_rvir, nH_0p5pkpc, nH_2pkpc, nH_25pkpc = np.loadtxt('Eris_gas_density.dat', unpack=True)
    ax.plot(1.+redshift, np.log10(nH_rvir), c='k', lw=2, alpha=1.0, label='$n (r<r_\mathrm{vir})$')
    ax.plot(1.+redshift, np.log10(nH_0p5pkpc), c='maroon', lw=2, alpha=1.0, label='$n (r<0.5~\mathrm{pkpc})$')
    ax.plot(1.+redshift, np.log10(nH_2pkpc), c='maroon', lw=2, dashes=[7,2], alpha=1.0, label='$n (r<2~\mathrm{pkpc})$')
    ax.plot(1.+redshift, np.log10(nH_25pkpc), c='maroon', lw=2, dashes=[7,2,2,2], alpha=1.0, label='$n (r<25~\mathrm{pkpc})$')
    
    ax.set_ylabel(r'$\log_{10}\left(n_\mathrm{H}/\mathrm{cm}^{-3}\right)$', fontsize=14)
    ax.yaxis.labelpad = 8
    ax.set_xticklabels('')

    plt.ylim(-5,4)
    
    plt.text(0.02, 0.1, '(a)', horizontalalignment='left', transform=ax.transAxes, fontsize='large')

    plt.legend(loc='upper right', fontsize=10, handlelength=3,
               frameon=False, framealpha=0.0, labelspacing=.1,
               handletextpad=0.4, borderpad=0.2, ncol=2, columnspacing=0)
    
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
    ax.set_xticks(zs)
    plt.minorticks_off()
    
    return

def plot_radial_distance(fig):

    mpl.rcParams['font.size'] = '14'

    ax = fig.add_subplot(nplots_x, nplots_y, plot_number+2)
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_label_position('right')
    ax.set_xscale('log')
    ax.set_xlim(zmin, zmax)
    
    ax.tick_params('both', which='major', length=5, width=1)

    with open('distance_history_median_halo_agecut.dat','r') as f:
        z, r_median, r_16, r_84 = np.loadtxt(f, usecols=(0,5,6,7), unpack=True)
    print z.max(), z.min()
    ax.fill_between(1.+z, np.log10(r_16), np.log10(r_84), color=colors[0], alpha=0.3, edgecolor='none') 
    ax.plot(1.+z, np.log10(0.64*r_median), color=colors[0], lw=1.5)

    with open('distance_history_median_halo_zcut.dat','r') as f:
        z, r_median, r_16, r_84 = np.loadtxt(f, usecols=(0,5,6,7), unpack=True)
    ax.fill_between(1.+z, np.log10(r_16), np.log10(r_84), color=colors[1], alpha=0.3, edgecolor='none') 
    ax.plot(1.+z, np.log10(0.64*r_median), color=colors[1], lw=1.5)

    with open('distance_history_median_bulge_zcut.dat','r') as f:
        z, r_median, r_16, r_84 = np.loadtxt(f, usecols=(0,5,6,7), unpack=True)
    ax.fill_between(1.+z, np.log10(r_16), np.log10(r_84), color=colors[2], alpha=0.3, edgecolor='none') 
    ax.plot(1.+z, np.log10(0.64*r_median), color=colors[2], lw=1.5)

    with open('distance_history_median_bulge_agecut.dat','r') as f:
        z, r_median, r_16, r_84 = np.loadtxt(f, usecols=(0,5,6,7), unpack=True)
    ax.fill_between(1.+z, np.log10(r_16), np.log10(r_84), color=colors[3], alpha=0.3, edgecolor='none') 
    ax.plot(1.+z, np.log10(0.64*r_median), color=colors[3], lw=1.5)

#    ax.set_ylim(-26, -16)
    
    ax.set_ylabel(r'$\log_{10}\left(r/r_\mathrm{vir}\right)$', fontsize=14)

    plt.axhline(0.0, c='k', lw=2, dashes=[7,2], alpha=1.0)

    ax.set_xticklabels('')

    plt.text(0.02, 0.1, '(b)', horizontalalignment='left', transform=ax.transAxes, fontsize='large')
    plt.text(2.0, 0.05, '$r=r_\mathrm{vir}$', fontsize=12)
    
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
    ax.set_xticks(zs)
    plt.minorticks_off()

    return

def plot_sound_speed(fig):

    mpl.rcParams['font.size'] = '14'

    ax = fig.add_subplot(nplots_x, nplots_y, plot_number+4)
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_label_position('right')
    ax.set_xscale('log')
    ax.set_xlim(zmin, zmax)
    
    ax.tick_params('both', which='major', length=5, width=1)

    with open('accretion_rate_halo_agecut.dat','r') as f:
        z, cs_median, cs_16, cs_84 = np.loadtxt(f, usecols=(0,20,21,22), unpack=True)
    ax.fill_between(1.+z, np.log10(cs_16), np.log10(cs_84), color=colors[0], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(cs_median), color=colors[0], lw=1.5, label='halo (age)')

    with open('accretion_rate_halo_zcut.dat','r') as f:
        z, cs_median, cs_16, cs_84 = np.loadtxt(f, usecols=(0,20,21,22), unpack=True)
    ax.fill_between(1.+z, np.log10(cs_16), np.log10(cs_84), color=colors[1], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(cs_median), color=colors[1], lw=1.5, label='halo (metallicity)')

    with open('accretion_rate_bulge_zcut.dat','r') as f:
        z, cs_median, cs_16, cs_84 = np.loadtxt(f, usecols=(0,20,21,22), unpack=True)
    ax.fill_between(1.+z, np.log10(cs_16), np.log10(cs_84), color=colors[2], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(cs_median), color=colors[2], lw=1.5, label='bulge (metallicity)')

    with open('accretion_rate_bulge_agecut.dat','r') as f:
        z, cs_median, cs_16, cs_84 = np.loadtxt(f, usecols=(0,20,21,22), unpack=True)
    ax.fill_between(1.+z, np.log10(cs_16), np.log10(cs_84), color=colors[3], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(cs_median), color=colors[3], lw=1.5, label='bulge (age)')
    
    ax.set_ylabel(r'$\log_{10}\left(c_s/\mathrm{km}\,\mathrm{s}^{-1}\right)$', fontsize=14)
    ax.yaxis.labelpad = 16

    plt.axhline(0.0, c='k', lw=1, dashes=[7,2], alpha=1.0)

    ax.set_xticklabels('')
    plt.ylim(0.5, 3.0)
    
    plt.text(0.02, 0.1, '(d)', horizontalalignment='left', transform=ax.transAxes, fontsize='large')
    
    # ages = np.array([13, 12, 10, 5, 1])*u.Gyr
    # ageticks = np.array([z_at_value(cosmo.lookback_time, age) for age in ages])
    # ageticks = (1.0+ageticks)
    # print ageticks, ages
    # ax2 = ax.twiny()
    # ax2.set_xscale('log')
    # plt.minorticks_off()
    # ax2.set_xlim(zmin, zmax)
    # ax2.tick_params('both', which='major', length=5, width=1)
    # ax2.set_xticks(ageticks)
    # ax2.set_xticklabels(['{:g}'.format(age) for age in ages.value])
    # ax2.set_xlabel(r'lookback time [Gyr]')

    zs = np.array([1,2,3,4,5,6,7,8,9,10,12])
    ax.set_xticks(zs)
    plt.minorticks_off()

    return


def plot_velocity(fig):

    mpl.rcParams['font.size'] = '14'

    ax = fig.add_subplot(nplots_x, nplots_y, plot_number+3)
    ax.set_xlim(zmin, zmax)
    ax.set_xscale('log')
    
    ax.tick_params('both', which='major', length=5, width=1)    
    
    with open('accretion_rate_halo_agecut.dat','r') as f:
        z, deltav_median, deltav_16, deltav_84 = np.loadtxt(f, usecols=(0,8,9,10), unpack=True)
    ax.fill_between(1.+z, np.log10(deltav_16), np.log10(deltav_84), color=colors[0], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(deltav_median), color=colors[0], lw=1.5)

    with open('accretion_rate_halo_zcut.dat','r') as f:
        z, deltav_median, deltav_16, deltav_84 = np.loadtxt(f, usecols=(0,8,9,10), unpack=True)
    ax.fill_between(1.+z, np.log10(deltav_16), np.log10(deltav_84), color=colors[1], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(deltav_median), color=colors[1], lw=1.5)

    with open('accretion_rate_bulge_zcut.dat','r') as f:
        z, deltav_median, deltav_16, deltav_84 = np.loadtxt(f, usecols=(0,8,9,10), unpack=True)
    ax.fill_between(1.+z, np.log10(deltav_16), np.log10(deltav_84), color=colors[2], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(deltav_median), color=colors[2], lw=1.5)

    with open('accretion_rate_bulge_agecut.dat','r') as f:
        z, deltav_median, deltav_16, deltav_84 = np.loadtxt(f, usecols=(0,8,9,10), unpack=True)
    ax.fill_between(1.+z, np.log10(deltav_16), np.log10(deltav_84), color=colors[3], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, np.log10(deltav_median), color=colors[3], lw=1.5)

    zv, tv, vmax, vcirc = np.loadtxt('Eris_vmax.dat', unpack=True)
    ax.plot(1.+zv, np.log10(vmax), c='k', lw=2, dashes=[7,2], alpha=1.0, label='$v_\mathrm{max}$')
    ax.plot(1.+zv, np.log10(vcirc), c='maroon', lw=2, dashes=[7,2], alpha=1.0, label='$v_\mathrm{circ}$')
    
    ax.set_ylabel(r'$\log_{10}\left(v_\mathrm{rel}/\mathrm{km}\,\mathrm{s}^{-1}\right)$', fontsize=14)
    ax.yaxis.labelpad = 9

    ax.set_xticklabels('')
    plt.ylim(0.5, 3.0)

    plt.text(0.02, 0.1, '(c)', horizontalalignment='left', transform=ax.transAxes, fontsize='large')

    plt.legend(loc='lower right', fontsize=10, handlelength=3,
               frameon=False, framealpha=0.0, labelspacing=.1,
               handletextpad=0.4, borderpad=0.2)
    
    zs = np.array([1,2,3,4,5,6,7,8,9,10,12])
    ax.set_xticks(zs)
    # ax.set_xticklabels(['{:g}'.format(z-1) for z in zs])
    plt.minorticks_off()

    return

def plot_gas_metallicity(fig):

    mpl.rcParams['font.size'] = '14'

    ax = fig.add_subplot(nplots_x, nplots_y, plot_number+5)
    ax.set_xlim(zmin, zmax)
    ax.set_xscale('log')

    ax.tick_params('both', which='major', length=5, width=1)    

    with open('accretion_rate_halo_agecut.dat','r') as f:
        z, fefrac_median, fefrac_16, fefrac_84 = np.loadtxt(f, usecols=(0,11,12,13), unpack=True)
    ax.fill_between(1.+z, fefrac_16, fefrac_84, color=colors[0], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, fefrac_median, color=colors[0], lw=1.5)

    with open('accretion_rate_halo_zcut.dat','r') as f:
        z, fefrac_median, fefrac_16, fefrac_84 = np.loadtxt(f, usecols=(0,11,12,13), unpack=True)
    ax.fill_between(1.+z, fefrac_16, fefrac_84, color=colors[1], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, fefrac_median, color=colors[1], lw=1.5)

    with open('accretion_rate_bulge_zcut.dat','r') as f:
        z, fefrac_median, fefrac_16, fefrac_84 = np.loadtxt(f, usecols=(0,11,12,13), unpack=True)
    ax.fill_between(1.+z, fefrac_16, fefrac_84, color=colors[2], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, fefrac_median, color=colors[2], lw=1.5)

    with open('accretion_rate_bulge_agecut.dat','r') as f:
        z, fefrac_median, fefrac_16, fefrac_84 = np.loadtxt(f, usecols=(0,11,12,13), unpack=True)
    ax.fill_between(1.+z, fefrac_16, fefrac_84, color=colors[3], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, fefrac_median, color=colors[3], lw=1.5)
    
    ax.set_ylabel(r'[Fe/H]$_\mathrm{gas}$', fontsize=14)
    ax.yaxis.labelpad = 6
    ax.set_xlabel('redshift')

    ax.set_xticklabels('')

    plt.text(0.02, 0.1, '(e)', horizontalalignment='left', transform=ax.transAxes, fontsize='large')

    zs = np.array([1,2,3,4,5,6,7,8,9,10,12])
    ax.set_xticks(zs)
    ax.set_xticklabels(['{:g}'.format(z-1) for z in zs])
    plt.minorticks_off()

    return

def plot_gas_obyfe(fig):

    mpl.rcParams['font.size'] = '14'

    ax = fig.add_subplot(nplots_x, nplots_y, plot_number+5)
    ax.set_xlim(zmin, zmax)
    ax.set_xscale('log')

    ax.tick_params('both', which='major', length=5, width=1)    

    with open('accretion_rate_halo_agecut.dat','r') as f:
        z, obyfe_median, obyfe_16, obyfe_84 = np.loadtxt(f, usecols=(0,17,18,19), unpack=True)
    ax.fill_between(1.+z, obyfe_16, obyfe_84, color=colors[0], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, obyfe_median, color=colors[0], lw=1.5, label='halo (age)')

    with open('accretion_rate_halo_zcut.dat','r') as f:
        z, obyfe_median, obyfe_16, obyfe_84 = np.loadtxt(f, usecols=(0,17,18,19), unpack=True)
    ax.fill_between(1.+z, obyfe_16, obyfe_84, color=colors[1], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, obyfe_median, color=colors[1], lw=1.5, label='halo (metallicity)')

    with open('accretion_rate_bulge_zcut.dat','r') as f:
        z, obyfe_median, obyfe_16, obyfe_84 = np.loadtxt(f, usecols=(0,17,18,19), unpack=True)
    ax.fill_between(1.+z, obyfe_16, obyfe_84, color=colors[2], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, obyfe_median, color=colors[2], lw=1.5, label='bulge (metallicity)')
        
    with open('accretion_rate_bulge_agecut.dat','r') as f:
        z, obyfe_median, obyfe_16, obyfe_84 = np.loadtxt(f, usecols=(0,17,18,19), unpack=True)
    ax.fill_between(1.+z, obyfe_16, obyfe_84, color=colors[3], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, obyfe_median, color=colors[3], lw=1.5, label='bulge (age)')
        
    ax.set_ylabel(r'[O/Fe]$_\mathrm{gas}$', fontsize=14)
    ax.yaxis.labelpad = -2
    ax.set_xlabel('redshift')

    plt.text(0.02, 0.1, '(e)', horizontalalignment='left', transform=ax.transAxes, fontsize='large')

    zs = np.array([1,2,3,4,5,6,7,8,9,10,12])
    ax.set_xticks(zs)
    ax.set_xticklabels(['{:g}'.format(z-1) for z in zs])
    plt.minorticks_off()

    return 


def plot_gas_oxygen(fig):

    mpl.rcParams['font.size'] = '14'

    ax = fig.add_subplot(nplots_x, nplots_y, plot_number+6)
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_label_position('right')
    ax.set_xlim(zmin, zmax)
    ax.set_xscale('log')

    ax.tick_params('both', which='major', length=5, width=1)    

    with open('accretion_rate_halo_agecut.dat','r') as f:
        z, ofrac_median, ofrac_16, ofrac_84 = np.loadtxt(f, usecols=(0,14,15,16), unpack=True)
    ax.fill_between(1.+z, ofrac_16, ofrac_84, color=colors[0], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, ofrac_median, color=colors[0], lw=1.5, label='halo (age-selected)')

    with open('accretion_rate_halo_zcut.dat','r') as f:
        z, ofrac_median, ofrac_16, ofrac_84 = np.loadtxt(f, usecols=(0,14,15,16), unpack=True)
    ax.fill_between(1.+z, ofrac_16, ofrac_84, color=colors[1], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, ofrac_median, color=colors[1], lw=1.5, label='halo (metallicity-selected)')

    with open('accretion_rate_bulge_zcut.dat','r') as f:
        z, ofrac_median, ofrac_16, ofrac_84 = np.loadtxt(f, usecols=(0,14,15,16), unpack=True)
    ax.fill_between(1.+z, ofrac_16, ofrac_84, color=colors[2], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, ofrac_median, color=colors[2], lw=1.5, label='bulge (metallicity-selected)')

    with open('accretion_rate_bulge_agecut.dat','r') as f:
        z, ofrac_median, ofrac_16, ofrac_84 = np.loadtxt(f, usecols=(0,14,15,16), unpack=True)
    ax.fill_between(1.+z, ofrac_16, ofrac_84, color=colors[3], alpha=0.3, edgecolor='none')
    ax.plot(1.+z, ofrac_median, color=colors[3], lw=1.5, label='bulge (age-selected)')
    
    ax.set_ylabel(r'[O/H]$_\mathrm{gas}$', fontsize=14)
    ax.set_xlabel('redshift')
    ax.set_ylim(-3.0,0.5)

    plt.text(0.02, 0.1, '(f)', horizontalalignment='left', transform=ax.transAxes, fontsize='large')

    zs = np.array([1,2,3,4,5,6,7,8,9,10,12])
    ax.set_xticks(zs)
    ax.set_xticklabels(['{:g}'.format(z-1) for z in zs])
    plt.minorticks_off()

    plt.legend(loc='lower right', fontsize=10, handlelength=3,
               frameon=False, framealpha=0.0, labelspacing=.1,
               handletextpad=0.4, borderpad=0.2, ncol=2, columnspacing=0)
    

    return 

def summary_plot():

    nx = 3
    ny = 2 
    factor_x = 3.0
    factor_y = 3.0
    ldim = 0.25*factor_x
    bdim = 0.25*factor_y
    rdim = 0.25*factor_x
    tdim = 0.3*factor_y
    wspace = 0.5
    hspace = 0.9 

    plotdim_x = factor_x*nx + (nx-1)*wspace
    plotdim_y = factor_y*ny + (ny-1)*hspace

    hdim = plotdim_x + ldim + rdim 
    vdim = plotdim_y + tdim + bdim 

    fig = plt.figure(figsize=(hdim, vdim), dpi=100)

    l = ldim/hdim
    b = bdim/vdim
    r = (ldim + plotdim_x)/hdim
    t = (bdim + plotdim_y)/vdim 
    fig.subplots_adjust(left=l, bottom=b, right=r, top=t, wspace=wspace/hdim, hspace=hspace/vdim)

    plot_gas_density(fig)
    plot_radial_distance(fig)
    plot_velocity(fig)
    plot_sound_speed(fig)
    plot_gas_metallicity(fig)
    plot_gas_oxygen(fig)

    plt.savefig('evolution2.pdf')

    return

summary_plot()

