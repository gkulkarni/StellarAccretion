import numpy as np
import matplotlib as mpl
mpl.use('Agg') 
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.stats import poisson_conf_interval as pci

fig = plt.figure(figsize =(6,6), dpi = 100)
gs = gridspec.GridSpec(4,1)
gs.update(hspace=0.0, wspace=0.0)

ax = plt.subplot(gs[0:3,:])

plt.ylabel(r'PDF')

ax.tick_params('both', which='major', length=5, width=1)
ax.tick_params('both', which='minor', length=2, width=1)
ax.tick_params('x', which='major', pad=6)

# feh, n = np.loadtxt('/data/shens/Eris_star_prop/Salvadori_2007_scaled.dat', unpack=True)
# plt.scatter(feh, n, s=32, edgecolor='none', c='dodgerblue', zorder=2, label='Salvadori et al.\ 2007')

zh, zl, n, n_orig = np.loadtxt('/data/shens/Eris_star_prop/Eris_mdf1_largebins.dat', unpack=True)

plt.step(zh, n_orig/0.2, c='forestgreen', lw=1.5, label='Eris without accretion', zorder=1)
plt.step(zh, n/0.2, c='red', lw=1.5, label='Eris with accretion', zorder=1)

# zh, zl, n, n_orig = np.loadtxt('/data/shens/Eris_star_prop/Eris_mdf1.dat', unpack=True)

# plt.step(zh, n_orig, c='black', lw=1.5, label='Eris without accretion $N_\mathrm{orig}$', zorder=1)
# # plt.step(zh, n, c='red', lw=1.5, label='Eris with accretion $N_\mathrm{acc}$', zorder=1)

ax.set_xticklabels('')

fedb, ndb, ndb_lerr, ndb_uerr = np.loadtxt('deBennasuti.dat', unpack=True)
totndb = np.sum(ndb)
pdb = ndb/totndb * 9686./937635
plt.scatter(fedb, pdb/0.1, s=32, edgecolor='none', c='k', zorder=3, label='de Bennasuti et al.\ 2016')

nlims = pci(ndb,interval='frequentist-confidence')
nlims = nlims/totndb * 9686./937635
uperr = (nlims[1] - pdb)/0.1
downerr = (pdb - nlims[0])/0.1
dpmine = uperr - downerr 

ndb_uerr = ndb_uerr/totndb * 9686./937635 / 0.1
ndb_lerr = -ndb_lerr/totndb * 9686./937635 / 0.1
dptheir = ndb_uerr + ndb_lerr

uerr = np.where(dpmine > dptheir, uperr, ndb_uerr)
lerr = np.where(dpmine > dptheir, downerr, ndb_lerr)

(_, caps, _) = plt.errorbar(fedb, pdb/0.1, ecolor='k', capsize=3, yerr=np.vstack((lerr, uerr)), fmt='None', zorder=3)

for cap in caps:
    cap.set_markeredgewidth(1)

plt.yscale('log')
plt.xlim(-7.5, -2)
plt.ylim(1.0e-6, 1.0)

ax.get_yticklabels()[-8].set_visible(False)

plt.legend(loc='upper left',fontsize=12,handlelength=3,frameon=False,framealpha=0.0,
           labelspacing=.1,handletextpad=0.4,borderpad=0.2,scatterpoints=1)


ax = plt.subplot(gs[3,0])

plt.ylabel(r'ratio')
plt.xlabel(r'$[\mathrm{Fe}/\mathrm{H}]$')

ax.tick_params('both', which='major', length=6, width=1)
ax.tick_params('x', which='major', pad=6)

zh, zl, n, n_orig = np.loadtxt('/data/shens/Eris_star_prop/Eris_mdf1_largebins.dat', unpack=True)
r = n/n_orig
# r2 = np.nan_to_num(r)
r2 = np.nan_to_num(np.where(r>1.0e10, 0.0, r))
plt.step(zh, r2, c='dodgerblue', lw=1.5, label='Eris without accretion', zorder=1)

plt.yticks((0,0.5,1,1.5))
plt.xlim(-7.5, -2)
plt.ylim(0.0, 1.5)

plt.savefig("mdf.pdf",bbox_inches='tight')
