import matplotlib as mpl
mpl.use('Agg') 
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy as sp

def fill_between_steps(ax, x, y1, y2=0, step_where='pre', **kwargs):
    ''' fill between a step plot and 

    Parameters
    ----------
    ax : Axes
       The axes to draw to

    x : array-like
        Array/vector of index values.

    y1 : array-like or float
        Array/vector of values to be filled under.
    y2 : array-Like or float, optional
        Array/vector or bottom values for filled area. Default is 0.

    step_where : {'pre', 'post', 'mid'}
        where the step happens, same meanings as for `step`

    **kwargs will be passed to the matplotlib fill_between() function.

    Returns
    -------
    ret : PolyCollection
       The added artist

    https://github.com/matplotlib/matplotlib/issues/643#issuecomment-63762568

    '''
    if step_where not in {'pre', 'post', 'mid'}:
        raise ValueError("where must be one of {{'pre', 'post', 'mid'}} "
                         "You passed in {wh}".format(wh=step_where))

    # make sure y values are up-converted to arrays 
    if np.isscalar(y1):
        y1 = np.ones_like(x) * y1

    if np.isscalar(y2):
        y2 = np.ones_like(x) * y2

    # temporary array for up-converting the values to step corners
    # 3 x 2N - 1 array 

    vertices = np.vstack((x, y1, y2))

    # this logic is lifted from lines.py
    # this should probably be centralized someplace
    if step_where == 'pre':
        steps = np.ma.zeros((3, 2 * len(x) - 1), np.float)
        steps[0, 0::2], steps[0, 1::2] = vertices[0, :], vertices[0, :-1]
        steps[1:, 0::2], steps[1:, 1:-1:2] = vertices[1:, :], vertices[1:, 1:]

    elif step_where == 'post':
        steps = np.ma.zeros((3, 2 * len(x) - 1), np.float)
        steps[0, ::2], steps[0, 1:-1:2] = vertices[0, :], vertices[0, 1:]
        steps[1:, 0::2], steps[1:, 1::2] = vertices[1:, :], vertices[1:, :-1]

    elif step_where == 'mid':
        steps = np.ma.zeros((3, 2 * len(x)), np.float)
        steps[0, 1:-1:2] = 0.5 * (vertices[0, :-1] + vertices[0, 1:])
        steps[0, 2::2] = 0.5 * (vertices[0, :-1] + vertices[0, 1:])
        steps[0, 0] = vertices[0, 0]
        steps[0, -1] = vertices[0, -1]
        steps[1:, 0::2], steps[1:, 1::2] = vertices[1:, :], vertices[1:, :]
    else:
        raise RuntimeError("should never hit end of if-elif block for validated input")

    # un-pack
    xx, yy1, yy2 = steps

    # now to the plotting part:
    return ax.fill_between(xx, yy1, y2=yy2, **kwargs)

cols = {'halo_agecut' : (0.0, 0.24705882352941178, 1.0),
        'halo_zcut' : (1.0, 0.7686274509803922, 0.0),
        'bulge_zcut' : (0.9098039215686274, 0.0, 0.043137254901960784),
        'bulge_agecut' : (0.5411764705882353, 0.16862745098039217, 0.8862745098039215)}

# fig = plt.figure(figsize=(6, 6/1.62), dpi=100)
fig = plt.figure(figsize=(6*1.62, 6), dpi=100)
ax = fig.add_subplot(1, 1, 1)
plt.ylabel(r'PDF')
plt.xlabel(r'$[\mathrm{Fe}/\mathrm{H}]_\mathrm{accreted}$')

ax.tick_params('both', which='major', length=7, width=1)
ax.tick_params('both', which='minor', length=3, width=1)
ax.tick_params('x', which='major', pad=6)

zh, zl, n = np.loadtxt('/data/shens/Eris_star_prop/FeH_histogram_halo_zcut_format.dat', unpack=True)
lines = plt.step(zh+np.log10(0.64/3.0), n/n.sum()/0.1, c=cols['halo_zcut'], lw=1.5, label='halo (metallicity-selected)')
fill_between_steps(ax, zh+np.log10(0.64/3.0), n/n.sum()/0.1, 0, color=cols['halo_zcut'], alpha=0.3)

zh, zl, n = np.loadtxt('/data/shens/Eris_star_prop/FeH_histogram_halo_agecut_format.dat', unpack=True)
lines = plt.step(zh+np.log10(0.64/3.0), n/n.sum()/0.1, c=cols['halo_agecut'], lw=1.5, label='halo (age-selected)')
fill_between_steps(ax, zh+np.log10(0.64/3.0), n/n.sum()/0.1, 0, color=cols['halo_agecut'], alpha=0.3)

zh, zl, n = np.loadtxt('/data/shens/Eris_star_prop/FeH_histogram_bulge_zcut_format.dat', unpack=True)
lines = plt.step(zh+np.log10(0.64/3.0), n/n.sum()/0.1, c=cols['bulge_zcut'], lw=1.5, label='bulge (metallicity-selected)')
fill_between_steps(ax, zh+np.log10(0.64/3.0), n/n.sum()/0.1, 0, color=cols['bulge_zcut'], alpha=0.3)

zh, zl, n = np.loadtxt('/data/shens/Eris_star_prop/FeH_histogram_bulge_agecut_format.dat', unpack=True)
lines = plt.step(zh+np.log10(0.64/3.0), n/n.sum()/0.1, c=cols['bulge_agecut'], lw=1.5, label='bulge (age-selected)')
fill_between_steps(ax, zh+np.log10(0.64/3.0), n/n.sum()/0.1, 0, color=cols['bulge_agecut'], alpha=0.3)

plt.xlim(-11, -1)

plt.legend(loc='upper left',fontsize=14,handlelength=3,frameon=False,framealpha=0.0,
           labelspacing=.1,handletextpad=0.4,borderpad=0.2)



plt.savefig("acchist.pdf",bbox_inches='tight')


