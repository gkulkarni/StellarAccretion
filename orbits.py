import numpy as np
import matplotlib as mpl
mpl.use('Agg') 
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage.filters import gaussian_filter as gf

def plot_halo(fig):

    nstars = 2000
    nfiles = 400 

    xs = np.loadtxt('zs.txt', unpack=True)
    xs = xs.reshape(2000, 400)

    ys = np.loadtxt('xs.txt', unpack=True)
    ys = ys.reshape(2000, 400)

    nc = 1000
    a = np.zeros((nc,nc),dtype=np.float32) 
    xl, xu = -50.0, 50.0
    yl, yu = -50.0, 50.0

    lx = xu - xl 
    ly = yu - yl 

    dx = lx/nc 
    dy = ly/nc

    print 'dx, dy=', dx, dy 

    def xloc(x):
        return int((x-xl)/dx)

    def yloc(x):
        return int((x-yl)/dy)

    def grid(xarr, yarr):

        for x, y in zip(xarr, yarr):

            if x < xl or x > xu: continue
            if y < yl or y > yu: continue 

            xpos = xloc(x)
            ypos = yloc(y)

            a[xpos, ypos] += 1.0

    for i in xrange(nstars):
        x = xs[i]
        x = x[x<1.0e6]

        y = ys[i]
        y = y[y<1.0e6]

        grid(x, y)

    a = a/a.sum()

    ax = fig.add_subplot(1, 2, 1) # Magic numbers
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.tick_params('x', which='major', pad=6)
    ax.xaxis.labelpad = 15

    ax.xaxis.tick_top()
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_label_position('top')

    s = plt.imshow(np.log10(a), cmap=cm.jet, vmin=-6, vmax=-3)
    print np.min(np.log10(a)), np.max(np.log10(a))
    plt.xlabel('$z$ [pkpc]')
    plt.ylabel('$x$ [pkpc]')

    ylabels = [-50, -25, 0, 25, 50]
    ylocs = [yloc(x) for x in ylabels]
    ylocs[-1] -= 1 
    ylabels = ['$'+str(x).strip()+'$' for x in ylabels]

    xlabels = [-50, -25, 0, 25, 50]
    xlocs = [xloc(x) for x in xlabels]
    xlocs[-1] -= 1 
    xlabels = ['$'+str(x).strip()+'$' for x in xlabels]

    plt.yticks(ylocs, ylabels)
    plt.xticks(xlocs, xlabels)

    return s

def plot_bulge(fig):

    nstars = 2000
    nfiles = 400 

    xs = np.loadtxt('zs_bulge.txt', unpack=True)
    xs = xs.reshape(2000, 400)

    ys = np.loadtxt('xs_bulge.txt', unpack=True)
    ys = ys.reshape(2000, 400)

    plt.plot(xs[350,50:], ys[350,50:])
    plt.savefig('xs.png')

    nc = 1000
    a = np.zeros((nc,nc),dtype=np.float32) 
    xl, xu = -4.0, 4.0
    yl, yu = -4.0, 4.0

    lx = xu - xl 
    ly = yu - yl 

    dx = lx/nc 
    dy = ly/nc

    print 'dx, dy=', dx, dy 

    def xloc(x):
        return int((x-xl)/dx)

    def yloc(x):
        return int((x-yl)/dy)

    def grid(xarr, yarr):

        for x, y in zip(xarr, yarr):

            if x < xl or x > xu: continue
            if y < yl or y > yu: continue 

            xpos = xloc(x)
            ypos = yloc(y)

            a[xpos, ypos] += 1.0

    for i in xrange(nstars):
        x = xs[i]
        x = x[x<1.0e6]

        y = ys[i]
        y = y[y<1.0e6]

        grid(x, y)

    a = a/a.sum()

    ax = fig.add_subplot(1, 2, 2) # Magic numbers
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.tick_params('x', which='major', pad=6)
    ax.xaxis.labelpad = 15
    
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_label_position('right')

    ax.xaxis.tick_top()
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_label_position('top')
    
    s = plt.imshow(np.log10(a), cmap=cm.jet, vmin=np.log10(6.4e-7), vmax=np.log10(6.4e-4))
    print np.min(np.log10(a)), np.max(np.log10(a))
    plt.xlabel('$z$ [pkpc]')
    plt.ylabel('$x$ [pkpc]')

    ylabels = [-4, -2, 0, 2, 4]
    ylocs = [yloc(x) for x in ylabels]
    ylabels = ['$'+str(x).strip()+'$' for x in ylabels]

    xlabels = [-4, -2, 0, 2, 4]
    xlocs = [xloc(x) for x in xlabels]
    xlabels = ['$'+str(x).strip()+'$' for x in xlabels]

    plt.yticks(ylocs, ylabels)
    plt.xticks(xlocs, xlabels)

    return s 
    
nx = 2
ny = 1 
factor_x = 4.0
factor_y = 4.0
ldim = 0.3*factor_x
bdim = 0.35*factor_y
rdim = 0.25*factor_x
tdim = 0.3*factor_y
wspace = 0.8
hspace = 0.5 

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

s1 = plot_halo(fig)
s2 = plot_bulge(fig)

wr = 0.4*plotdim_x/hdim
lcorr = 0.03*factor_x/hdim
cbaxes = fig.add_axes([l+lcorr, 0.25*factor_y/vdim, factor_x/hdim, 0.04])
cbaxes.tick_params('both', which='major', length=7, width=1)
cbaxes.tick_params('x', which='major', pad=6)
cb = plt.colorbar(s1, cax=cbaxes, orientation='horizontal')
cb.set_label(r'PDF [pkpc$^{-2}$]', labelpad=10)
cb.solids.set_edgecolor("face")
cb.set_ticks([-6, -5, -4, -3])
cb.set_ticklabels(['$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$'])

wr = 0.4*plotdim_x/hdim
lcorr = 0.03*factor_x/hdim
cbaxes = fig.add_axes([l+(factor_x+wspace)/hdim-lcorr, 0.25*factor_y/vdim, factor_x/hdim, 0.04])
cbaxes.tick_params('both', which='major', length=7, width=1)
cbaxes.tick_params('x', which='major', pad=6)
cb = plt.colorbar(s2, cax=cbaxes, orientation='horizontal')
cb.set_label(r'PDF [pkpc$^{-2}$]', labelpad=10)
cb.solids.set_edgecolor("face")
t = np.log10(np.logspace(-2, 1, 4)* 6.4e-05)
cb.set_ticks(t)
cb.set_ticklabels(['$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$'])

plt.savefig('orbits.pdf')


