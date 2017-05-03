import ResolvedStellarPops as rsp
import ResolvedStellarPops.PadovaTracks as pc
import ResolvedStellarPops.graphics.GraphicsUtils as rg
import brewer2mpl
import numpy as np
import os
import matplotlib.pylab as plt
from interpolate_tracks import axes_limits_dicts, put_an_error_pt_down, load_pagb_dict

def color_tracks(photsys='UVbright'):
    filenames = ['/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.02.dat',
                 '/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.03.dat',
                 '/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.04.dat',
                 '/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.05.dat',
                 '/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.06.dat']
    for filename in filenames:
        z = float(filename.split('z')[1].replace('.dat', ''))
        pc.color_tracks(filename, logl=3, logte=4, mass=5, z=z, comments=False,
                        photsys=photsys)


def color_by_arg(bins=None, cmap=None, ax=None, fig=None, labelfmt='$%.3f$',
                 xdata=None, ydata=None, coldata=None, xlim=None, ylim=None,
                 slice_inds=None, legend=True, filter1=None, filter2=None,
                 ymag=None):
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = plt.axes()
    if bins is None:
        bins = 10
    if slice_inds is not None:
        xdata = xdata[slice_inds]
        ydata = ydata[slice_inds]
        coldata = coldata[slice_inds]
    # need the bins to be an array to use digitize.
    if type(bins) == int:
        hist, bins = np.histogram(coldata, bins=bins)
    inds = np.digitize(coldata, bins)
    uinds = np.unique(inds)
    # digitize sticks all points that aren't in bins in the final bin
    # cut that bin, or plot will be meaningless..
    if uinds[-1] == len(bins):
        uinds = uinds[:-1]
    if cmap is None:
        if 3 <= len(uinds) <= 11:
            #bmap = brewer2mpl.get_map('Spectral', 'Diverging', len(uinds))
            bmap = brewer2mpl.get_map('Paired', 'Qualitative', len(uinds))
            cols = bmap.mpl_colors
        else:
            cols = rg.discrete_colors(len(uinds), colormap='RdYlGn')
    else:
        cols = rg.discrete_colors(len(uinds), colormap=cmap)
    sub_inds = np.array([])
    nc = len(cols[0])
    colors = np.ndarray(shape=(len(xdata), nc), dtype=float)
    labs = []
    for j, i in enumerate(uinds):
        sinds, = np.nonzero(inds == i)
        N = len(sinds)
        if N == 0:
            continue
        if labelfmt != '':
            labs.append(labelfmt % bins[i])  # bins are left bin edges.
        colors[sinds] = cols[j]
        sub_inds = np.append(sub_inds, sinds)
    inds = map(int, sub_inds[:])
    np.random.shuffle(inds)

    ax.scatter(xdata[inds], ydata[inds], marker='o', s=15, edgecolors='none',
                color=colors[inds])

    if xlim is None:
        xlim = (xdata.min(), xdata.max())
    if ylim is None:
        ylim = (ydata.max(), ydata.min())
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if ymag is None:
        ymag = filter2
    if filter1 is not None and filter2 is not None:
        ax.set_xlabel('$%s-%s$' % (filter1, filter2), fontsize=20)
        ax.set_ylabel('$%s$' % ymag, fontsize=20)
        ax.tick_params(labelsize=16)
    
    if legend is True:
        # fake out the legend...
        [ax.plot(999, 999, 'o', color=cols[i], mec=cols[i], label=labs[i])
         for i in range(len(labs))]
        ax.legend(loc=0, numpoints=1, frameon=False)

    return fig, ax


def read_isoch(filename):
    '''
    files provided by Leo:
    model HB_age logL logTe Mass stage
    why the fuck do I lose the last line?!
    '''
    data = rsp.fileIO.readfile(filename)

    starts, = np.nonzero(data['age']==0)
    ends = np.roll(starts, -1) - 1

    isos = np.array([data[starts[i]: ends[i]] for i in range(len(starts))])

    return isos


#def put_an_error_pt_down(ax, x, y, xerr=0.07, yerr=0.05):
#    ''' call ax.errorbar '''
#    ax.errorbar(x, y, xerr=xerr, yerr=xerr, ecolor='red', lw=2)


def put_zy_on_plot(ax, Z, xloc=0.05, yloc=0.95):
    ''' call ax.text with axes fraction'''
    ax.text(xloc, yloc, '$Z=%g$' % Z, fontsize=20, ha='left',
            va='center', transform=ax.transAxes)


def tods_data(isoss=None):
    data = rsp.fileIO.readfile('bi_master_phil.dat')
    filter1 = 'F330W'
    filter2 = 'F435W'
    photsys = 'acs_hrc'
    all_tracks = ['/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.02.dat',
                  '/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.03.dat',
                  '/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.04.dat',
                  '/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.05.dat',
                  '/home/rosenfield/research/Proposals/Cycle22/HPHB/models/alltracks_z0.06.dat']
    zs = [float(a.split('z')[1].split('.dat')[0]) for a in all_tracks]
    if isoss is None:
        isoss = [read_isoch('%s.%s' % (tracks,photsys)) for tracks in all_tracks]
    for i, isos in enumerate(isoss):
        z = zs[i]
        snr = np.sqrt(data['SNR_F330W']**2 + data['SNR_F435W']**2)
        inds, = np.nonzero(snr > 40.)
        fig, ax = call_color_by_arg(isos, 50, filter1, filter2, photsys=photsys,
                                    cold='logte', cba_kw={'legend': True,
                                                          'xlim': (-1.5, 2.5),
                                                          'ylim': (26.4, 20.6),
                                                          'ymag': filter2,
                                                          'bins':10})
        ax.plot(data['F330W'][inds]-data['F435W'][inds], data['F330W'][inds], '.',
                color='black', zorder=0)
        put_zy_on_plot(ax, z)
        plt.savefig('m31_%s_%s_%.4f.png' % (filter1, filter2, z))
    return isoss


def call_color_by_arg(isos, nstars, filter1, filter2,
                      distAv={'dmod':24.47, 'Av': 0.206},
                      photsys='UVbright', xlim=None, ylim=None,
                      fig=None, ax=None, ymag=None, cold='mass', cba_kw={}):
    if ymag is None:
        ymag = filter2
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    color = np.array([])
    mag = np.array([])
    mass = np.array([])
    nages = len(isos)
    for i in range(nages):
        # select a few randos.
        old, = np.nonzero(isos[i]['age'] < 1.9e8)
        rands = np.random.random(len(isos[i]))
        rinds, = np.nonzero(rands < float(nstars)/float(len(isos[i])))
        inds = list(set(rinds) & set(old))
        if len(inds) == 0:
            print 'skipping an iso...'
            continue

        Mag1 = isos[i][filter1][inds]
        Mag2 = isos[i][filter2][inds]
        Mass = isos[i][cold][inds]

        kw = dict(distAv.items() + {'filter1': filter1,
                                    'filter2': filter2}.items())

        mag1 = rsp.astronomy_utils.Mag2mag(Mag1, filter1, photsys, **kw)
        mag2 = rsp.astronomy_utils.Mag2mag(Mag2, filter2, photsys, **kw)

        color = np.append(color, mag1 - mag2)
        if ymag == filter2:
            mag = np.append(mag, mag2)
        if ymag == filter1:
            mag = np.append(mag, mag1)
        mass = np.append(mass, Mass)

        def_kw = {'xdata': color, 'ydata': mag, 'coldata': mass,
                  'xlim': xlim,
                  'ylim': ylim,
                  'bins': np.arange(0.470, 0.57, step=0.01),
                  'ax': ax,
                  'fig': fig,
                  'filter1': filter1,
                  'filter2': filter2,
                  'ymag': ymag}
    cba_kw = dict(def_kw.items() + cba_kw.items())
    fig, ax = color_by_arg(**cba_kw)

    return fig, ax

def populate_hphb_cmd(zs=[0.02, 0.03, 0.04, 0.05, 0.06]):
    # hard coding!
    track_src = '/Users/phil/research/Proposals/Cycle22/HPHB/models/'
    photsys = 'UVbright'

    # nstars per isochrone to be picked
    nstars = 10

    # plot limits
    color_lims, mag_lims, line_lims = axes_limits_dicts()
    
    # push to m31
    distAv = {'dmod':24.47,
              'Av': 0.206}
    
    # for the uncer example
    errpt = {'F225W': (0, 23.9),
             'F140LP': (-2.05, 24.9)}
    # cmds to make (will die if not 225 140 vs 336...)
    filter1s = ['F225W', 'F140LP']
    filter2 = 'F336W'

    # pagb info
    pagb1_cmd = load_pagb_dict(1.0, *['F225W', 'F140LP', 'F336W'], **distAv)
    pagb1_5_cmd = load_pagb_dict(1.5, *['F225W', 'F140LP', 'F336W'], **distAv)

    for j, z in enumerate(zs):
        track_file, = rsp.fileIO.get_files(track_src, '*%.2f.dat.%s' % (z, photsys))
        isos = read_isoch(track_file)
        for filter1 in filter1s:
            # for 3 panel fig.
            if filter1 == 'F140LP':
                labelfmt = ''
                legend = False
            else:
                labelfmt = '$M=%.3fM_\odot$'
                legend = True

            # plot models
            fig, ax = call_color_by_arg(isos, nstars, filter1, filter2,
                                        xlim=color_lims['%s-%s' % (filter1, filter2)],
                                        ylim=mag_lims[filter2],
                                        photsys=photsys, cba_kw={'legend': legend})

            # add p-agb
            for cmd in [pagb1_cmd, pagb1_5_cmd]:
                ax.plot(cmd[filter1] - cmd[filter2], cmd[filter2], lw=2,
                        color='black')
            
            # add phot uncert example
            put_an_error_pt_down(ax, *errpt[filter1])
            
            # add snr lines
            put_a_line_on_it(ax, line_lims[filter1], consty=False)
            put_a_line_on_it(ax, line_lims[filter2], consty=True)

            # add Z to the plot
            put_zy_on_plot(ax, z)
            plt.savefig('cmd_%s_%s_%.4f.png' % (filter1, filter2, z))
    return fig, ax

def put_a_line_on_it(ax, val, consty=True, color='black',
                     ls='--', lw=2, filter1=None):
    """
    if consty is True: plots a constant y value across ax.xlims().
    if consty is False: plots a constant x on a plot of y vs x-y
    """
    (xmin, xmax) = ax.get_xlim()
    (ymin, ymax) = ax.get_ylim()
    xarr = np.linspace(xmin, xmax, 20)
    # y axis is magnitude...
    yarr = np.linspace(ymin, ymax, 20)
    if consty is True:
        # just a contsant y value over the plot range of x.
        ax.hlines(val, xmin, xmax, color=color, lw=lw)
        new_xarr = xarr
    if consty is False:
        # a plot of y vs x-y and we want to mark
        # where a constant value of x is
        # e.g, f814w vs f555-f814; val is f555
        new_xarr = val - yarr
        # e.g, f555w vs f555-f814; val is f814
        if filter1 is not None:
            yarr = xarr + val
            new_xarr = xarr
        ax.plot(new_xarr, yarr, ls, color=color, lw=lw)