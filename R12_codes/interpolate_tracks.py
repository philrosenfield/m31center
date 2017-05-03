import os
import numpy as np
from scipy.interpolate import interp1d, splev, splrep, splprep

import ResolvedStellarPops as rsp
import ResolvedStellarPops.PadovaTracks as pc

import re
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter
params = {'legend.fontsize': 20}
plt.rcParams.update(params)
nullfmt = NullFormatter()

def color_teffs(ts):
    logtes = np.concatenate([t.data.LOG_TE for t in ts.tracks])
    F140s = np.concatenate([t.data.F140LP for t in ts.tracks])
    Us = np.concatenate([t.data.U for t in ts.tracks])
    F225s = np.concatenate([t.data.F225W for t in ts.tracks])
    F336s = np.concatenate([t.data.F336W for t in ts.tracks])
    ord = np.argsort(logtes)
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(logtes[ord], F140s[ord]-Us[ord], label='F140W-U')
    ax.plot(logtes[ord], F225s[ord]-Us[ord], label='F225W-U')
    ax.plot(logtes[ord], F336s[ord]-Us[ord], label='F336W-U')
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_xlabel('Log Te')
    ax.set_ylabel('Color')
    ax.legend(loc=0)
    return ax

def quick_color_em(tracks_base, prefix, photsys='UVbright'):
    '''
    I wanted to do these awesome interpolations, but really needed the
    tracks to be colored. So this goes quickly through each directory
    and adds a HB.dat file that has # in the header, a HB.dat.UVbright
    file that is the output of Leo's fromHR2mags.

    sometimes leo's code has bad line endings or skips lines, i donno. so
    when reading in as TrackSet, you'll get loads of warnings...

    if this is useful, perhaps generalize the photsys, and stick it in
    galaxies (so it works on PMS, PMS.HB)
    '''
    def add_comments_to_HB_header(tracks_base, prefix):
        tracks = os.path.join(tracks_base, prefix)
        track_names = rsp.fileIO.get_files(tracks, '*F7_*PMS.HB')

        for name in track_names:
            with open(name, 'r') as t:
                lines = t.readlines()
            try:
                imode, = [i for (i, l) in enumerate(lines) if l.strip().startswith('MODE ')]
            except ValueError:
                '\n %s \n' % name

            lines[:imode + 1] = ['# ' + l for l in lines[:imode + 1]]

            oname = '%s.dat' % name

            with open(oname, 'w') as t:
                t.writelines(lines)


    def color_tracks(tracks_base, prefix, cmd):
        tracks = os.path.join(tracks_base, prefix)
        track_names = rsp.fileIO.get_files(tracks, '*F7_*PMS.HB.dat')

        for name in track_names:
            z = float(name.split('Z')[1].split('_Y')[0])
            os.system(cmd % (name, z))

    cmd = '/Users/phil/research/Italy/fromHR2mags/fromHR2mags %s ' % photsys
    # this is set for .PMS and .PMS.HB tracks
    cmd += '%s 5 6 2 %.3f'
    add_comments_to_HB_header(tracks_base, prefix)
    color_tracks(tracks_base, prefix, cmd)


def call_color_em_alltracks(tracks_base=None, prefixs=None, photsys='UVbright'):
    if tracks_base is None:
        tracks_base = '/Users/phil/research/parsec2match/S12_set/CAF09_S12D_NS/'

    if prefixs is None:
        prefixs = ['S12D_NS_Z0.0001_Y0.249',
                   'S12D_NS_Z0.0002_Y0.249',
                   'S12D_NS_Z0.0005_Y0.249',
                   'S12D_NS_Z0.001_Y0.25',
                   'S12D_NS_Z0.002_Y0.252',
                   'S12D_NS_Z0.004_Y0.256',
                   'S12D_NS_Z0.006_Y0.259',
                   'S12D_NS_Z0.008_Y0.263',
                   'S12D_NS_Z0.014_Y0.273',
                   'S12D_NS_Z0.017_Y0.279',
                   'S12D_NS_Z0.01_Y0.267',
                   'S12D_NS_Z0.02_Y0.284',
                   'S12D_NS_Z0.03_Y0.302',
                   'S12D_NS_Z0.04_Y0.321',
                   'S12D_NS_Z0.05_Y0.339',
                   'S12D_NS_Z0.06_Y0.356']

    for p in prefixs:
        quick_color_em(tracks_base, p, photsys=photsys)


def TodsData():
    gal = rsp.galaxies.galaxy(gal_file, filetype='fitstable', filter1='F300W',
                              filter2='F435W', hla=False, angst=False,
                              photsys='std')
    photsys = 'acs_hrc'
    tracks_base = '/Users/phil/research/Proposals/Cycle21/HPHB/data/'
    prefixs= ['C09_AF02I1_Z0.07_Y0.375', 'C09_AF02I1_Z0.04_Y0.321']
    #call_color_em_alltracks(tracks_base=tracks_base, prefixs=prefixs,
    #                        photsys=photsys)
    ts_indict = {'tracks_dir': tracks_base,
                 'track_search_term': '*F7_*HB*%s' % photsys}

    for prefix in prefixs:
        ts_indict['prefix'] = prefix
        ts = pc.TrackSet(**ts_indict)
        gal.plot_cmd(gal.color, gal.mag2, scatter_off=True)
        ts.plot_all_tracks('F330W', 'F435W', annotate=False, cmd=True,
                           convert_mag_kw={'photsys':'acs_hrc','dmod':24.47},
                           ax=gal.ax)
        interpolate_tracks.put_zy_on_plot(gal.ax, prefix)
        plt.savefig('%s_%s.png' % (prefix, gal.name.replace('fits','')))


def bulge_asts2matchfake():
    '''
    Turns the AST file I made for PHAT paper 1 into a matchfake format file.
    This only needs to be done once, but if it has to be repeated...
    its suuuper slow....
    '''
    # open ast file
    ast_file='/Users/phil/research/Proposals/Cycle21/HPHB/data/B01_ASTs.dat'
    with open(ast_file, 'r') as f:
        line = f.readline()
        line = f.readline()

    # grab columns
    colheads = line.replace('#', '').strip().split(', ')
    get_cols = ['mag1in', 'mag2in', 'mag1out', 'mag2out']
    use_cols = [colheads.index(i) for i in get_cols]
    asts = np.genfromtxt(ast_file, usecols=use_cols, names=get_cols)

    # matchfake file goes mag1in mag2in mag1diff mag2diff

    mag1diff = asts['mag1out'] - asts['mag1in']
    mag2diff = asts['mag2out'] - asts['mag2in']
    # why keep shit?
    mask, = np.nonzero((mag1diff < 20.) & (mag2diff < 20.))

    matchdata = np.array([asts['mag1in'], asts['mag2in'], mag1diff, mag2diff]).T[mask]

    # write out matchfake
    matchfake = '/Users/phil/research/Proposals/Cycle21/HPHB/data/B01_F275W_F336W.matchfake'
    np.savetxt(matchfake, matchdata, fmt='%.4f')
    return asts


def compare_hess(mbinsize, cbinsize):
    prefixs = ['S12D_NS_Z0.04_Y0.321',
               'S12D_NS_Z0.05_Y0.339',
               'S12D_NS_Z0.06_Y0.356']

    ts_indict = {'tracks_dir': '/Users/phil/research/parsec2match/S12_set/CAF09_S12D_NS/',
                 'track_search_term': '*F7_*HB*UVbright'}

    itss = []
    for i in range(len(prefixs)):
        ts_indict['prefix'] = prefixs[i]
        itss.append(InterpolateTracks(ts_indict, max_mass=0.6, min_mass=0.465))

    [make_uvbright_sim(its=itss[i]) for i in range(len(itss))]

    color_lims, mag_lims, line_lims = axes_limits_dicts()
    filters = ['F275W', 'F336W', 'F225W', 'F140LP']
    filters = ['F225W', 'F140LP']
    filter2 = 'F336W'

    for filter1 in filters:
        if filter1 == filter2:
            continue
        xlim = color_lims['%s-%s' % (filter1, filter2)]
        ylim = mag_lims[filter2]
        # want centers of bins, and want the same bins for each filter set
        # so we can see differences.
        cbins = np.arange(*np.asarray(xlim) + cbinsize, step=cbinsize)
        mbins = np.arange(*np.asarray(ylim[::-1]) + mbinsize, step=mbinsize)
        sim = simulation()
        sim.filter1 = filter1
        sim.filter2 = filter2

        for its in itss:
            its.ast_corrections(filter1, filter2)
            recs = its.recovered_stars(filter1, filter2)
            its.filter1 = filter1
            its.filter2 = filter2
            its.mag1 = its.__dict__['m%ss_cor' % filter1]
            its.mag2 = its.__dict__['m%ss_cor' % filter2]
            its.color = its.mag1 - its.mag2
            incmd = within_cmd(its.color, its.mag2, xlim, ylim)
            slice_inds = list(set(recs) & set(incmd))
            sim.mag1 = np.append(sim.mag1, its.mag1[slice_inds])
            sim.mag2 = np.append(sim.mag2, its.mag2[slice_inds])

            its.make_hess(mbinsize, slice_inds=slice_inds,
                          hess_kw={'cbin': cbins, 'mbin': mbins})

            # plot them ... why not ...
            # ax = its.hess_plot(imshow_kw={'interpolation': 'gaussian'})
            #put_zy_on_plot(ax, its.prefix)
            #plt.savefig('%s_%s_%s_hess.png' % (its.prefix, filter1, filter2), dpi=300)
            #plt.close()
            fig = plt.figure()
            ax = plt.axes()
            its.color_by_arg(0,0,0, xdata=its.color, ydata=its.mag2,
                             coldata=its.mass_interps, xlim=xlim, ylim=ylim,
                             ax=ax, fig=fig, slice_inds=slice_inds,
                             bins=its.masses)
            plt.savefig('%s_%s_%s_cmd_bymass.png' % (its.prefix, filter1, filter2), dpi=300)
            plt.close()


        ihesses = [its.hess[2] for its in itss]
        nhesses = len(ihesses)
        # fill big sim..
        sim.color = sim.mag1 - sim.mag2
        sim.hess = rsp.astronomy_utils.hess(sim.color, sim.mag2, mbinsize,
                                            **{'cbin':cbins, 'mbin':mbins})
        ax = rsp.astronomy_utils.hess_plot(sim.hess)

        mean_hess = (ihesses[0] + ihesses[1] + ihesses[2])/3.
        frac_hesses = [(ihesses[i]-mean_hess)/mean_hess for i in range(nhesses)]
        imshow_kw = {'norm': None,
                     'cmap': cm.RdBu}
        for i, frac_hess in enumerate(frac_hesses):
            hess = (cbins, mbins, frac_hess)
            imshow_kw['vmax'] = np.nanmax(np.abs(frac_hess))
            imshow_kw['vmin'] = -1. * imshow_kw['vmax']
            ax = rsp.astronomy_utils.hess_plot(hess, filter1=filter1,
                                               filter2=filter2, imshow_kw=imshow_kw,
                                               colorbar=True)
            put_zy_on_plot(ax, itss[i].prefix)
            plt.savefig('%s_%s_%s_hess_fracdiff.png' % (itss[i].prefix, filter1, filter2), dpi=300)
            plt.close()

    return itss


def HRDfigure(ts):
    '''
    This is for the HRD figure in the intro but saved at each track draw.
    '''
    def read_pagb(pagb_file):
        with open(pagb_file) as f:
            colhead = f.readline()
            col_keys = colhead.replace('#', '').split()
            pagb_data = np.genfromtxt(f, names=col_keys)
        return pagb_data

    smooth = rsp.math_utils.smooth
    plt_kw = {'lw': 3}

    peagb = {'file': '/Users/phil/research/PHAT/tracks/hbz07y389/Z0.07Y0.389OUTA1.72_F7_M.580.PMS.HB.phat'}
    agbm = {'file': '/Users/phil/research/Proposals/Cycle21/HPHB/data/C09_AF02I1_Z0.04_Y0.321/Z0.04Y0.321OUTA1.74_F7_M.470.PMS.HB.res0.00.UVbright.dat'}
    agb = {'file': '/Users/phil/research/PHAT/tracks/hbz07y389/Z0.07Y0.389OUTA1.72_F7_M1.60.PMS.HB.phat'}
    pagb = {'file': '/Users/phil/research/PHAT/tracks/vw_extended/m0.597z016H.dat.phat'}
    zahb = {'file': None, 'arrow': None}

    peagb['pcolor'] = 'red'
    agbm['pcolor'] = 'blue'
    agb['pcolor'] = 'black'
    pagb['pcolor'] = 'purple'
    zahb['pcolor'] = 'darkgreen'

    peagb['arrow'] = [150, 600, 800]
    agbm['arrow'] = [150, 500, 700]
    agb['arrow'] = [700, 1200]
    pagb['arrow'] = [10, 20]

    peagb['xy'] = (4.5, 3.17)
    agbm['xy'] = (4.9, 2.2)
    agb['xy'] = (3.7, 3.4)
    pagb['xy'] = (5, 3.8)
    zahb['xy'] = (4, 1.48)

    peagb['name'] = 'PE-AGB'
    agbm['name'] = "AGB-Manqu\\'e"
    agb['name'] = 'AGB'
    pagb['name'] = 'P-AGB'
    zahb['name'] = 'ZAHB'

    # ZAHB
    zahb['logte'] = smooth(np.array([t.data.LOG_TE[1:4] for t in ts.tracks]).ravel())
    zahb['logl'] = smooth(np.array([t.data.LOG_L[1:4] for t in ts.tracks]).ravel())

    # PAGB
    pagb_data = read_pagb(pagb_file)
    pagb['logte'] = pagb_data['logte']
    pagb['logl'] = pagb_data['logl']

    # the rest
    track_dict = {'peagb': peagb, 'agbm': agbm, 'agb': agb}
    pctrack_kw = {'min_lage': 0, 'cut_long': False}
    for key, tdict in track_dict.items():
        track = pc.Track(tdict['file'], **pctrack_kw)
        inds = None
        if key == 'agbm':
            inds = np.arange(len(track.data.LOG_TE))[:-170]
        elif key == 'peagb':
            inds = np.arange(600, len(track.data.LOG_TE))
        elif key == 'agb':
            inds = np.arange(rsp.math_utils.closest_match(3.74, track.data.LOG_L)[0])
        tdict['logl'] = track.data.LOG_L[inds]
        tdict['logte'] = track.data.LOG_TE[inds]
        tdict['track'] = track
    track_dict['pagb'] = pagb
    track_dict['zahb'] = zahb

    # plots.
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes()
    ax.set_ylim(1., 4.1)
    ax.set_xlim(5.2, 3.3)
    for key, tdict in track_dict.items():
        logte = tdict['logte']
        logl = tdict['logl']
        ax.plot(logte, logl, color=tdict['pcolor'], **plt_kw)
        if tdict['arrow'] is not None:
            [rsp.graphics.GraphicsUtils.arrow_on_line(ax, logte, logl, arr_pt,
                                                      tdict['pcolor'])
                                                  for arr_pt in tdict['arrow']]
        ax.annotate("$%s$" % tdict['name'], xy=tdict['xy'],
                    color=tdict['pcolor'], **ann_kwargs)

    xlab ='$\log \ T_{{\\rm eff}} \ (K)$'
    ylab = '$\log \ L \ (L_\odot$)'
    ax.set_xlabel(xlab, fontsize=20)
    ax.set_ylabel(ylab, fontsize=20)
    plt.text(0.5, 1.06, '$F140LP-F336W$',
         horizontalalignment='center',
         fontsize=20,
         transform = ax.transAxes)

    ax.tick_params(labelsize=16)
    ax.set_ylim(1., 4.1)
    put_filters_top(ax, ts)
    exts = ('.png', '.eps')
    [plt.savefig('hrd%s' % ext, dpi=300, bbox_to_inches='tight') for ext in exts]

    return


def put_filters_top(ax, ts):
    F140s = np.concatenate([t.data.F140LP for t in ts.tracks])
    logls = np.concatenate([t.data.LOG_L for t in ts.tracks])
    logtes = np.concatenate([t.data.LOG_TE for t in ts.tracks])
    F336s = np.concatenate([t.data.F336W for t in ts.tracks])
    ord = np.argsort(logtes)
    logtes = logtes[ord]
    F336s = F336s[ord]
    F140s = F140s[ord]
    cut, = np.nonzero(logtes <= 4.8)
    nocut, = np.nonzero(logtes >= 4.8)
    color = F140s-F336s
    nocut2, = np.nonzero(color[nocut] <=-2)
    cut = np.concatenate([cut, nocut[nocut2]])
    logtes = logtes[cut]
    F336s = F336s[cut]
    F140s = F140s[cut]
    non_dupes = ts.tracks[0].remove_dupes(logtes, F336s, F140s)

    logtes = logtes[non_dupes]
    F336s = F336s[non_dupes]
    F140s = F140s[non_dupes]
    color = F140s - F336s
    interp_arr = np.linspace(3.3, 5.2, 10000)
    interp = interp1d(logtes, color, bounds_error=0)
    model_arr = interp(interp_arr)

    colvals = [8., 4., 0., -1, -2, -2.386]

    minds = [np.nanargmin(np.abs(c-model_arr)) for c in colvals]
    xdata  = interp_arr[minds]
    coords = 'figure fraction'
    text=(0.5, 0.95)
    offset_text=0.05,
    offset_tick=-0.08,
    shk=0.4
    for i in range(len(colvals)):
        ax.annotate(r'$%.2f$' % colvals[i],
                    xytext=(xdata[i], ax.get_ylim()[1] + offset_text),
                    xy=(xdata[i], ax.get_ylim()[1] + offset_tick),
                    arrowprops=dict(facecolor='grey', width=1, shrink=shk,
                                    lw=0.0, headwidth=0), fontsize=16,
                                    ha='center')
    for tick in np.concatenate((ax.xaxis.get_major_ticks(),ax.xaxis.get_minor_ticks())):
        tick.tick2On = False
    return ax


def cmd_sim_figs(filter1='F225W', filter2='F336W', prefixs=None,
                 ts_indict=None):
    '''
    this makes fig2s

    # use S12??
    prefixs = ['S12D_NS_Z0.04_Y0.321',#],
               #'S12D_NS_Z0.05_Y0.339' ,
               'S12D_NS_Z0.06_Y0.356'][::-1]

    ts_indict = {'tracks_dir': '/Users/phil/research/parsec2match/S12_set/CAF09_S12D_NS/',
                 'track_search_term': '*F7_*HB*UVbright'}
    '''
    exts = ['png', 'pdf']
    filter1s = ['F225W', 'F140LP']
    filter2 = 'F336W'

    color_lims, mag_lims, line_lims = axes_limits_dicts()
    errpt = {filter1s[0]: (0, 23.9), filter1s[1]: (-2.05, 24.9)}

    distAv = {'dmod':24.47, 'Av': 0.206}

    if prefixs is None:
        prefixs=['C09_AF02I1_Z0.07_Y0.375', 'C09_AF02I1_Z0.04_Y0.321']

    if ts_indict is None:
        ts_indict = {'tracks_dir': '/Users/phil/research/Proposals/Cycle21/HPHB/data/',
                     'track_search_term': '*F7_*HB*UVbright'}

    itss = []
    mass_array = np.arange(.465, .5, .005)
    for i in range(len(prefixs)):
        ts_indict['prefix'] = prefixs[i]
        itss.append(InterpolateTracks(ts_indict, mass_array=mass_array))

    [make_uvbright_sim(its=itss[i]) for i in range(len(itss))]

    for filter1 in filter1s:
        xlim = color_lims['%s-%s' % (filter1, filter2)]
        ylim = mag_lims[filter2]
        pagb1_cmd = load_pagb_dict(1.0, *[filter1, filter2], **distAv)
        pagb1_5_cmd = load_pagb_dict(1.5, *[filter1, filter2], **distAv)

        for its in itss:
            figfmt = '_'.join(map(str, (its.prefix, filter1, filter2)))
            its.ast_corrections(filter1, filter2)
            recs = its.recovered_stars(filter1, filter2)
            its.filter1 = filter1
            its.filter2 = filter2
            its.mag1 = its.__dict__['m%ss_cor' % filter1]
            its.mag2 = its.__dict__['m%ss_cor' % filter2]
            its.color = its.mag1 - its.mag2
            incmd = within_cmd(its.color, its.mag2, xlim, ylim)
            slice_inds = list(set(recs) & set(incmd))

            fig = plt.figure(figsize=(8, 8))
            ax = plt.axes()

            if filter1 == filter1s[1]:
                labelfmt = ''
                legend = False
            else:
                labelfmt = '$M=%.3fM_\odot$'
                legend = True

            cba_kw = {'xdata': its.color,
                      'ydata': its.mag2,
                      'coldata': its.mass_interps,
                      'xlim': xlim,
                      'ylim': ylim,
                      'bins': np.unique(its.mass_interps),
                      'ax': ax,
                      'fig': fig,
                      'labelfmt': labelfmt,
                      'legend': legend,
                      'slice_inds': slice_inds}

            ax = its.color_by_arg(0, 0, 0, **cba_kw)

            for cmd in [pagb1_cmd, pagb1_5_cmd]:
                ax.plot(cmd[filter1] - cmd[filter2], cmd[filter2], lw=2,
                        color='black')

            put_an_error_pt_down(ax, *errpt[filter1])
            put_zy_on_plot(ax, its.prefix)
            its.put_a_line_on_it(ax, line_lims[filter1], consty=False,
                                 annotate=False)
            its.put_a_line_on_it(ax, line_lims[filter2], consty=True,
                                 annotate=False)
            ax.set_xlabel('$%s-%s$' % (filter1, filter2), fontsize=20)
            ax.set_ylabel('$%s$' % filter2, fontsize=20)

            for ext in exts:
                figname = '%s.%s' % (figfmt, ext)
                plt.savefig(figname, dpi=300, bbox_inches='tight')
                print 'wrote %s' % figname
    return


def put_an_error_pt_down(ax, x, y, xerr=0.07, yerr=0.05):
    ax.errorbar(x, y, xerr=xerr, yerr=xerr, ecolor='red', lw=2)


def put_zy_on_plot(ax, prefix):
    Z, Y = map(float, prefix.replace('Z','').replace('Y','').split('_')[2:])
    ax.text(0.05, 0.95, '$Z=%.3f\ Y=%.3f$' % (Z, Y), fontsize=20, ha='left',
            va='center', transform=ax.transAxes)


def within_cmd(color, mag, clim, mlim):
    # these are magnitudes...
    if mlim[0] < mlim[1]:
        mlim = mlim[::-1]

    inds, = np.nonzero((color < clim[1]) & (color > clim[0]) &
                       (mag < mlim[0]) & (mag > mlim[1]))
    return inds


def axes_limits_dicts():
    '''
    place to house axes limits.
    '''
    color_lims = {'F275W-F336W': (-1, 2),
                  'F225W-F336W': (-1, 1.5),
                  'F140LP-F336W':  (-2.5, 0),
                  'F330W-F435W': (-1.75, 1.3)}

    mag_lims = {'F336W': (26.5, 19.5),
                'F435W': (26., 18.3)}

    line_lims = {'F336W': 26.0, 'F140LP': 24.8, 'F225W': 25.1,
                 'F435W': 99, 'F330W': 99}
    return color_lims, mag_lims, line_lims


def new_tracks():
    prefix = 'D13_Z0.03_Y0.34'
    tracks_dir = '/Users/phil/research/Proposals/Cycle21/HPHB/data/'
    #quick_color_em(tracks_dir, prefix)

    default = {'prefix': prefix,
               'tracks_dir': tracks_dir,
               'track_search_term': '*F7_*HB*UVbright'}
    make_uvbright_sim(in_dict=default)


def make_uvbright_sim(total_mass=2e4, sim_cmd_kw={},
                      ts_dict={}, its=None, max_hbmass=0.6):

    sim_cmd_default = {'dmod': 24.47, 'Av': 0.206, 'delta_t': 1e4}
    sim_cmd_kw = dict(sim_cmd_kw.items() + sim_cmd_default.items())


    if its is None:
        ts_default = {'tracks_dir': '/Users/phil/research/parsec2match/S12_set/CAF09_S12D_NS/',
                      'track_search_term': '*F7_*HB*UVbright'}
        ts_indict = dict(ts_default.items() + in_dict.items())
        its = InterpolateTracks(ts_indict, max_mass=max_hbmass)

    act_mass = its.simulate_cmd(**sim_cmd_kw)

    if act_mass < total_mass:
        print act_mass
        sys.exit()
    return act_mass

def plot_uvbright_sim(filter1s, filter2, xlims=None, ylim=None, plt_kw={},
                      sim_cmd_kw={}):
    print 'i broke this.'
    #filter1s = ['F275W', 'F225W', 'F140LP']
    filter1s = ['F275W']
    filter2 = 'F336W'

    ylim = None
    xlims = None
    plt_kw = {'color': 'grey', 'linestyle': '-', 'label': None}

    for filter1, xlim in zip(filter1s, xlims):
        its.ast_corrections(filter1, filter2)
        its.filter1 = filter1
        its.filter2 = filter2
        its.mag2 = its.__dict__['m%ss_cor' % filter2]
        its.color = its.__dict__['m%ss_cor' % filter1] - its.mag2

        #inds, = np.nonzero(np.isfinite(its.mag2) & np.isfinite(its.color) &
        #                   (np.abs(its.mag2) < 40) & (np.abs(its.color) < 40))

        inds, = np.nonzero((its.color < xlim[1]) & (its.color > xlim[0]) &
                           (its.mag2 < ylim[0]) & (its.mag2 > ylim[1]))

        fig = plt.figure()
        ax = plt.axes()
        its.plot_cmd(its.color, its.mag2, xlim=xlim, ylim=ylim, ax=ax, fig=fig,
                     scatter_off=True, slice_inds=inds)
        #its.make_hess(0.03, hess_kw={'cbinsize':0.03}, slice_inds=inds)
        #imshow_kw = {'interpolation': 'gaussian'}
        #its.hess_plot(fig=fig, ax=ax, imshow_kw=imshow_kw)


        for t in its.tracks:
            filt_Mag_dict = {filter1: t.data['%s' % filter1],
                             filter2: t.data['%s' % filter2]}
            its.correct_extinction_dist(t, t.data.LOG_TE, filt_Mag_dict, **sim_cmd_kw)

            mag2 = t.__dict__['m%s' % filter2]
            color = t.__dict__['m%s' % filter1] - mag2
            its._plot_cmd(color, mag2, filter1, filter2, xlim=xlim, ylim=ylim,
                          ax=ax, fig=fig, plt_kw=plt_kw)
        '''
        figname = '%s.png' % '_'.join((its.prefix, filter1, filter2))
        plt.savefig(figname, dpi=300)
        print 'wrote %s' % figname
        '''

    return act_mass

#BROKEN WITH MORGANS REFACTOR
#class InterpolateTracks(rsp.galaxies.starpop, pc.TrackSet):
#    def __init__(self, ts_indict, max_mass=120., min_mass=0., mass_array=None):
#        pc.TrackSet.__init__(self, **ts_indict)
#        rsp.galaxies.star_pop.__init__(self)
#
#        filename = '/Users/phil/research/Proposals/Cycle21/HPHB/data/ext_UVbright.dat'
#        with open(filename,'r') as f:
#            line = f.readline()
#            col_keys = line.replace('#','').strip().split()
#        self.ext_tab = np.genfromtxt(filename, names=col_keys)
#        self.cut_masses(max_mass=max_mass, min_mass=min_mass, mass_array=mass_array)
#
#
#    def cut_masses(self, max_mass=120., min_mass=0., mass_array=None):
#        if mass_array is None:
#            self.tracks = [t for t in self.tracks if t.mass > min_mass and
#                                                     t.mass < max_mass]
#        else:
#            self.tracks = [t for t in self.tracks if '%.3f' % t.mass in ['%.3f' % m for m in mass_array]]
#
#    def cut_wd(self):
#        '''
#        adds attribute to track with inds of non-wd.
#        '''
#        for t in self.tracks:
#            logg = -10.616 + np.log10(t.mass) + 4.0 * t.data.LOG_TE - t.data.LOG_L
#            t.not_wd, = np.nonzero(logg < 6.5)
#
#    def get_extinctions(self, logte, logg=4.5, Rv=3.1, filters=None):
#        A2 = 'no data'
#        approx_te = np.round(10**logte, decimals=-3)
#        table_te = np.unique(self.ext_tab['Teff'])
#        (matched_inds, diffs) = zip(*[rsp.math_utils.closest_match(rt, table_te) for rt in approx_te])
#        inds = np.array([np.nonzero((self.ext_tab['Teff'] == table_te[i]) &
#                                    (self.ext_tab['logg'] == logg) &
#                                    (self.ext_tab['Rv'] == Rv))[0]
#                                    for i in matched_inds]).ravel()
#
#        for filt in filters:
#            try:
#                self.__setattr__('A%s' % filt, self.ext_tab['A%s' % filt][inds])
#            except IndexError:
#                print '%s not found.' % filt
#
#    def recovered_stars(self, filter1, filter2, mlim=90):
#        m1 = self.__dict__['m%ss_cor' % filter1]
#        m2 = self.__dict__['m%ss_cor' % filter2]
#        recs, = np.nonzero((m1 < mlim) & (m2 < mlim))
#        return recs
#
#    def correct_extinction_dist(self, track, logte, filt_Mag_dict, dmod=None, Av=None):
#        self.get_extinctions(logte, filters=filt_Mag_dict.keys())
#        for filt, Mag in filt_Mag_dict.items():
#            Al = self.__dict__['A%s' % filt]
#            track.__setattr__('m%s' % filt, Mag + dmod + Al * Av)
#
#    def squish(self, *attrs, **kwargs):
#        '''
#        concatenates an attribute or many attributes and adds them to galaxies
#        instance -- with an 's' at the end to pluralize them... that might
#        be stupid.
#        ex
#        t.squish('mF275W', 'mag2', 'ra', 'dec')
#        gs.ras ...
#        '''
#        slice_attr = kwargs.get('slice_attr')
#        if slice_attr is not None:
#            for attr in attrs:
#                try:
#                    new_list = [t.__getattribute__(attr)[t.__dict__[slice_attr]] for t in self.tracks]
#                except AttributeError:
#                    new_list = [t.data[attr][t.__dict__[slice_attr]] for t in self.tracks]
#                new_val = np.concatenate(new_list)
#                self.__setattr__('%ss' % attr, new_val)
#        else:
#            for attr in attrs:
#                try:
#                    new_list = [t.__getattribute__(attr) for t in self.tracks]
#                except AttributeError:
#                    new_list = [t.data[attr] for t in self.tracks]
#                new_val = np.concatenate(new_list)
#                self.__setattr__('%ss' % attr, new_val)
#
#    def fixed_dtime_interpolation(self, delta_t=5e4):
#        '''
#        should this automatically do all the tracks? or just one...
#        '''
#        for track in self.tracks:
#            age = track.data.AGE
#            mode = track.data.MODE
#
#            max_age = np.max(age)
#            relative_age = age/max_age
#            nages = max_age/delta_t
#
#            track.interp_arr = np.linspace(0., 1, nages)
#            track.interp = interp1d(relative_age, mode, bounds_error=0)
#            track.model_arr = track.interp(track.interp_arr)
#        return
#
#    def interpolate_from_mode(self, *cols):
#        for track in self.tracks:
#            mode = track.data.MODE
#            for col in cols:
#                track.__setattr__('%s_interp' % col, interp1d(mode,
#                                                              track.data[col],
#                                                              bounds_error=0))
#        return
#
#    def simulate_cmd(self, filters=['F275W', 'F336W', 'F225W', 'F140LP'],
#                     Av=None, dmod=None, delta_t=5e4):
#        '''
#        interpolates tracks in fixed time steps, delta_t
#        creates interpolated cmd as self.m<filter> or self.M<filter> if no
#        dmod and Av supplied.
#        '''
#        tot_mass = 0.
#
#        self.fixed_dtime_interpolation(delta_t=delta_t)
#        self.interpolate_from_mode('LOG_TE', *filters)
#
#        for track in self.tracks:
#            fofte = track.LOG_TE_interp(track.model_arr)
#            filt_mag_dict = {}
#            for filt in filters:
#                filt_mag_dict[filt] = track.__dict__['%s_interp' % filt](track.model_arr)
#
#            if Av is not None or dmod is not None:
#                self.correct_extinction_dist(track, fofte, filt_mag_dict, Av=Av,
#                                             dmod=dmod)
#            else:
#                [track.__setattr__('M%s' % filt, Mag) for filt, Mag in
#                                 filt_Mag_dict.items()]
#
#            tot_mass += float(len(filt_mag_dict[filt])) * track.mass
#            track.__setattr__('mass_interp', np.repeat(track.mass, len(filt_mag_dict[filt])))
#
#        if Av is not None or dmod is not None:
#            self.squish(*['m%s' % f for f in filters])
#        else:
#            self.squish(*['M%s' % f for f in filters])
#
#        self.squish('mass_interp')
#        return tot_mass
#
#    def ast_corrections(self, filter1, filter2, mag1_shift=1., mag2_shift=1.):
#
#        ast_file = '/Users/phil/research/Proposals/Cycle21/HPHB/data/B01_%s_%s.matchfake'
#        mag1 = self.__dict__['m%ss' % filter1]
#        mag2 = self.__dict__['m%ss' % filter2]
#
#        b01_asts = rsp.galaxies.artificial_star_tests(ast_file % (filter1,
#                                                                  filter2))
#        b01_asts.mag1 += mag1_shift
#        b01_asts.mag2 += mag2_shift
#
#        m1_cor, m2_cor = b01_asts.ast_correction(mag1, mag2)
#        self.__setattr__('m%ss_cor' % filter1, m1_cor)
#        self.__setattr__('m%ss_cor' % filter2, m2_cor)
#
#    def _plot_cmd(self, color, mag, filter1, filter2, xlim=None, ylim=None,
#                  fig=None, ax=None, plt_kw={}, save_plot=False):
#        if fig is None:
#            fig = plt.figure()
#
#        if ax is None:
#            ax = plt.axes()
#
#        if xlim is None:
#            xlim = (-1, 2)
#
#        if ylim is None:
#            ylim = (25.4, 22)
#
#        label = '\ '.join(self.prefix.replace('Z','Z=').replace('Y','Y=').split('_')[2:])
#        default = {'color': 'black', 'marker': '.', 'linestyle': 'none', 'label': label}
#        if plt_kw.get('linestyle') is not None:
#            del default['marker']
#        plt_kw = dict(default.items() + plt_kw.items())
#
#        ax.plot(color, mag, **plt_kw)
#
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#        ax.set_xlabel('$%s-%s$' % (filter1, filter2), fontsize=20)
#        ax.set_ylabel('$%s$' % filter2, fontsize=20)
#        ax.tick_params(labelsize=16)
#        ax.legend(numpoints=1, loc=0, frameon=0)
#
#        if save_plot is True:
#            figname = '%s.png' % '_'.join((self.prefix, filter1, filter2))
#            plt.savefig(figname, dpi=300)
#            print 'wrote %s' % figname
#
#    def asts2new_mags(self):
#        '''
#        This interpolates to find an estimate of the ast corrections for
#        filters we propose to use.
#
#        First it shifts all mags to m31 distance and av.
#        Then it interpolates a mean color vs mean color line for each filter
#        combo (275-336 vs x-336)
#        Then it cycles through the matchfake file and applies the interpolation
#        to f140lp and 225. It's pretty slow, and makes a silly assumption
#        in the end that's kinda hidden.
#        f140lp-f336 - (f275-f336) leaves f140lp-f275. To get f275, I just
#        take the astmag1 from the matchfake file, not part of the interpolation
#        It's not cool, but this is a rough estimate, who is to say the
#        photometric conditions will be similar in these filters anyway...
#        so maybe a shitty premise.
#        '''
#        # ready the filters
#        filters= ['F275W', 'F140LP', 'F225W', 'F336W']
#        self.squish(*filters)
#        self.squish('LOG_TE')
#
#        filt_mag_dict = dict(((f, self.__dict__['%ss' % f]) for f in filters))
#
#        self.correct_extinction_dist(self, self.LOG_TEs, filt_mag_dict,
#                                     Av=0.206, dmod=24.47)
#
#        colors = np.array([self.__dict__['m%s' % f] -
#                           self.__dict__['m%s' % filters[-1]]
#                           for f in filters[:-1]])
#
#        # cut in color space (could be done by astmag range)
#        inds, = np.nonzero((colors[0] <= 2.5) & (colors[0] >= -1.))
#        tmp = colors.T[inds]
#        colors = tmp.T
#
#        # bin up the 275-336 color and cycle through adding means
#        bins0 = np.arange(np.min(colors[0]), np.max(colors[0]), 0.1)
#        inds0 = np.digitize(colors[0], bins0)
#        mcolor0 = np.array([])
#        mcolor1 = np.array([])
#        mcolor2 = np.array([])
#        for i in np.unique(inds0):
#            if i == len(bins0):
#                continue
#            inds = np.nonzero(inds0 == i)
#            mcolor2 = np.append(mcolor2, np.mean(colors[2][inds]))
#            mcolor1 = np.append(mcolor1, np.mean(colors[1][inds]))
#            mcolor0 = np.append(mcolor0, bins0[i])
#
#        # load the matchfake file.
#        matchfake = '/Users/phil/research/Proposals/Cycle21/HPHB/data/B01_F275W_F336W.matchfake'
#        mag1in, mag2in, mag1diff, mag2diff = np.loadtxt(matchfake, unpack=True)
#        mag1out = mag1diff + mag1in
#        mag2out = mag2diff + mag2in
#
#        colorin = mag1in - mag2in
#        colorout = mag1out - mag2out
#
#        match1fake = '/Users/phil/research/Proposals/Cycle21/HPHB/data/B01_F140LP_F336W.matchfake'
#        match2fake = '/Users/phil/research/Proposals/Cycle21/HPHB/data/B01_F225W_F336W.matchfake'
#
#        for mcolor, mf in zip((mcolor1, mcolor2), (match1fake, match2fake)):
#            tckp, u = splprep([mcolor0, mcolor], s=0)
#            col0new, colnew = splev(np.linspace(0, 1, 1000), tckp)
#
#            magin = np.array([])
#            magout = np.array([])
#            corresponding_mag2in = np.array([])
#            corresponding_mag2out = np.array([])
#
#            for i in range(len(colorin)):
#                if colorin[i] < col0new.min() or colorin[i] > col0new.max():
#                    continue
#                in_ind, diff = rsp.math_utils.closest_match(colorin[i], col0new)
#                if diff > 0.01:
#                    print 'bad', diff, colorin[i], col0new[in_ind]
#                magin = np.append(magin, colnew[in_ind] - col0new[in_ind] + mag1in[i])
#                corresponding_mag2in = np.append(corresponding_mag2in, mag2in[i])
#
#                out_ind, diff = rsp.math_utils.closest_match(colorout[i], col0new)
#                if diff > 0.01:
#                    #print 'ok', diff, colorout[i], col0new[out_ind]
#                    magout = np.append(magout, 99.)
#                else:
#                    magout = np.append(magout, colnew[out_ind] - col0new[out_ind] + mag1out[i])
#                corresponding_mag2out = np.append(corresponding_mag2out, mag2out[i])
#
#            magdiff = magout - magin
#            cmag2diff = corresponding_mag2out - corresponding_mag2in
#            mdata = np.array([magin, corresponding_mag2in, magdiff, cmag2diff]).T
#
#            np.savetxt(mf, mdata, fmt='%.4f')
#            print 'wrote %s' % mf
#        return
#
#    def check_leos_model(self, dmod=24.47, Av=0.206, af275w=2.050, af336w=1.671):
#        filename = '/Users/phil/research/PHAT/UVBrightpaper/supporting_files/lf_C_Z0.07_Y0.389.dat'
#        Mag, n275, n336 = np.loadtxt(filename, unpack=True)
#        mF336W = Mag + dmod + Av * af336w
#        mF275W = Mag + dmod + Av * af275w
#
#        bins = np.arange(16,32,0.1)
#        hist, bins = np.histogram(self.mF275Ws, bins=bins)
#        norm = n275.max()/hist.max()
#        plt.plot(bins[1:], hist*norm)
#
#        ax = plt.axes()
#        ax.plot(bins[1:], hist*norm, label='new\ sim')
#        ax.plot(mF275W, n275, label='leo\ sim')
#        ax.set_xlabel('F275W',fontsize=20)
#        ax.legend(loc=0, numpoints=1)
#        plt.savefig('comp_lf_275.png')
#
#        # found max by hand...
#        # hist[np.nonzero(np.round(bins, 1) == 26.4)]
#        # np.nonzero(np.round(bins, 1) == 26.4)
#        hist, bins = np.histogram(self.mF336Ws, bins=bins)
#        norm = n336.max()/hist[104]
#
#        plt.figure()
#        ax = plt.axes()
#        ax.plot(mF336W, n336, label='leo\ sim',color='green')
#        ax.plot(bins[1:], hist*norm, label='new\ sim', color='blue')
#        ax.set_xlabel('F336W',fontsize=20)
#        ax.legend(loc=0, numpoints=1)
#        plt.savefig('comp_lf_336.png')
#

# From cycle 20
def load_vw(H=True, mass=None):
    '''
    If H == True: gives hydrogen pagb track, else he pagb track
    if mass specified, returns the table and the indices of said mass.
    if mass not specified, returns the table and list of indices of each mass.
    '''
    names = ['mass', 'Y', 'Z', 'time', 'Teff', 'logl','mbol', 'F140LP',
             'F175Wfoc', 'F275Wfoc', 'F200LP', 'F218W', 'F225W', 'F275W',
             'F300X', 'F336W', 'F25QTZ', 'FUV', 'NUV', 'U', 'B', 'V']
    if H==True:
        filename = '/Users/phil/research/Proposals/Cycle21/HPHB/data/H-pnn.tab.UVbright'
    else:
        filename = '/Users/phil/research/Proposals/Cycle21/HPHB/data/He-pnn.tab.UVbright'

    vw_tab = np.genfromtxt(filename, names=names)

    if mass is None:
        un = np.unique(vw_tab['mass'])
        # this is the first one
        iMs = [list(vw_tab['mass']).index(u) for u in un]
        # the indices of each
        Ms = [np.arange(iMs[i], iMs[i+1]) for i in range(len(iMs)-1)]
        # the last one
        Ms.append(np.arange(iMs[i+1], len(vw_tab['mass'])))
    else:
        Ms, = np.nonzero(vw_tab['mass']==mass)
    return vw_tab, Ms


def load_pagb_dict(pagb_mass, *filters, **distAv_kw):
    '''
    Loads pagb tracks of a given mass to plot with HPHBs
    returns dict of the cmd in *filters
    '''
    dmod = distAv_kw.get('dmod', 0.)
    Av =  distAv_kw.get('Av', 0.)
    def get_extinctions(logte, logg=4.5, Rv=3.1, filt_Mag_dict={},
                        dmod=0., Av=0.):
        filename = '/Users/phil/research/Proposals/Cycle21/HPHB/data/ext_UVbright.dat'
        with open(filename,'r') as f:
            line = f.readline()
        col_keys = line.replace('#', '').strip().split()
        ext_tab = np.genfromtxt(filename, names=col_keys)
        approx_te = np.round(10 ** logte, decimals=-3)
        tab_te = np.unique(ext_tab['Teff'])
        (minds, diffs) = zip(*[rsp.math_utils.closest_match(rt, tab_te)
                               for rt in approx_te])
        inds = np.array([np.nonzero((ext_tab['Teff'] == tab_te[i]) &
                                    (ext_tab['logg'] == logg) &
                                    (ext_tab['Rv'] == Rv))[0]
                                    for i in minds]).ravel()
        filt_mag_dict = {}
        for filt, Mag in filt_Mag_dict.items():
            Al = ext_tab['A%s' % filt][inds]
            filt_mag_dict[filt] = Mag + dmod + Al * Av

        return filt_mag_dict

    pagbs, mass_inds = load_vw(mass=pagb_mass)
    pagbs = pagbs[mass_inds]
    zinds = np.nonzero(pagbs['Z'] == 0.016)[0]
    pagbs = pagbs[zinds]

    logte = pagbs['Teff']

    pagb_abs_cmd = {}
    for filt in filters:
        pagb_abs_cmd[filt] = pagbs[filt]

    pagb_cmd = get_extinctions(logte, filt_Mag_dict=pagb_abs_cmd, dmod=dmod, Av=Av)
    pagb_cmd['logte'] = logte
    pagb_cmd['logl'] = pagbs['logl']
    return pagb_cmd


if __name__ == '__main__':
    import pdb
    pdb.set_trace()
    cmd_sim_figs()