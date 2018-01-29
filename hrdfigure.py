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
    zahb['logte'] = smooth(np.array([t.data.LOG_TE[1:4]
                                     for t in ts.tracks]).ravel())
    zahb['logl'] = smooth(np.array([t.data.LOG_L[1:4]
                                    for t in ts.tracks]).ravel())

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
            inds = np.arange(rsp.math_utils.closest_match(
                3.74, track.data.LOG_L)[0])
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

    xlab = '$\log \ T_{{\\rm eff}} \ (K)$'
    ylab = '$\log \ L \ (L_\odot$)'
    ax.set_xlabel(xlab, fontsize=20)
    ax.set_ylabel(ylab, fontsize=20)
    plt.text(0.5, 1.06, '$F140LP-F336W$',
             horizontalalignment='center',
             fontsize=20,
             transform=ax.transAxes)

    ax.tick_params(labelsize=16)
    ax.set_ylim(1., 4.1)
    put_filters_top(ax, ts)
    exts = ('.png', '.eps')
    [plt.savefig('hrd%s' % ext, dpi=300, bbox_to_inches='tight')
     for ext in exts]

    return

        def check_leos_model(self, dmod=24.47, Av=0.206, af275w=2.050, af336w=1.671):
            filename = '/Users/phil/research/PHAT/UVBrightpaper/supporting_files/lf_C_Z0.07_Y0.389.dat'
            Mag, n275, n336 = np.loadtxt(filename, unpack=True)
            mF336W = Mag + dmod + Av * af336w
            mF275W = Mag + dmod + Av * af275w

            bins = np.arange(16, 32, 0.1)
            hist, bins = np.histogram(self.mF275Ws, bins=bins)
            norm = n275.max() / hist.max()
            plt.plot(bins[1:], hist * norm)

            ax = plt.axes()
            ax.plot(bins[1:], hist * norm, label='new\ sim')
            ax.plot(mF275W, n275, label='leo\ sim')
            ax.set_xlabel('F275W', fontsize=20)
            ax.legend(loc=0, numpoints=1)
            plt.savefig('comp_lf_275.png')

            # found max by hand...
            # hist[np.nonzero(np.round(bins, 1) == 26.4)]
            # np.nonzero(np.round(bins, 1) == 26.4)
            hist, bins = np.histogram(self.mF336Ws, bins=bins)
            norm = n336.max() / hist[104]

            plt.figure()
            ax = plt.axes()
            ax.plot(mF336W, n336, label='leo\ sim', color='green')
            ax.plot(bins[1:], hist * norm, label='new\ sim', color='blue')
            ax.set_xlabel('F336W', fontsize=20)
            ax.legend(loc=0, numpoints=1)
            plt.savefig('comp_lf_336.png')

    # From cycle 20
    def load_vw(H=True, mass=None):
        '''
        If H == True: gives hydrogen pagb track, else he pagb track
        if mass specified, returns the table and the indices of said mass.
        if mass not specified, returns the table and list of indices of each mass.
        '''
        names = ['mass', 'Y', 'Z', 'time', 'Teff', 'logl', 'mbol', 'F140LP',
                 'F175Wfoc', 'F275Wfoc', 'F200LP', 'F218W', 'F225W', 'F275W',
                 'F300X', 'F336W', 'F25QTZ', 'FUV', 'NUV', 'U', 'B', 'V']
        if H == True:
            filename = '/Users/phil/research/Proposals/Cycle21/HPHB/data/H-pnn.tab.UVbright'
        else:
            filename = '/Users/phil/research/Proposals/Cycle21/HPHB/data/He-pnn.tab.UVbright'

        vw_tab = np.genfromtxt(filename, names=names)

        if mass is None:
            un = np.unique(vw_tab['mass'])
            # this is the first one
            iMs = [list(vw_tab['mass']).index(u) for u in un]
            # the indices of each
            Ms = [np.arange(iMs[i], iMs[i + 1]) for i in range(len(iMs) - 1)]
            # the last one
            Ms.append(np.arange(iMs[i + 1], len(vw_tab['mass'])))
        else:
            Ms, = np.nonzero(vw_tab['mass'] == mass)
        return vw_tab, Ms

    def get_extinctions(logte, logg=4.5, Rv=3.1, filt_Mag_dict={},
                        dmod=0., Av=0.):
        filename = '/Users/phil/research/Proposals/Cycle21/HPHB/data/ext_UVbright.dat'
        with open(filename, 'r') as f:
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

    def load_pagb_dict(pagb_mass, *filters, **distAv_kw):
        '''
        Loads pagb tracks of a given mass to plot with HPHBs
        returns dict of the cmd in *filters
        '''
        dmod = distAv_kw.get('dmod', 0.)
        Av = distAv_kw.get('Av', 0.)

        pagbs, mass_inds = load_vw(mass=pagb_mass)
        pagbs = pagbs[mass_inds]
        zinds = np.nonzero(pagbs['Z'] == 0.016)[0]
        pagbs = pagbs[zinds]

        logte = pagbs['Teff']

        pagb_abs_cmd = {}
        for filt in filters:
            pagb_abs_cmd[filt] = pagbs[filt]

        pagb_cmd = get_extinctions(logte, filt_Mag_dict=pagb_abs_cmd,
                                   dmod=dmod, Av=Av)
        pagb_cmd['logte'] = logte
        pagb_cmd['logl'] = pagbs['logl']
        return pagb_cmd

    if __name__ == '__main__':
        import pdb
        pdb.set_trace()
        cmd_sim_figs()
