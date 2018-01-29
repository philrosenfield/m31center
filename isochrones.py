import os
import numpy as np
from scipy.interpolate import interp1d, splev, splrep, splprep

import re
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter
params = {'legend.fontsize': 20}
plt.rcParams.update(params)
nullfmt = NullFormatter()


center = [10.68470833,  41.26875]

def circ(r, n=100, x0=0, y0=0):
    x, y = zip(*[(np.cos(2 * np.pi / n * x) * r, np.sin(2 * np.pi / n * x) * r)
                 for x in np.arange(0,n+1)])
    return np.column_stack([np.array(x) + x0, np.array(y) + y0])

uvb = load_model('models/uvbright/parsec_v1.2s_Z0.06_Y0.356.HB.dat'

rad_arcsec = np.linspace(4.3, 300, 10)
rads = (rad_arcsec * u.arcsec).to(u.degree).value
circs = [circ(r, x0=center[0], y0=center[1]) for r in rads]

filter1s = ['F225W_VEGA', 'F275W_VEGA', 'F140LP_VEGA']
filter2 = 'F336W_VEGA'
filt2 = 'm{}'.format(filter2.replace('_VEGA', ''))


# "isochrones"
for filter1 in filter1s:
    umass = np.unique(uvb[mass])
    filt1 = 'm{}'.format(filter1.replace('_VEGA', ''))
    umass = umass[umass < .5]
    if len(umass) == 0:
        continue

    tracks = [uvb.loc[uvb[mass] == um].copy() for um in umass]
    ages = [0, 1e8, 1.4e8, 1.6e8, 1.7e8, 1.8e8, 1.9e8, 2e8]
    cols, mags = [], []

    fig, ax = plt.subplots()
    ax.plot(uvis[filter1] - uvis[filter2], uvis[filter2], '.')
    ax.set_ylim(28,20)
    ax.set_xlim(-2.5, 2)

    for track in tracks:
        #print('{:g}'.format(track[age].max()))
        interp = interp1d(track[age], track['Mode'], bounds_error=0)
        modes = np.array(np.round(interp(ages)), dtype=int)
        inds = np.concatenate([np.nonzero(track['Mode'] == m)[0] for m in modes])

        col = track[filt1].iloc[inds] - track[filt2].iloc[inds]
        mag = track[filt2].iloc[inds]
        ax.plot(track[filt1] - track[filt2], track[filt2],
                label=r'$M={:.3f}M\odot$'.format(track[mass].iloc[0]))
        # ax.plot(col, mag , 'o')
        dif = len(ages) - len(inds)
        mag[mag > 28] = np.nan
        if dif > 0:
            col = np.append(col, np.zeros(dif) * np.nan)
            mag = np.append(mag, np.zeros(dif) * np.nan)

        cols.append(col)
        mags.append(mag)
    # ax.legend()
    col_seq = np.array(cols).T
    mag_seq = np.array(mags).T

    for i, _ in enumerate(col_seq):
        m = np.argsort(col_seq[i])
        lw = 1
        if i == 0:
            lw = 2
        elif i == len(col_seq) - 1:
            lw = 2
        else:
            color = 'grey'
        ax.plot(col_seq[i][m], mag_seq[i][m], color='k', lw=lw)
        try:
            j = np.nanargmin(col_seq[i])
        except:
            continue
        ax.text(col_seq[i][j] - np.abs(col_seq[i][j]) * 0.01, mag_seq[i][j],
                '{:g}'.format(ages[i]/1e8), ha='right', fontsize=8)
    ax.set_title('{:g}'.format(track[z].iloc[0]))
    # ax.plot(track[logt], track[logl], label=track[mass].iloc[0])
        # ax.plot(track[logt].iloc[inds], track[logl].iloc[inds], 'o')


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


def load_extinction_table(extfn=None):
    extfn = extfn or 'models/ext_UVbright.dat'
    return np.genfromtxt(extfn, names=True)

class InterpolateTracks(object):
   def __init__(self, ts_indict, max_mass=120., min_mass=0., mass_array=None):
       # load tracks?
       self.ext_tab = load_extinction_table(extfn)
       # self.cut_masses(max_mass=max_mass, min_mass=min_mass, mass_array=mass_array)


   def cut_masses(self, max_mass=120., min_mass=0., mass_array=None):
       if mass_array is None:
           self.tracks = [t for t in self.tracks if t.mass >= min_mass and
                                                    t.mass =< max_mass]
       else:
           self.tracks = [t for t in self.tracks if '%.3f' % t.mass in ['%.3f' % m for m in mass_array]]

   def cut_wd(self):
       '''
       adds attribute to track with inds of non-wd.
       '''
       for t in self.tracks:
           logg = -10.616 + np.log10(t.mass) + 4.0 * t.data.LOG_TE - t.data.LOG_L
           t.not_wd, = np.nonzero(logg < 6.5)

   def get_extinctions(self, logte, logg=4.5, Rv=3.1, filters=None):
       A2 = 'no data'
       approx_te = np.round(10**logte, decimals=-3)
       table_te = np.unique(self.ext_tab['Teff'])
       (matched_inds, diffs) = zip(*[rsp.math_utils.closest_match(rt, table_te) for rt in approx_te])
       inds = np.array([np.nonzero((self.ext_tab['Teff'] == table_te[i]) &
                                   (self.ext_tab['logg'] == logg) &
                                   (self.ext_tab['Rv'] == Rv))[0]
                                   for i in matched_inds]).ravel()

       for filt in filters:
           try:
               self.__setattr__('A%s' % filt, self.ext_tab['A%s' % filt][inds])
           except IndexError:
               print '%s not found.' % filt

   def recovered_stars(self, filter1, filter2, mlim=90):
       m1 = self.__dict__['m%ss_cor' % filter1]
       m2 = self.__dict__['m%ss_cor' % filter2]
       recs, = np.nonzero((m1 < mlim) & (m2 < mlim))
       return recs

   def correct_extinction_dist(self, track, logte, filt_Mag_dict, dmod=None, Av=None):
       self.get_extinctions(logte, filters=filt_Mag_dict.keys())
       for filt, Mag in filt_Mag_dict.items():
           Al = self.__dict__['A%s' % filt]
           track.__setattr__('m%s' % filt, Mag + dmod + Al * Av)

   def squish(self, *attrs, **kwargs):
       '''
       concatenates an attribute or many attributes and adds them to galaxies
       instance -- with an 's' at the end to pluralize them... that might
       be stupid.
       ex
       t.squish('mF275W', 'mag2', 'ra', 'dec')
       gs.ras ...
       '''
       slice_attr = kwargs.get('slice_attr')
       if slice_attr is not None:
           for attr in attrs:
               try:
                   new_list = [t.__getattribute__(attr)[t.__dict__[slice_attr]] for t in self.tracks]
               except AttributeError:
                   new_list = [t.data[attr][t.__dict__[slice_attr]] for t in self.tracks]
               new_val = np.concatenate(new_list)
               self.__setattr__('%ss' % attr, new_val)
       else:
           for attr in attrs:
               try:
                   new_list = [t.__getattribute__(attr) for t in self.tracks]
               except AttributeError:
                   new_list = [t.data[attr] for t in self.tracks]
               new_val = np.concatenate(new_list)
               self.__setattr__('%ss' % attr, new_val)


   def interpolate_from_mode(self, *cols):
       for track in self.tracks:
           mode = track.data.MODE
           for col in cols:
               track.__setattr__('%s_interp' % col, interp1d(mode,
                                                             track.data[col],
                                                             bounds_error=0))
       return

   def simulate_cmd(self, filters=['F275W', 'F336W', 'F225W', 'F140LP'],
                    Av=None, dmod=None, delta_t=5e4):
       '''
       interpolates tracks in fixed time steps, delta_t
       creates interpolated cmd as self.m<filter> or self.M<filter> if no
       dmod and Av supplied.
       '''
       tot_mass = 0.

       self.fixed_dtime_interpolation(delta_t=delta_t)
       self.interpolate_from_mode('LOG_TE', *filters)

       for track in self.tracks:
           fofte = track.LOG_TE_interp(track.model_arr)
           filt_mag_dict = {}
           for filt in filters:
               filt_mag_dict[filt] = track.__dict__['%s_interp' % filt](track.model_arr)

           if Av is not None or dmod is not None:
               self.correct_extinction_dist(track, fofte, filt_mag_dict, Av=Av,
                                            dmod=dmod)
           else:
               [track.__setattr__('M%s' % filt, Mag) for filt, Mag in
                                filt_Mag_dict.items()]

           tot_mass += float(len(filt_mag_dict[filt])) * track.mass
           track.__setattr__('mass_interp', np.repeat(track.mass, len(filt_mag_dict[filt])))

       if Av is not None or dmod is not None:
           self.squish(*['m%s' % f for f in filters])
       else:
           self.squish(*['M%s' % f for f in filters])

       self.squish('mass_interp')
       return tot_mass
