import os
import glob
import numpy as np
import matplotlib.pyplot as plt


def parse_mag_tab(photsys, filt, bcdir=None, all=False):
    if not bcdir:
        try:
            bcdir = os.environ['BCDIR']
        except KeyError:
            print('Error: need to pass BCDIR or exist as environmental variable')
            return 1

    tab_mag_dir = os.path.join(bcdir, 'tab_mag_odfnew/')
    tab_mag, = glob.glob1(tab_mag_dir, 'tab_mag_{:s}.dat'.format(photsys))

    tab = open(os.path.join(tab_mag_dir, tab_mag), 'r').readlines()
    mags = tab[1].strip().split()
    Alam_Avs = np.array(tab[3].strip().split(), dtype=float)
    Alam_Av = {k:v for k,v in zip(mags, Alam_Avs)}
    if not all:
        retv = Alam_Av[filt]
    else:
        retv = Alam_Av
    return retv

Alam_Av = {'F225W': 2.32794, 'F336W': 1.65798, 'F275W': 1.94440,
           'F140LP': 1.82410}

dmod = 24.47
Av = 0.206

fname = 'Z0.03Y0.302/Z0.03Y0.302OUTA1.74_F7_M0.467.HB.DAT.UVbright'
filt1 = 'F225W'
filt2 = 'F336W'
for fname in fnames:
    data =  np.genfromtxt(fname, names=True)
    mag1 = data[filt1] + dmod + Av * Alam_Av[filt1]
    mag2 = data[filt2] + dmod + Av * Alam_Av[filt2]
    color = mag1 - mag2
    mag = mag2
    ax.plot(color, mag, alpha=0.3)
