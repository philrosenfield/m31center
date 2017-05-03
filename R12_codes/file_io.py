import os
import pyfits
from GenUtils import get_afile
import numpy as np
from PadovaTracksUtils import get_tracks
import matplotlib.nxutils as nxutils
from params import *
def load_integrated_phot():
    from params import PHAT_DIR_EXTRA2
    filename = PHAT_DIR_EXTRA2+'phatuv_intphot_v2.dat'
    names = ['f275w_cts',
             'f275w_mag',
             'f275w_sb',
             'f336w_cts',
             'f336w_mag',
             'f336w_sb',
             'int_color']
    int_flx = np.genfromtxt(filename,names=names)
    return int_flx

def load_contamination(filename):
    lines = open(filename).readlines()    
    d = {}
    for line in lines:
        if line.startswith('#'): continue
        line = line.strip()
        filter,region,val = line.split()
        if not filter in d.keys(): d[filter] = {}
        d[filter][region] = float(val)
    return d
    
def load_defaults(filename):
    try:
        literal_eval
    except NameError:
        from ast import literal_eval
    lines = open(filename).readlines()
    d = {}
    for line in lines:
        if line.startswith('#'): continue
        if len(line.strip()) == 0: continue
        key,val = line.strip().split()
        try:
            d[key] = float(val)
        except:
            if val.startswith('UVBrights'): val = os.path.join(PHAT_CODE_DIR,val)
            if val.upper().startswith('TR') or val.upper().startswith('FA'): val = literal_eval(val)
            d[key] = val
    return d

def load_hphb_model(**kwargs):
    from params import PHAT_TRACKS_DIR
    dmod = kwargs.get('dmod')
    Av = kwargs.get('Av_bulge')
    af275w = kwargs.get('af275w')
    af336w = kwargs.get('af336w')
    filename = PHAT_TRACKS_DIR+'lf_C_Z0.07_Y0.389.dat'
    model_absMag, n275, n336 = np.loadtxt(filename,unpack=True)
    model_F336W = model_absMag+dmod+Av*af336w
    model_F275W = model_absMag+dmod+Av*af275w
    return model_F336W, n336, model_F275W, n275

def load_comp_data(filename):
    names = ['region', 'comp50mag1','comp50mag2', 'comp90mag1','comp90mag2']
    data = np.genfromtxt(filename,names=names)
    return data
    
def get_comp50mags(filename):
    region, comp50mag1,comp50mag2, comp90mag1,comp90mag2 = np.loadtxt(filename,unpack=True)        
    return comp50mag1,comp50mag2

def get_comp90mags(filename):
    region, comp50mag1,comp50mag2, comp90mag1,comp90mag2 = np.loadtxt(filename,unpack=True)        
    return comp90mag1,comp90mag2

def read_m31_bulge_phot(**kwargs):
    from params import PHAT_DIR_EXTRA
    names = ['region',
             'semimajor_axis',
             'sb_36',
             'rif_36',
             'rasb_36',
             'rif_fuv',
             'rasb_fuv',
             'rif_nuv',
             'rasb_nuv',
             'fuv_nuv_color',
             'nPNe',
             'bulge_fraction']
    dat = np.genfromtxt(PHAT_DIR_EXTRA+'m31bulge_phot_v4.dat',names=names)
    
    overlap = kwargs.get('overlap')
    if overlap == True:
        # see overlap_fields.get_overlap_areas()[1]
        frac_areas = np.array([ 1., 0.98457183,  0.84610758,  0.29658597,  0.16800527,  0.10820975])
        dat['rif_36'] = dat['rif_36']*frac_areas
        print 'changing radial intergated flux to scale to overlap areas'
    return dat

def load_completeness_corrections(filename):
    frac1s = frac2s = []
    lines = open(filename,'r').readlines()
    magbins = map(float,lines[1].strip().split())
    for line in lines[2:]:
        if line.startswith('F336W'):
            dat = line.strip().replace('F336W','').split()
            region = int(dat[0])
            frac2s.append(np.array(map(float,dat[1:])))
        if line.startswith('F275W'):
            dat = line.strip().replace('F275W','').split()
            region = int(dat[0])
            frac1s.append(np.array(map(float,dat[1:])))
    return frac1s,frac2s,magbins

def load_all_tracks(file_ext='.phat',loud=False):
    '''
    loads as lists of track objects.
    '''
    from params import PHAT_TRACKS_BASE,VW_BASE,TRACKS_BASE
    track_bulge = PHAT_TRACKS_BASE+'hbz07y389/'
    btrack = [get_tracks(phat,loud=loud) for phat in get_afile(track_bulge,'*HB'+file_ext)]
    
    track_bulge_more = PHAT_TRACKS_BASE+'hbz07y389/toomanytoplot/'
    btrackmore = [get_tracks(phat,loud=loud) for phat in get_afile(track_bulge_more,'*'+file_ext)]

    track_disk = PHAT_TRACKS_BASE+'S11_Z0.017_Y0.279/'
    dtrack = [get_tracks(phat,loud=loud) for phat in get_afile(track_disk,'*'+file_ext)]
    
    # The HB tracks don't have MS...
    track_msbulge = TRACKS_BASE+'ell00/C_Z0.07_Y_0.389/'
    bmstrack = [get_tracks(phat,loud=loud) for phat in get_afile(track_msbulge,'*ENV_0.7*'+file_ext)]
    
    return btrack, btrackmore, dtrack, bmstrack

def read_reg(reg_name,shape='polygon'):
    '''
    Takes a ds9 reg file and loads ra dec into arrays. 
    Only tested on polygon shape. 
    returns ras, decs
    '''
    xys = [map(float,line.replace(shape+'(','').replace(')','').split('#')[0].split(',')) 
               for line in open(reg_name,'r').readlines() if line.startswith(shape)]
    ras = [xy[0::2] for xy in xys]
    decs =[xy[1::2] for xy in xys]
    return ras,decs

def get_bulge_file(single=True,overlap=False,st=False):
    from params import PHAT_DATA_BASE
    if single==True:
        bulge_fit = PHAT_DATA_BASE+'uv-brick01-diy-v4st.fits'
    if overlap == True:
        bulge_fit = '/astro/net/angst2/philrose/PHAT/misc/overlap_test/12058_M31-B01-F040510111617-UVIS_F275W_F336W.st.fits'
    if st == True:
        bulge_fit = PHAT_DATA_BASE + 'uv-brick01-v4st.fits'
    return bulge_fit

def read_data_diy(bulge_fit_sd):
    """
    For dustin's DIY catalog.
    """
    # Single Detections
    print 'using',bulge_fit_sd.split('/')[-1],'for single detections'
    b1_sd = pyfits.open(bulge_fit_sd)[1].data
    
    primary = b1_sd.field('primary')
    non_dupes = np.nonzero(primary)[0]
    mag1  = b1_sd.field('mag1_uvis')[non_dupes]
    mag2  = b1_sd.field('mag2_uvis')[non_dupes]
    snr1    = b1_sd.field('snr1')[non_dupes]
    snr2    = b1_sd.field('snr2')[non_dupes]
    sharp1  = b1_sd.field('sharp1')[non_dupes]
    sharp2  = b1_sd.field('sharp2')[non_dupes]
    crowd1  = b1_sd.field('crowd1')[non_dupes]
    crowd2  = b1_sd.field('crowd2')[non_dupes]
    round1  = b1_sd.field('round1')[non_dupes]
    round2  = b1_sd.field('round2')[non_dupes]
    mag1err = b1_sd.field('mag1_err')[non_dupes]
    mag2err = b1_sd.field('mag2_err')[non_dupes]
    ra    = b1_sd.field('ra')[non_dupes]
    dec   = b1_sd.field('dec')[non_dupes]
    
    phot = {'mag1': mag1,
            'mag2': mag2,
            'snr1': snr1,
            'snr2': snr2,
            'sharp1': sharp1,
            'sharp2': sharp2,
            'crowd1': crowd1,
            'crowd2': crowd2,
            'round1': round1,
            'round2': round2,
            'mag1err': mag1err,
            'mag2err': mag2err,
            'ra': ra,
            'dec': dec}
    return phot

def read_overlap_data():
    from overlap_fields import get_overlap_regions
    fits_file = os.path.join(PHAT_DATA_BASE,'12058_M31-B01-F040510111617-UVIS_F275W_F336W.st.fits')
    print 'using',fits_file,'for overlap region analysis'
    ovlp = pyfits.open(fits_file)[1].data
    
    mag1  = ovlp.field('mag1_uvis')
    mag2  = ovlp.field('mag2_uvis')
    snr1    = ovlp.field('snr1')
    snr2    = ovlp.field('snr2')
    sharp1  = ovlp.field('sharp1')
    sharp2  = ovlp.field('sharp2')
    crowd1  = ovlp.field('crowd1')
    crowd2  = ovlp.field('crowd2')
    round1  = ovlp.field('round1')
    round2  = ovlp.field('round2')
    mag1err = ovlp.field('mag1_err')
    mag2err = ovlp.field('mag2_err')
    ra    = ovlp.field('ra')
    dec   = ovlp.field('dec')
    
    radec = np.column_stack((ra,dec))
    verts = get_overlap_regions()[-1]
    inds = np.nonzero(nxutils.points_inside_poly(radec, verts))[0]
    
    phot = {'mag1': mag1[inds],
            'mag2': mag2[inds],
            'snr1': snr1[inds],
            'snr2': snr2[inds],
            'sharp1': sharp1[inds],
            'sharp2': sharp2[inds],
            'crowd1': crowd1[inds],
            'crowd2': crowd2[inds],
            'round1': round1[inds],
            'round2': round2[inds],
            'mag1err': mag1err[inds],
            'mag2err': mag2err[inds],
            'ra': ra[inds],
            'dec': dec[inds]}
    return phot
     
def read_data_st(bulge_fit):
    # Read Data
    b1 =  pyfits.open(bulge_fit)[1].data
    print 'using',bulge_fit.split('/')[-1]
    try:
        mag1 = b1.field('mag1')
        mag2 = b1.field('mag2') 
        ra = b1.field('ra')
        dec = b1.field('dec')
    except KeyError:
        mag1 = b1.field('MAG1_UVIS')
        mag2 = b1.field('MAG2_UVIS') 
        ra = b1.field('RA')
        dec = b1.field('DEC')

    phot = {'mag1': mag1,
            'mag2': mag2,
            'ra': ra,
            'dec': dec}
    return phot
