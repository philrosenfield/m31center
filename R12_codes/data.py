import os
import numpy as np
import graphics
import matplotlib.pyplot as plt
import GraphicsUtils
import file_io

def reg_to_array(reg_file,commentchar='#',tag=None):
    import re
    ra,dec,tags = np.array([]),np.array([]),np.array([])
    for r in open(reg_file,'r').readlines():
        r = r.strip()
        if r.startswith(commentchar): continue
        if r.startswith('circle'):
            ra = np.append(ra,float(r.split(',')[0].split('(')[-1]))
            dec = np.append(dec,float(r.split(',')[1]))
            tags = np.append(tag,str(r.split('{')[-1].replace('}','')))
        if r.startswith('polygon'):
            if tag==None:
                data = map(float,r.split('(')[-1].split(')')[0].split(','))
                ra = np.append(ra,data[::2])
                dec = np.append(dec,data[1::2])
                tags = np.append(tag,str(r.split('{')[-1].replace('}','')))
            else:
                if re.search(tag,r):
                    data = map(float,r.split('(')[-1].split(')')[0].split(','))
                    ra = data[::2]
                    dec = data[1::2]
    return ra,dec,tags
    
def quality_cut(sharp,snr,crowd,round,loud=False,**kwargs):
    sharp_cut = kwargs.get('sharp_cut',0.075)
    snr_cut = kwargs.get('snr_cut',4.)
    crowd_cut = kwargs.get('crowd_cut','None')
    round_cut = kwargs.get('round_cut',1.5)
    if loud == True:
        print '# abs(sharp)<%f snr>%f abs(round)<%f crowd<%s'%(sharp_cut,
                                                              snr_cut,
                                                              round_cut,
                                                              str(crowd_cut))
    
    snrc = np.nonzero(snr>snr_cut)[0]
    sharpc = np.nonzero(abs(sharp)<sharp_cut)[0]
    roundc = np.nonzero(abs(round)<round_cut)[0]
    
    if crowd_cut != 'None':
        crowdc = np.nonzero(abs(crowd)<crowd_cut)[0]
        # crowd could be 0.1
        c = list(set(snrc) & set(sharpc) & set(crowdc) & set(roundc))
    else:
        c = list(set(snrc) & set(sharpc) & set(roundc))
    
    return c

def quality_cut_double_st(sharp1,snr1,sharp2,snr2,defaults=True):
    if defaults == True:
        snr_cut = 4.
        sharp_cut = 0.075
        snrc1 = np.nonzero(snr1>snr_cut)[0]
        snrc2 = np.nonzero(snr2>snr_cut)[0]
        sharpc = np.nonzero((sharp1+sharp2)**2<sharp_cut)[0]
    else:
        'hard code mutherfucker.'
    return list(set(snrc1) & set(sharpc) & set(snrc2))
    
def do_quality_cuts(phot,**kwargs):
    '''
    single detections are combined to get color cut.
    '''
    single_det1 =  quality_cut(phot['sharp1'],
                               phot['snr1'],
                               phot['crowd1'],
                               phot['round1'],
                               **kwargs)
    single_det2 =  quality_cut(phot['sharp2'],
                               phot['snr2'],
                               phot['crowd2'],
                               phot['round2'],
                               **kwargs)
    
    both = list(set(single_det2) & set(single_det1))
    return single_det1, single_det2, both

def make_dolphot_residual_plot(phot, both, rad_inds, phot_cat='',**kwargs):
    
    all_inds = range(0,len(phot['mag1']))
    rejects = list(set(all_inds) - set(both))
    areas = kwargs.get('areas')
    
    sfx = '_rj'
    stitle = 'rejects'
    
    cmdkwargs = {'filter1': kwargs.get('Filter1'),
                 'filter2': kwargs.get('Filter2'),
                 'abs_mag': False,
                 'save_plot': True,
                 'xlim': (-1,2.5),
                 'ylim': (26,19)
                 }

    print 'region rejects/arcsec^2 rejects/sources in region rejects in cmd stars in cmd' 
    for i in range(len(rad_inds)):
        inds = list(set(rejects) & set(rad_inds[i]))
        nrejects = float(len(inds))
        nsources = float(len(rad_inds[i]))
        color = phot['mag1'] - phot['mag2']
        mag2 = phot['mag2']
        xlim = cmdkwargs['xlim']
        ylim = cmdkwargs['ylim']
        sources_in_cmd = in_cmd_box(color[rad_inds[i]],mag2[rad_inds[i]],
                                    xlim,ylim)
        rej_in_cmd = in_cmd_box(color[inds],mag2[inds],xlim,ylim)
        
        n_scmd = float(len(sources_in_cmd))
        n_rcmd = float(len(rej_in_cmd))
        print i, n_rcmd/areas[i], n_rcmd/n_scmd, n_rcmd, n_scmd
        cmdkwargs['mag1'] = phot['mag1'][inds]
        cmdkwargs['mag2'] = phot['mag2'][inds]
        cmdkwargs['title'] = phot_cat+' '+stitle,
        cmdkwargs['outfile'] = 'cmd_rejects_region%i.png'%i
        print cmdkwargs['outfile']
        ax = GraphicsUtils.plot_cmd(None,**cmdkwargs)

    return

def make_reg_of_shit(phot,both,**kwargs):
    shit = 0.0
    crap = np.nonzero(phot['sharp1'][both] < -0.5)[0]
    PHATDataUtils.radec_to_reg(phot['ra'][both][crap],phot['dec'][both][crap],outfile='crapsharp.reg')

    crap = np.nonzero(phot['sharp1'] < -0.5)[0]
    PHATDataUtils.radec_to_reg(phot['ra'][crap],phot['dec'][crap],outfile='crapsharp.reg')
   
def in_cmd_box(color,mag2,xlim,ylim):
    return np.nonzero( (color>xlim[0]) & (color < xlim[1]) & (mag2<ylim[0]) & (mag2>ylim[1]))[0]

def make_dolphot_cut_plots(phot,sd1, sd2, both, comp, good, rad_inds,
                          phot_cat='',**kwargs):
    '''
    makes plots of x vs mag1 mag2
    makes plots of x vs x
    for x= sharp, crowd, round
    
    phot_cat is an optional title string and prefix to saved plots.
    
    '''
    cols = kwargs.get('cols')
    phot_dir = kwargs.get('phot_dir','.')
    startdir = os.getcwd()
    os.chdir(phot_dir)
    indstoplot = [(sd1,sd2),(both,both),(good,good)]
    sfx = ['_sd','_dd','_good']
    stitle = ['single detections','color cut with both detections',
              'color, mag, and red cut']
    dolphots = ['sharp','crowd','round']
    lims = [(-1,1), (0,3.5),(-1.5,6.5)]
    f1 = 'F275W'
    f2 = 'F336W'
    xlim = (19,28)
    if len(phot_cat)>0: pc = phot_cat+'_'
    for i,inds in enumerate([(sd1,sd2),(both,both),(good,good)]):
        for k,dolphot in enumerate(dolphots):
            y1 = phot['%s1'%dolphot]
            y2 = phot['%s2'%dolphot]
            x1 = phot['mag1']
            x2 = phot['mag2']
            ylab = dolphot.title()
            ylim = lims[k]
            ax1,ax2 = graphics.make_double_plot(x1,y1,x2,y2,inds[0],inds[1],
                                                rad_inds,f1,f2,ylab,xlim,ylim,
                                                cols,fig=k+1)
            [ax.set_xlim(xlim) for ax in (ax1,ax2)]
            ax1.set_title(phot_cat+' '+stitle[i])
            plt.savefig('%s%s_mag%s.png'%(pc,dolphot,sfx[i]))
            
        plt.close('all')
        # these plots don't make sense for single detections.
        if inds[0] != inds[1]: continue 
        for k,dolphot in enumerate(dolphots):
            x = phot['%s1'%dolphot]
            y = phot['%s2'%dolphot]
            ylim = lims[k]
            ylab = dolphot.title()
            ax = graphics.make_single_plot(x,y,inds[0],rad_inds,f1,f2,
                                           ylab,ylim,ylim,cols,fig=k+1)
            ax.set_title(phot_cat+' '+stitle[i])
            plt.savefig('%s%s%s.png'%(pc,dolphot,sfx[i]))
        plt.close('all')
    print 'wrote photometry diagnostic plots to %s'%phot_dir
    os.chdir(startdir)
    return

def photometry_cuts(phot,**kwargs):
    '''
    phot is the full catalog dictionary, no cuts. 
    returns indicies of phot as follows:
    
    For Table 1, and integrated LF:
    1. Single detections in each filter: sd1, sd2
    
    For CMDs which make the HP-HB/PE-AGB separation, which are then used to make
    LFs of HP-HB stars, 
    the Surface Brightness as a function of radius, 
    the HP-HB/PEAGB ratio,
    the surface density as a function of radius
    the fuel consumption thereom.
    
    2. Single detection in both filters combined to make a color cut for CMDs: both
    2a. 2. + Corrected for 50% completeness limit (in both filters): comp
    2b. 2a. + Corrected for red stars: good
    '''
    
    single = kwargs.get('single')
    
    mag1_cut = kwargs.get('mag1_cut',99.)
    mag2_cut = kwargs.get('mag2_cut',99.)
    
    if mag2_cut == 99. or mag1_cut == 99.: print 'no mag cuts.'
    
    comp1 = np.nonzero(phot['mag1']<mag1_cut)[0]
    comp2 = np.nonzero(phot['mag2']<mag2_cut)[0]

    if single == True:
        sd1, sd2, both = do_quality_cuts(phot,**kwargs)
        comp = list(set(comp1) & set(comp2) & set(both))
    else:
        comp = list(set(comp1) & set(comp2))
        both = comp
        sd1,sd2 = -1,-1
    red = np.nonzero((phot['mag1']-phot['mag2']>0.8) &
                     (phot['mag2']<19))[0]
    print 'red stars cut out:',len(red)
    good = list(set(comp)-set(red)) 
    return sd1, sd2, both, comp, good

def completeness_cut(phot,**kwargs):
    mag1_cut = kwargs.get('mag1_cut',99.)
    mag2_cut = kwargs.get('mag2_cut',99.)
    comp1 = np.nonzero(phot['mag1']<mag1_cut)[0]
    comp2 = np.nonzero(phot['mag2']<mag2_cut)[0]
    comp =  list(set(comp1) & set(comp2))
    red = np.nonzero((phot['mag1']-phot['mag2']>0.8) &
                     (phot['mag2']<19))[0]
    good = list(set(comp)-set(red))
    both = comp
    
    return comp1,comp2,both,comp,good

def get_mag_cuts(**kwargs):
    comp_data = kwargs['comp_data']
    mag_cut = kwargs.get('mag_cut',90)
    if mag_cut == 50:
        print 'using 50% completeness mag cut'
        mag1_cut,mag2_cut =  file_io.get_comp50mags(comp_data)[0][0],file_io.get_comp50mags(comp_data)[1][0]
    elif mag_cut == 90:
        print 'using 90% completeness mag cut'
        mag1_cut,mag2_cut =  file_io.get_comp90mags(comp_data)[0][0],file_io.get_comp90mags(comp_data)[1][0]
    else:
        print 'no mag cut'
        mag1_cut, mag2_cut = None, None
    print 'mag1_cut: ',mag1_cut
    print 'mag2_cut: ',mag2_cut    
    return mag1_cut, mag2_cutimport PadovaTracksUtils
import numpy as np
import GenUtils
import os

def info_from_track_filename(filename):
    '''
    ex
    info_from_track_filename('/astro/net/angst2/philrose/tracks/C09_AF02_I1_Z0.07_Y0.375/Z0.07Y0.375OUTA1.74_F7_M.485.PMS.HB.res0.00.UVbright')
    {'F': 7.0,
     'M': 0.48499999999999999,
     'OUTA': 1.74,
     'Y': 0.375,
     'Z': 0.070000000000000007,
     'photsys': 'UVbright',
     'res': 0.0}
    '''
    d = {}
    
    # remove dir info
    filename=filename.split('/')[-1]
    
    # split on PMS (tail is usually color_tracks info, but can be just HB)
    a,z = filename.split('.PMS')
    
    # these two will be nonsense if color_tracks.py hasn't been run:
    d['res'] = float('.'.join(z.split('res')[-1].split('.')[0:2]))
    d['photsys'] = z.split('res')[-1].split('.')[-1]
    
    # take the _ of the filename and give spaces to split
    a = a.replace('_',' ').replace('R1','')
    
    # get all the letters together for dict keys:
    s = ''.join(c for c in a if not c.isdigit())
    s = s.replace('.',' ').split()
    
    # copy the filename
    x = a[:]
    s.append(' ')
    # for each dict key
    for i in range(len(s)-1):
        # if dict key is in the filename
        if s[i] in a:
            # remove the key from filename copy
            x = x.replace(s[i],'')
            # grab the value next to the dict key
            y = float(x.split(s[i+1])[0])
            # update x for the next value
            x = x.split(s[i+1])[-1]
            d[s[i]]= y
    return d    

def load_track(track_file):
    f=open(track_file,'r')
    lines = f.readlines()
    f.close()
    
    for i in range(len(lines)):
        if lines[i].strip().startswith('BEGIN TRACK'): break
    try:
        col_keys = lines[i+1].split()
    except UnboundLocalError:
        print 'No BEGIN TRACK ', filename
        
    
    track = np.genfromtxt(track_file,names=col_keys,skip_header=i+2)
    name = os.path.split(track_file)[1]
    info = info_from_track_filename(track_file)
    Track = PadovaTrack(name,
                        track,
                        info['M'],
                        info['Z'],
                        info['Y'],
                        info['res'],
                        info['photsys'])
    return Track

def load_tracks(track_dir,track_ext,mass=None):
    track_files = GenUtils.get_afile(track_dir,'*'+track_ext)
    tracks = []
    for track_file in track_files:
        tracks.append( load_track(track_file))
    if mass != None:
        tracks = [track for track in tracks if track.mass == mass]
    return tracks

def get_tracks_col(tracks,col):
    return np.array([t.track[col] for t in tracks])

class PadovaTrack(object):
    def __init__(self,
                 name,
                 track,
                 mass,
                 Z,
                 Y,
                 res,
                 photsys):
        self.track = track
        self.name = name
        self.mass = mass
        self.Z = Z
        self.Y = Y
        self.res = res
        self.photsys = photsys
    
    def get_stable_col(self,key):
        lnuc = self.track['LNUC']
        ind = np.nonzero((lnuc < 1.2) & (lnuc > 0.8))[0]
        try:
            col = self.track[key]
            scol = col[ind]
        except KeyError:
            print key,'not found'
            scol = 0
        return scol
        
    def get_shell_heb_col(self,key):
        ycen = self.track['YCEN']
        ly = self.track['LY']
        lnuc = self.track['LNUC']
        ind = np.nonzero((ycen==0) & (ly>0) & (lnuc < 1.2) & (lnuc > 0.8))[0]
        try:
            col = self.track[key]
            scol = col[ind]
        except KeyError:
            print key,'not found'
            scol = 0
        return scol

def load_vw(H=True,mass=None):
    '''
    If H == True: gives hydrogen pagb track, else he pagb track
    if mass specified, returns the table and the indices of said mass.
    if mass not specified, returns the table and list of indices of each mass.
    '''
    names = ['mass','Y', 'Z', 'time', 'Teff', 'logl','mbol', 'F140LP', 'F175Wfoc', 'F275Wfoc', 'F200LP', 'F218W', 'F225W', 'F275W', 'F300X', 'F336W', 'F25QTZ', 'FUV', 'NUV', 'U','B','V']
    if H==True:
        filename = '/astro/net/angst2/philrose/tracks/vw/H-pnn.tab.UVbright'
    else:
        filename = '/astro/net/angst2/philrose/tracks/vw/He-pnn.tab.UVbright'
    vw_tab = np.genfromtxt(filename,names=names)

    if mass == None:
        un = np.unique(vw_tab['mass'])
        iMs = [list(vw_tab['mass']).index(u) for u in un] # this is the first one
        Ms = [np.arange(iMs[i],iMs[i+1]) for i in range(len(iMs)-1)] # The indices of each
        Ms.append(np.arange(iMs[i+1],len(vw_tab['mass'])))
    else:
        Ms = np.nonzero(vw_tab['mass']==mass)[0]
    return vw_tab,Ms

def load_ext_table(filename='ext_UVbright.dat'):
    colheads = open(filename,'r').readline().replace('#','').strip().split()
    ext_tab = np.genfromtxt(filename,names=colheads)
    return ext_tab

def load_tracks_old(track_location,track_ext):
    '''
    returns a list of PadovaTracksUtils.Track objects
    for all files that match track_ext in directory track_location.
    ''' 
    track_files = GenUtils.get_afile(track_location,'*'+track_ext)
    tracks = [PadovaTracksUtils.get_tracks(t) for t in track_files]
    return tracks

def load_pagb_plot_kwargs(pagb_mass,filt1,filt2,filt3,filt4,age_fact=50.,**kwargs):
    '''
    Loads pagb tracks to plot with HPHBs
    delta_t is time spacing to plot
    age_fact is some multiple of detla_t to increase the sampling. 
    The idear is delta_t is same for HPHBs, but is changed by age_fact for pagbs.
    This is really only set up for one mass, and if you want the legend to work, 
    need to set masses by hand.
    
    returns dict of the cmd in four filters, and the kwargs needed in
    plot_interp_cmd_once.
    '''
    Av = kwargs.get('Av')
    dmod = kwargs.get('dmod')
    delta_t = kwargs.get('delta_t')
    pagbs,mass_inds = load_vw(mass=pagb_mass)
    pagbs = pagbs[mass_inds]
    zinds = np.nonzero(pagbs['Z'] == 0.016)[0]
    pagbs = pagbs[zinds]
    # model number is going to be a uniform array.
    models = np.arange(len(pagbs))
    # masses in the plot code are used for legend and sort plot.
    if pagb_mass == 1.: masses = 0.569
    if pagb_mass == 1.5: masses = 0.597
    # hot chooses which to plot in color... this is all.
    ages = pagbs['time']
    line_colors = 'black'
    teffs = pagbs['Teff']
    Mag1 = pagbs[filt1]
    Mag2 = pagbs[filt2]
    Mag3 = pagbs[filt3]
    Mag4 = pagbs[filt4]
    pagb_cmd = { filt1: Mag1, filt2: Mag2, filt3: Mag3, filt4: Mag4}
    pagb_kwargs = { 'Av': Av,
                    'dmod': dmod,
                    'ages': ages,
                    'masses': masses,
                    'teffs': teffs,
                    'line_colors': line_colors,
                    'models': models,
                    'delta_t': delta_t,
                    'age_fact': age_fact}

    return pagb_cmd, pagb_kwargs

from astroML.stats import binned_statistic_2d
        for j in range(len(types)):
            for i, mass_met in enumerate([opt_mass_met_file,
                                          ir_mass_met_file]):
                k += 1
                with open(mass_met, 'r') as mmf:
                    lines = [l.strip() for l in mmf.readlines()
                             if not l.startswith('#')]

                mag = np.concatenate([np.array(l.split(), dtype=float)
                                      for l in lines[0::3]])
                mass = np.concatenate([np.array(l.split(), dtype=float)
                                       for l in lines[1::3]])
                mh = np.concatenate([np.array(l.split(), dtype=float)
                                     for l in lines[2::3]])

                N, xedges, yedges = binned_statistic_2d(mag, mass, mh,
                                                        types[j], bins=50)
                im = grid[k].imshow(N.T, origin='lower',
                               extent=[xedges[0], xedges[-1], yedges[0],
                                       yedges[-1]],
                               aspect='auto', interpolation='nearest',
                               cmap=cmaps[j])
                grid[k].cax.colorbar(im)
                #grid[i].cax.set_label('$[M/H]$')

