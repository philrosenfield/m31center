import os
import GenUtils 
import Astronomy
import file_io
import graphics
import PadovaTracksUtils

from matplotlib.ticker import NullFormatter, MultipleLocator,ScalarFormatter
import pyfits
import numpy as np
import matplotlib.nxutils as nxutils
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text',usetex=True)
from params import *

def area_within_regs(reg_name,closed=False):
    '''
    if closed = false, I'll append the first x,y to the list.
    '''
    tot_areas = np.array([])
    ras,decs = file_io.read_reg(reg_name,shape='polygon')
    for ra,dec in zip(ras,decs):
        ra.append(ra[0])
        dec.append(dec[0])
        verts = np.column_stack((ra,dec))
        tot_areas = np.append(tot_areas,GenUtils.area(verts))
    return tot_areas

def area_within_annulus(reg_name):
    '''
    returns deg^s of each annulus area, that is is subtracts out the inner ones.
    deg^2 to arcsec^2: areas*3600.*3600.
    deg^2 to arcmin^2: areas*3600.
    '''
    # finds area of concentric reg shapes
    # use GenUtils.area for the area of a polygon.
    tot_areas = area_within_regs(reg_name)
        
    areas = np.array([])
    for i in range(len(tot_areas)):
        if i==0:
            areas = np.append(areas,tot_areas[i])
        else:
            areas = np.append(areas,tot_areas[i]-tot_areas[i-1])

    # deg^2 to arcsec^2: areas*3600.*3600.
    # deg^2 to arcmin^2: areas*3600.
    return np.array(areas)

def points_inside_ds9_polygon(reg_name,ra_points,dec_points):
    '''
    uses read_reg
    '''
    try:
        file_io
    except NameError:
        import file_io
    radec = np.column_stack((ra_points,dec_points))
    ras,decs = file_io.read_reg(reg_name,shape='polygon')
    verts = [np.column_stack((ras[i],decs[i])) for i in range(len(ras))]
    masks = [nxutils.points_inside_poly(radec, vert) for vert in verts]
    
    inds = [np.nonzero(masks[i]-masks[i+1])[0] for i in range(len(masks)-1)]
    return inds

def get_hphb_cmd_verts(tracks,**kwargs):
    # Critical mass
    Mass = kwargs.get('m_crit',0.58)
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    Av = kwargs.get('Av_bulge')
    af275w = kwargs.get('af275w')
    af336w = kwargs.get('af336w')
    dmod = kwargs.get('dmod')
    
    T = PadovaTracksUtils.get_track_by_mass(tracks,Mass)
    try:
        tMag1 = PadovaTracksUtils.Tracks.get_stable_col(T,Filter1+'1')
        tMag2 = PadovaTracksUtils.Tracks.get_stable_col(T,Filter2+'1')
    except:
        tMag1 = T.data_array[Filter1+'1']
        tMag2 = T.data_array[Filter2+'1']
    tmag1 = tMag1+dmod+Av*af275w
    tmag2 = tMag2+dmod+Av*af336w
    
    # verts for one track
    verts = PadovaTracksUtils.track2polygon(tmag1,tmag2)
    return verts
    
def get_HPHB_PAGB_separation(phot,rad_inds,cut_inds,**kwargs):
    '''
    rad_inds are inds of phot as a function of radius
    cut_inds are inds of phot for some quality cut.
    '''
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')

    hphb_inds = pagb_inds = nhphbs = npagbs = np.array([])
    flux1_pagbs = flux2_pagbs = flux1_hphbs = flux2_hphbs = np.array([])
    btrack = kwargs.get('btrack')
    verts = get_hphb_cmd_verts(btrack,**kwargs)
    
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[i]) & set(cut_inds))
        color = phot['mag1'][inds]-phot['mag2'][inds]
        points = np.column_stack((color,phot['mag2'][inds]))
        mask = nxutils.points_inside_poly(points, verts)
        hphb = np.nonzero(mask)[0]
        pagb = np.nonzero(abs(mask*1.-1.))[0]
        
        hphb_inds = np.append(hphb_inds,hphb)
        pagb_inds = np.append(pagb_inds,pagb)
        
        nhphbs = np.append(nhphbs,len(hphb))
        npagbs = np.append(npagbs,len(pagb))
        
        # mags
        hphb_mag1 = phot['mag1'][inds][hphb]
        hphb_mag2 = phot['mag2'][inds][hphb]
        pagb_mag1 = phot['mag1'][inds][pagb]
        pagb_mag2 = phot['mag2'][inds][pagb]
        
        # flux
        flux1_hphb = np.sum(Astronomy.vegamag2flux(
            hphb_mag1,Filter1,'wfc3uvis'))
        flux2_hphb = np.sum(Astronomy.vegamag2flux(
            hphb_mag2,Filter2,'wfc3uvis'))        
        flux1_pagb = np.sum(Astronomy.vegamag2flux(
            pagb_mag1,Filter1,'wfc3uvis'))
        flux2_pagb = np.sum(Astronomy.vegamag2flux(
            pagb_mag2,Filter2,'wfc3uvis'))

        flux1_pagbs = np.append(flux1_pagbs,flux1_pagb)
        flux2_pagbs = np.append(flux2_pagbs,flux2_pagb)
        flux1_hphbs = np.append(flux1_hphbs,flux1_hphb)
        flux2_hphbs = np.append(flux2_hphbs,flux2_hphb)
    
    return nhphbs, npagbs, flux1_hphbs, flux2_hphbs, flux1_pagbs, flux2_pagbs

def LF(phot, inds_cut1, inds_cut2, rad_inds, **kwargs):
    '''
    The cut inds should be quality cut without completeness corrections.
    
    Makes a LF and normalizes to the 275 peak (see code).
    
    There is a completeness correction on this and it only works if it uses
    the same magbins as the completeness corrections file. 
    This is hard coded into reading the completeness corrections file.
    '''
    sup_phot = kwargs.get('sup_phot')
    rif_36 = sup_phot['rif_36']
    comp_file = kwargs.get('comp_file')
    comp_data = kwargs.get('comp_data')
    frac1s, frac2s, magbins =  file_io.load_completeness_corrections(comp_file)
    comp50mag1, comp50mag2 = file_io.get_comp90mags(comp_data)
    ylab =r'$\rm{Number}/10^7 M_\odot$'
    xlab1 = '$%s$'%kwargs.get('Filter1')
    xlab2 = '$%s$'%kwargs.get('Filter2')
    cols = kwargs.get('cols')
    save_lf = kwargs.get('save_lf',False)
    if save_lf == True:
        raw_lf = kwargs.get('raw_lf',False)
        lf_file = open('lf_data.dat','w')
        print 'writing lf_data.dat'
        lf_file.write('# Filter Region magbin ; hist\n')
    ax1,ax2 = graphics.two_panel_plot(11,5,xlab1,xlab2,ylab)
    for i in range(len(rad_inds)):
        # M/L ratio 1e7 scaling is for AXIS
        norm = rif_36[i+1]/2.315e-9*0.95*1e-7
        inds1 = list(set(inds_cut1) & set(rad_inds[i]))
        hist1, bins = np.histogram(phot['mag1'][inds1],bins=magbins)
        known = np.nonzero(bins[1:] <= comp50mag1[i])[0]
        poor =  np.nonzero(bins[1:] >= comp50mag1[i])[0]
        if save_lf==True:
            lf_file.write('%s %i '%(kwargs.get('Filter1'),i))
            [lf_file.write('%f '%b) for b in bins[1:]]
            lf_file.write('; ')
            if raw_lf == True:
                print 'writing the raw LF, no norm.'
                [lf_file.write('%f '%h) for h in hist1]
            else:
                [lf_file.write('%f '%h) for h in hist1*frac1s[i]/norm]
            lf_file.write('\n')
        ax1.plot(bins[1:],hist1*frac1s[i]/norm,
                 drawstyle='steps',color=cols[i],lw=1)
        ax1.plot(bins[1:][known],hist1[known]*frac1s[i][known]/norm,
                 drawstyle='steps',color=cols[i],lw=2)
        ax1.set_xlim(26,19)
        
        inds2 = list(set(inds_cut2) & set(rad_inds[i]))
        hist2, bins = np.histogram(phot['mag2'][inds2],bins=magbins)
        known = np.nonzero(bins[1:] <= comp50mag2[i])[0]
        poor =  np.nonzero(bins[1:] >= comp50mag2[i])[0]
        if save_lf==True:
            lf_file.write('%s %i '%(kwargs.get('Filter2'),i))
            [lf_file.write('%f '%b) for b in bins[1:]]
            lf_file.write('; ')
            if raw_lf == True:
                [lf_file.write('%f '%h) for h in hist2]
            else:
                [lf_file.write('%f '%h) for h in hist2*frac2s[i]/norm]
            lf_file.write('\n')
        ax2.plot(bins[1:],hist2*frac2s[i]/norm,
                 drawstyle='steps',color=cols[i],lw=1)        
        ax2.plot(bins[1:][known],hist2[known]*frac2s[i][known]/norm,
                 drawstyle='steps',color=cols[i],lw=2)
        ax2.set_xlim(26,19)
        
    if 'save_lf'==True:
        ls_file.close()
    nomodel = kwargs.get('nomodel',False)
    if nomodel == False:
        model_F336W, n336, model_F275W, n275 = file_io.load_hphb_model(**kwargs)
        model_norm = (n336.max())/float(hist1.max())
        # hard coded model normalization below
        ax2.plot(model_F336W,n336/23000.,':',lw=2,label='HP-HB',color='black')
        ax1.plot(model_F275W,n275/23000.,':',lw=2,label='HP-HB',color='black')
    ax1.set_ylim(0,6)
    ax2.set_ylim(0,6)
    plt.savefig('radialLF%s'%EXT)
    plt.close('all')
    return

def integrated_light_data(nhphbs,npagbs,flux1_hphbs,flux2_hphbs,flux1_pagbs,flux2_pagbs,**kwargs):
    # integrated flux counts to cgs units:
    int_flux = kwargs.get('int_flux')

    int_counts1 = int_flux['f275w_cts']
    int_counts2 = int_flux['f336w_cts']
    
    vega1counts = 1.1413e9
    vega1flx = 3.727e-9
    
    vega2counts = 2.4318e9
    vega2flx = 3.245e-9
    
    intflx1_cgs = int_counts1*vega1flx/vega1counts
    intflx2_cgs = int_counts2*vega2flx/vega2counts
    
    print 'number of pagbs in inner two bins compared to total stars in inner two bins.'
    frac_num_pagbs_inner2 = np.sum(npagbs[0:2])/(np.sum(npagbs[0:2])+
                                                 np.sum(nhphbs[0:2]))
    print frac_num_pagbs_inner2
    
    print 'resolved flux from hphb stars to total resolved flux in the inner two bins'
    tot_flux1 = flux1_hphbs+flux1_pagbs
    tot_flux2 = flux2_hphbs+flux2_pagbs
    frac_flux1_hphb_inner2 = np.sum(flux1_hphbs[0:2])/np.sum(tot_flux1[0:2])
    frac_flux2_hphb_inner2 = np.sum(flux2_hphbs[0:2])/np.sum(tot_flux2[0:2])
    print 'f275:',frac_flux1_hphb_inner2
    print 'f336:',frac_flux2_hphb_inner2
    
    print 'fraction of resolved flux to total integrated flux'
    resolved_to_integrated1 = np.sum(tot_flux1)/np.sum(intflx1_cgs)
    resolved_to_integrated2 = np.sum(tot_flux2)/np.sum(intflx2_cgs)
    print 'f275:',resolved_to_integrated1
    print 'f335:',resolved_to_integrated2
    
    return
    
def hphb_pagb_comp_plot(nhphbs, npagbs, flux1_hphbs, flux2_hphbs,
                        flux1_pagbs, flux2_pagbs,**kwargs):
    dmod = kwargs.get('dmod')
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    ftsize = kwargs.get('ftsize',20)
    # loading kwargs without innermost region.
    sup_phot = kwargs.get('sup_phot')[1:]
    int_flux = kwargs.get('int_flux')
    areas = kwargs.get('areas')[1:]
    
    # distance from center
    contour_arcsec = sup_phot['semimajor_axis']
    contour_pc = Astronomy.arcsec2kpc(contour_arcsec,dmod=dmod)*1000.
    
    tot_flux1 = flux1_hphbs + flux1_pagbs
    tot_flux2 = flux2_hphbs + flux2_pagbs


    # surface brightness
    sb1_hphbs = Astronomy.flux2vegamag(flux1_hphbs/areas,
                                       Filter1,'wfc3uvis')
    sb2_hphbs = Astronomy.flux2vegamag(flux2_hphbs/areas,
                                       Filter2,'wfc3uvis')
    sb1_pagbs = Astronomy.flux2vegamag(flux1_pagbs/areas,
                                       Filter1,'wfc3uvis')
    sb2_pagbs = Astronomy.flux2vegamag(flux2_pagbs/areas,
                                       Filter2,'wfc3uvis')
    
    # integrated surface brightness
    int_sb1 = int_flux['f275w_sb']
    int_sb2 = int_flux['f336w_sb']
    
    # Sersic plot see function for where numbers came from
    R_range,mu = sersic_profile(21.5,2.2,1.)

    # hphb/pagb ratio and errors
    h2p_ratio2 = flux2_hphbs/flux2_pagbs
    h2p_ratio1 = flux1_hphbs/flux1_pagbs
    
    # error estimate
    ratio_err = np.sqrt((np.sqrt(nhphbs)/nhphbs)**2+
                        (np.sqrt(npagbs)/npagbs)**2)
    
    f1_err = h2p_ratio1*ratio_err
    f2_err = h2p_ratio2*ratio_err
    
    f2_uerror = h2p_ratio2+f2_err
    f2_lerror = h2p_ratio2-f2_err
    f1_uerror = h2p_ratio1+f1_err
    f1_lerror = h2p_ratio1-f1_err

    # plots
    f1_color = 'black'
    f2_color = 'red'
    hphb_sym = ':'
    pagb_sym = '--'

    ax2,ax1 = graphics.two_panel_plot_v2(12,6,'x','x','y')

    ax1.set_ylabel(
        r'$\rm{Resolved\ Flux\ Ratio\ } \frac{\rm{HP-HB}}{\rm{P-AGB}}$',
        fontsize=15)

    graphics.error_wedge(ax1,f2_color,contour_pc,f2_uerror,f2_lerror,
                         alpha=0.2,lw=0)
    #graphics.error_wedge(ax1,f1_color,contour_pc,f1_uerror,f1_lerror,
    #                     alpha=0.2,lw=0)
    
    ax1.plot(contour_pc,h2p_ratio2,color=f2_color,lw=2)
    #ax1.plot(contour_pc,h2p_ratio1,color=f1_color,lw=2)
    
    # Surface brightness plot 

    ax2.set_ylabel(
        r'$\rm{Surface\ Brightness\ (mag\ arcsec^{-2})}$',
        fontsize=ftsize-3)
    ax2.plot(contour_pc,sb2_hphbs,hphb_sym,color=f2_color,lw=2)
    ax2.plot(contour_pc,sb2_pagbs,pagb_sym,color=f2_color,lw=2)
    ax2.plot(contour_pc,sb1_hphbs,hphb_sym,color=f1_color,lw=2)
    ax2.plot(contour_pc,sb1_pagbs,pagb_sym,color=f1_color,lw=2)
    
    ax2.plot(contour_pc,int_sb1,color=f1_color,lw=4)
    ax2.plot(contour_pc,int_sb2,color=f2_color,lw=4)
    # Legends
    ax2.plot(1e6,1e6,color='black',lw=4,label='Integrated')
    ax2.plot(1e6,1e6,hphb_sym,color='black',label='HP-HB',lw=2)
    ax2.plot(1e6,1e6,pagb_sym,color='black',label='P-AGB',lw=2)
    ax2.plot(R_range*1000.,mu,color='gray',label='n=2.2; Courteau et al. 2011')
    ax2.set_ylim(28,17)
    
    for ax in [ax1,ax2]:
        ax.set_xlabel(r'$\rm{Radius\ (pc)}$',fontsize=ftsize)
        ax.set_xlim(40,1000)
        ax.set_xscale('log')
    
    # align y axis labels
    ax1.yaxis.set_label_coords(-0.09,0.5)
    ax2.yaxis.set_label_coords(-0.07,0.5)

    #ax1.annotate(r'$F275W$',xy=(50,.55),color=f1_color,
    #             ha='left',fontsize=ftsize,**graphics.load_ann_kwargs())
    #ax1.annotate(r'$F336W$',xy=(50,.45),color=f2_color,
    #             ha='left',fontsize=ftsize,**graphics.load_ann_kwargs())

    ax2.legend(frameon=False,loc=3)
    graphics.put_feh_on_top(contour_pc,ax1,twoplot=True,offset_text =0.01,
                   offset_tick = -0.02,shk = 0.2)
    graphics.put_feh_on_top(contour_pc,ax2,twoplot=True,offset_text = -.1,
                   offset_tick = 0.2,shk = 0.0)
    plt.savefig('res_fluxes%s'%EXT)
    plt.close('all')
    print 'wrote res_fluxes%s'%EXT
    return

def mp_sersic(p, fjac = None, x = None, y = None, err = None):
    '''
    sersic profile for mpfit
    '''
    model = sersic(x, p)
    status = 0
    return [status, (y - model) / err]

def get_sersic_limits(**kwargs):
    from scipy import optimize
    from mpfit import mpfit
    int_flux = kwargs.get('int_flux')
    int_sb1 = int_flux['f275w_sb']
    int_sb2 = int_flux['f336w_sb']
    sup_phot = kwargs.get('sup_phot')[1:]
    dmod = kwargs.get('dmod')
    # distance from center
    contour_arcsec = sup_phot['semimajor_axis']
    contour_pc = Astronomy.arcsec2kpc(contour_arcsec,dmod=dmod)*1000.
    '''
    # uniform errors
    err = np.zeros(len(int_sb1))+.0000001
    
    # set up inputs
    #R = np.arange(0,12000,0.1)
    R = contour_pc
    p_in = {'x':R, 'y':int_sb1, 'err':err}
    
    # set up initial parameters:
    # norm = max(hist),
    # mean set to be half mean, and 3/2 mean,
    # sigma set to be same as dcol spacing...
    p0 = [25,2.2,1000.]
    parinfo = [{'limited': [False,False], 'parname':'mu_e'},
               {'limited': [False,True], 'limits': (1.,3.), 'parname':'n'},
               {'limited': [False,True], 'limits': (0,3000.), 'parname':'R_e'}]
    
    mp_dg = mpfit(mp_sersic, p0,functkw = p_in, quiet = 1, parinfo=parinfo) 
    '''
    sersic = lambda p,x: p[0]+\
             2.5*(1.992*p[1]-0.3271)/np.log(10.)*((x/p[2])**(1./p[1])-1.)
    errfunc = lambda p, x, y: sersic(p,x)-y
    
    p0 = [21,2.2,1000.]
    p1, success  = optimize.leastsq(errfunc, p0[:], args=(contour_pc, int_sb1),full_output=0)
    print 're-do if more than 4:',success
    p0 = [21,2.2,1000.]
    p2, success = optimize.leastsq(errfunc, p0[:], args=(contour_pc, int_sb2),full_output=0)
    print 'get_sersic_limits: re-do if more than 4:',success
    print '275 sersic: mu_e,n,Re (pc):',p1
    print '336 sersic: mu_e,n,Re (pc):',p2
    # in case you wanna plot... 
    radius = np.arange(40,13000)
    corrfit1 = sersic(p1,radius)
    corrfit2 = sersic(p2,radius)
    corrfit1 = sersic(mp_dg.params,radius)
    plt.plot(radius,corrfit2)   
    plt.plot(radius,corrfit1) 
    corrfit1 = sersic(p1,contour_pc)
    corrfit2 = sersic(p2,contour_pc)
    plt.plot(contour_pc,corrfit2,'black')   
    plt.plot(contour_pc,corrfit1,'black')    
    return

def sersic(x,p):
    '''
    same as sersic proile, just with better passing variables.
    p = [mu_e, n, R_e]
    '''
    #if (p[1] > 0.5) and (p[1] < 10.):
    bn = 1.992*p[1]-0.3271
    #else:
    #    print 'sersic profile: please choose 0.5<n<10 or do some damn coding'
        
    y = p[0]+2.5*bn/np.log(10.)*((x/p[2])**(1/p[1])-1)
    return y
    
def sersic_profile(mu_e,n,Re,R_range=None):
    '''
    Get limits using this:
    # Fit sersic profile
    sersic = lambda p,x: p[0]+\
             2.5*(1.992*p[1]-0.3271)/np.log(10.)*((x/p[2])**(1./p[1])-1.)
    errfunc = lambda p, x, y: sersic(p,x)-y
    
    p0 = [20,2.2,1.]
    p1, success = optimize.leastsq(errfunc, p0[:], args=(rpc[::-1], intsb1))
    print 're-do if more than 4:',success
    p0 = [25,2.2,1.]
    p2, success = optimize.leastsq(errfunc, p0[:], args=(rpc[::-1], intsb2))
    print 're-do if more than 4:',success
    corrfit1 = sersic(p1,rpc[::-1])
    corrfit2 = sersic(p2,rpc[::-1])
    print '275 sersic: mu_e,n,Re',p1
    print '336 sersic: mu_e,n,Re',p2
    '''
    if (n > 0.5) and (n < 10.):
        bn = 1.992*n-0.3271
    else:
        print 'sersic profile: please choose 0.5<n<10 or do some damn coding'
        
    if R_range == None: R_range = np.linspace(0,10,1000) # assuming kpc
    
    mu = mu_e+2.5*bn/np.log(10.)*((R_range/Re)**(1/n)-1)
    return R_range,mu

def frac_flux_plot(phot,rad_inds,cut_inds,**kwargs):
    areas = kwargs.get('areas')
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    cols = kwargs.get('cols')

    ylab = 'Flux ($10^{-3}$ erg s$^{-1}$ cm$^{-2}$ 0.1 mag$^{-1}$)'
    xlab1 = '$F275W$'
    xlab2 = '$F336W$'
    ax1,ax2 = graphics.two_panel_plot(10,5,xlab1,xlab2,ylab)
    mstep = 0.4
    ax1.set_ylabel('Flux Density ($10^{-19}$erg s$^{-1}$ cm$^{-2}$ %.1fmag$^{-1}$ $\prime\prime^{-2}$)'%mstep,fontsize=15)
    magbins = np.arange(19,26,mstep)
    Tot_flux1,Tot_flux2 = [],[]
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[i]) & set(cut_inds))
        m1 = phot['mag1'][inds]
        m2 = phot['mag2'][inds]
        
        flux1 = Astronomy.vegamag2flux(m1,Filter1,'wfc3uvis')
        flux2 = Astronomy.vegamag2flux(m2,Filter2,'wfc3uvis')
        
        binds1 = np.digitize(m1,magbins)-1
        binds2 = np.digitize(m2,magbins)-1
        
        added_flux1 = []
        added_flux2 = []
        for j in range(len(magbins)):
            added_flux1.append(np.sum(flux1[np.nonzero(binds1==j)[0]]))
            added_flux2.append(np.sum(flux2[np.nonzero(binds2==j)[0]]))

        Tot_flux1.append(np.sum(added_flux1))
        Tot_flux2.append(np.sum(added_flux2))
        
        comp_data = kwargs.get('comp_data')
        comp50mag1, comp50mag2 = file_io.get_comp50mags(filename=comp_data)
        known1 = np.nonzero(magbins <= comp50mag1[i])[0]
        poor1 =  np.nonzero(magbins >= comp50mag1[i])[0]
        known2 = np.nonzero(magbins <= comp50mag2[i])[0]
        poor2 =  np.nonzero(magbins >= comp50mag2[i])[0]

        norm = 1e19/areas[i]
        ax1.plot(magbins[known1],np.array(added_flux1)[known1]*norm,
                 ls='steps',color=cols[i],lw=2)
        ax2.plot(magbins[known2],np.array(added_flux2)[known2]*norm,
                 ls='steps',color=cols[i],lw=2)
        ax1.plot(magbins[poor1],np.array(added_flux1)[poor1]*norm,
                 ls='steps--',color=cols[i],lw=1.5)
        ax2.plot(magbins[poor2],np.array(added_flux2)[poor2]*norm,
                 ls='steps--',color=cols[i],lw=1.5)
        ax1.set_xlim((25,19))
        ax2.set_xlim((26,19))

    plt.savefig('integrated_LF_area%s'%EXT)
    plt.close('all')
    print 'wrote integrated_LF_area%s'%EXT            
    return Tot_flux1,Tot_flux2

def brown_fct(t_act_hphb,**kwargs):
    '''
    This was made to check the fct around brown98 area.
    '''
    sup_phot = kwargs['sup_phot']
    btrack = kwargs.get('btrack')
    mag1_cut = kwargs.get('mag1_cut')
    mag2_cut = kwargs.get('mag2_cut')
    frac_bulge = sup_phot['bulge_fraction'][1] # Bulge fraction B/(B+D)
    
    # file locations
    file_src = '/astro/net/angst2/philrose/PHAT/paper/supporting_files/'
    if not os.path.isdir(file_src):
        file_src = '/Users/phil/research/PHAT//paper/supporting_files/'
    b98region = os.path.join(file_src,'brown98_region_masked.reg')
    b98fits = os.path.join(file_src,'phat_photcat_b98region.fits')
    
    # load b98 photometry file
    b98_phat = pyfits.open(b98fits)[1].data
    b98_ra = b98_phat.field('ra')
    b98_dec = b98_phat.field('dec')
    b98_mag1 = b98_phat.field('mag1_uvis')
    b98_mag2 = b98_phat.field('mag2_uvis')
    b98_color = b98_mag1-b98_mag2
    
    # load b98 region
    ras,decs = file_io.read_reg(b98region,shape='polygon')
    foc_fov = np.column_stack((ras[0],decs[0]))
    
    # b98 called uv-brights stars 1 mag brighter than completeness limit
    b98_mag1_cut = mag1_cut - 1.
    b98_mag2_cut = mag2_cut - 1.
    inds1 = np.nonzero(b98_mag1 < b98_mag1_cut)[0]
    inds2 = np.nonzero(b98_mag2 < b98_mag2_cut)[0]
    b98_cut = list(set(inds1) & set(inds2))
    
    # load the hphb color-mag limits
    verts = get_hphb_cmd_verts(btrack,**kwargs)
    
    # which b98 stars are in the hphb color-mag limit
    points = np.column_stack((b98_color[b98_cut],b98_mag2[b98_cut]))
    #points = np.column_stack((b98_color,b98_mag2))
    inds = np.nonzero(nxutils.points_inside_poly(points, verts))[0]
    
    nhphb_foc = float(len(inds))
    
    rif_36_foc = 0.527
    L_T = rif_36_foc*1.61e8*frac_bulge #230893.52211437
    Bt = 2.2e-11 # From Spectral Evolution of Galaxies ASSL 122 p195+
    t_hphb = nhphb_foc/(Bt*L_T)
        
    # hard coded expected lifetimes.
    # median time_spent_on_cmd for the three largest mass tracks
    #t_act_hphb =  2e6#1e7 # 10694600 # 994000.0
    print 'Brown Check:'
    print 'NHPHB expected:',L_T*Bt*t_act_hphb
    print 'percent msto that become x as a function of radius:'
    print 'HPHB',t_hphb/t_act_hphb*100
    print 'nhbhbs',len(inds)
    return

def bc_calc():
    
    # this was the old way to do it. Trying to get the BC ended up 
    # underpredicting L_T of the bulge.
    # calculation to get total luminosity see lab book 2.
    ssp_files = ['ssp_1e10_zrange.dat',
                 'ssp_1e9_zrange.dat',
                 'ssp_solar_old.dat']

    names = open(ssp_file,'r').readlines()[9].replace('#','').split()
    ssp = np.genfromtxt(ssp_file,names=names)
    age = ssp['Age']
    BC = ssp['mbol']-ssp['IRAC_36']
    z = ssp['Z']

    plt.plot(BC[np.argsort(age)],np.sort(age),'o')
    plt.ylabel('Age (log yr)')
    plt.yscale('log')
    plt.axis([2.49,2.7,9e8,2e10])
    
    plt.plot(np.array(BC)[np.argsort(z)],np.sort(z))
    plt.plot(BC,z,'o',color='blue')
    plt.ylabel('z')
    
    plt.xlabel('BC')
    plt.title(ssp_file.replace('_',' '))
    plt.savefig('BC'+ssp_file.replace('dat','png'))

    
    BC = 2.60
    Mbol_sun = 4.77
    vega_zero = 280.9
    f_f36 = 10**(0.4*(BC-Mbol_sun))

    Jy_at_10pc = 1/(f_f36*vega_zero) # L_sun
    Jy_at_m31 = 10**(-0.4*(-2.5*np.log10(Jy_at_10pc)-dmod))
    return

def saglia_zh():
    '''
    reads data from saglia et al 2010, returns radius in arcsec and ave z/h. 
    rewrite to return only sag_data if you want the entire data set.
    '''
    sagf = '/Users/phil/research/PHAT/UVBrightpaper/supporting_files/aveSSP_Dalcanton.dat'
    header = open(sagf).readlines()[1].replace('#','').strip().split()
    sag_data = np.genfromtxt(sagf,names=header)
    return sag_data
 
def anti_cor_plot(**kwargs):
    ftsize = kwargs.get('ftsize')
    dmod = kwargs.get('dmod')
    sup_phot = kwargs['sup_phot']
    # excluding inner regions
    fuv_nuv = sup_phot['fuv_nuv_color']
    
    # distance from center
    contour_arcsec = sup_phot['semimajor_axis']
    contour_pc = Astronomy.arcsec2kpc(contour_arcsec,dmod=dmod)*1000.
    sag_data = saglia_zh()
    sag_r = sag_data['aveR']
    sag_zh = sag_data['aveZH']
    sag_pzh = sag_data['paveZH']
    sag_mzh = sag_data['maveZH']
    sag_pc = Astronomy.arcsec2kpc(sag_r,dmod=dmod)*1000.
    fig = plt.figure(1, figsize=(10,8))
    ax = plt.axes()
    lw =4 
    ax.plot(contour_pc,fuv_nuv,color='blue',lw=lw)
    ax.set_ylabel(r'$FUV-NUV$',fontsize=ftsize)
    ax.set_xlabel(r'$\rm{Radius\ (pc)}$',fontsize=ftsize)
    #ax.set_xscale('log')
    ax1 = ax.twinx()
    ax1.plot(sag_pc,sag_zh,'--',color='black',lw=lw)
    graphics.error_wedge(ax1,'black',sag_pc,sag_pzh,sag_mzh,lw=0,alpha=0.1)
    ax1.set_ylim(-0.16,0.33)
    ax1.set_ylabel(r'$\rm{[Z/H]}$',fontsize=ftsize)
    #ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.set_xlim(10,850)
    ax.set_ylim(.76,1.31)
    
    
def time_spent(npagbs,nhphbs,**kwargs):
    ftsize = kwargs.get('ftsize')
    dmod = kwargs.get('dmod')
    sup_phot = kwargs['sup_phot']
    # excluding inner regions
    rif_36 = sup_phot['rif_36'][1:]
    noerror = kwargs.get('noerror',False)
    frac_bulge = sup_phot['bulge_fraction'][1:] # Bulge fraction B/(B+D)
    
    # distance from center
    contour_arcsec = sup_phot['semimajor_axis'][1:]
    contour_pc = Astronomy.arcsec2kpc(contour_arcsec,dmod=dmod)*1000.
    
    # From Spectral Evolution of Galaxies ASSL 122 p195+
    Bt = 2.2e-11
    # From cliff's mass/light ratio calc using Oh 98.
    L_T = rif_36/2.315e-9
    L_T = rif_36*1.61e8*frac_bulge
    t_pagb = npagbs/(Bt*L_T)
    t_hphb = nhphbs/(Bt*L_T)
    
    # hard coded expected lifetimes.
    # mean time_spent_on_cmd for the three largest mass tracks
    t_act_hphb =  2e6#1740425.2#8685282.75# 36466400.0 #816000.0#1e7 # 10694600 ## 816000.0 36466400.0
    t_act_pagb = 1e4
    print 'NPAGB expected:',L_T*Bt*t_act_pagb
    print 'NHPHB expected:',L_T*Bt*t_act_hphb
    print 'NPAGB',npagbs
    print 'NHPHB',nhphbs
    print 'percent msto that become x as a function of radius:'
    print 'HPHB',t_hphb/t_act_hphb*100
    print 'PAGB', t_pagb/t_act_pagb*100
    print 'increase HPHB:', t_hphb[0]/t_act_hphb/(t_hphb[-1]/t_act_hphb)

    
    # error for plot
    err = t_hphb/t_act_hphb*100*1./np.sqrt(nhphbs)
    p_err = t_hphb/t_act_hphb*100+err
    m_err = t_hphb/t_act_hphb*100-err
    
    # plots
    fig = plt.figure(1, figsize=(10,8))
    ax = plt.axes()
    lw =4 
    if noerror == False:
        graphics.error_wedge(ax,'blue',contour_pc,p_err,m_err,lw=0,alpha=0.2)
        lw = 2
    ax.plot(contour_pc,t_hphb/t_act_hphb*100,color='blue',lw=lw)
    ax.set_ylabel(r'$\rm{\%\ of\ MSTO\ stars\ that\ pass\ though\ HP-HB}$',
                  fontsize=ftsize)
    ax.set_xlabel(r'$\rm{Radius\ (pc)}$',
                  fontsize=ftsize)
    ax.set_xscale('log')
    #ax.set_ylim(0.25,1.75)
    #graphics.put_feh_on_top(contour_pc,ax,offset_text=0.01,
    #                        offset_tick=-0.05,shk=0.3)
    
    
    # Add saglia
    ax1 = ax.twinx()

    sag_data = saglia_zh()
    sag_r = sag_data['aveR']
    sag_zh = sag_data['aveZH']
    sag_pzh = sag_data['paveZH']
    sag_mzh = sag_data['maveZH']
    
    sag_pc = Astronomy.arcsec2kpc(sag_r,dmod=dmod)*1000.
    ax1.plot(sag_pc,sag_zh,'--',color='black',lw=lw)
    #ax1.semilogy(sag_pc,sag_pzh,'--',color='black',lw=lw)
    #ax1.semilogy(sag_pc,sag_mzh,'--',color='black',lw=lw)
    graphics.error_wedge(ax1,'black',sag_pc,sag_pzh,sag_mzh,lw=0,alpha=0.1)
    ax1.set_ylim(-0.16,0.2)
    ax1.set_ylabel(r'$\rm{[Z/H]}$',fontsize=ftsize)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.set_xlim(50,820)
    plt.savefig('FCT%s'%EXT)
    plt.close('all')
    print 'wrote FCT%s'%EXT

def numbers_of_pops(npagbs, nhphbs,**kwargs):
    ftsize = kwargs.get('ftsize')
    dmod = kwargs.get('dmod')
    sup_phot = kwargs['sup_phot']
    colors = {'hphb': 'red',
              'pagb': 'purple',
              'lmxb': 'darkgreen',
              '3.6um': 'black'}
    # excluding inner regions
    rif_36 = sup_phot['rif_36'][1:]
    frac_bulge = sup_phot['bulge_fraction'][1:] # Bulge fraction B/(B+D)
    areas = kwargs.get('areas')[1:]
    
    # distance from center
    contour_arcsec = sup_phot['semimajor_axis'][1:]
    contour_pc = Astronomy.arcsec2kpc(contour_arcsec,dmod=dmod)*1000.
    
    # calculation to get total luminosity
    dmod2pc = 10**(0.2*dmod+1.)*3.085e18 # distance to M31 in cm 
    fourpiRsq = 4*np.pi*dmod2pc**2
    jy2ergscm2 = 1e-23*3e18*0.75e4/(3.550e4)**2
    Lsun = 3.85e33
    factor = jy2ergscm2*fourpiRsq/Lsun*280.9
    L = rif_36*factor*frac_bulge
    
    # can't remember which file I stole these from...
    Npne = [1,5,19,38,93]
    lmxbs = np.array([5.,9.,14.,17.,20.])

    # error calcs
    hphb_uerr,hphb_lerr = poisson_error(nhphbs,norm=areas)
    pagb_uerr,pagb_lerr = poisson_error(npagbs,norm=areas)
    lmxb_uerr,lmxb_lerr = poisson_error(lmxbs,norm=areas)
    
    fig = plt.figure(1, figsize=(10,10))
    ax = plt.axes()
    ax1 = ax.twinx()
    graphics.error_wedge(ax,colors['hphb'],contour_pc,hphb_uerr,hphb_lerr,
                         lw=0,alpha=0.2)
    ax.plot(contour_pc,nhphbs/areas,color=colors['hphb'],
            lw=2,label=r'HP-HB (%i)'%int(np.sum(nhphbs)))
    
    graphics.error_wedge(ax,colors['pagb'],contour_pc,pagb_uerr,pagb_lerr,
                         lw=0,alpha=0.2)
    
    ax.plot(contour_pc,npagbs/areas,color=colors['pagb'],
            lw=2,label=r'P-AGB (%i)'%int(np.sum(npagbs)))
    
    graphics.error_wedge(ax,colors['lmxb'],contour_pc,lmxb_uerr,lmxb_lerr,
                         lw=0,alpha=0.2)
    
    ax.plot(contour_pc,lmxbs/areas,color=colors['lmxb'],
            lw=2,label=r'LMXB (%i)'%int(np.sum(lmxbs)))
    ax.set_ylabel(r'$\rm{Resolved\ Stellar\ Population\ Surface\ Density\ (arcsec}^{-2}\rm{)}$',fontsize=ftsize)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(40,1000)
    ax.set_ylim(1e-4,2)
    ax1.semilogy(contour_pc,L*frac_bulge/areas,'--',color=colors['3.6um'],lw=3,
                 label='$L_{3.6\mu m}/$area $(L_{\odot}/\prime\prime^2)$')
    
    # trick the legend
    ax.plot(1e6,1e6,'--',lw=3,color=colors['3.6um'],
            label='$L_{3.6\mu m}/$area $(L_{\odot}/\prime\prime^2)$')
    ax1.set_yscale(ax.get_yscale())
    ax1.set_xscale(ax.get_xscale())
    ax1.set_xlim(ax.get_xlim())
    ax1.set_ylabel(r'$L_{3.6\mu m}/\rm{area}\ (L_{\odot}/\rm{arcsec}^2)$',
                   fontsize=ftsize)
    ax1.set_ylim(1.2e5,2.5e9)
    ax.set_xlim(40,1000.)
    ax.set_xlabel(r'$\rm{Radius\ (pc)}$',fontsize=ftsize)
    ax.legend(loc=0,frameon=False)
    
    # annotations
    pli_hphb,res_hphb = power_law_slopes_hand(contour_pc,nhphbs/areas)
    pli_pagb,res_pagb = power_law_slopes_hand(contour_pc,npagbs/areas)
    pli_lmxb,res_lmxb = power_law_slopes_hand(contour_pc,lmxbs/areas)
    pli_l36,res_l36 = power_law_slopes_hand(contour_pc,L/areas)
    
    numpmerr = r'$n=%.2f\pm%.2f$'
    ax.annotate(numpmerr%(pli_hphb,res_hphb),xy=(70,0.35),
                color=colors['hphb'],ha='left',fontsize=ftsize-3,
                **graphics.load_ann_kwargs())
    ax.annotate(numpmerr%(pli_pagb,res_pagb),xy=(70,0.09),
                color=colors['pagb'],ha='left',fontsize=ftsize-3,
                **graphics.load_ann_kwargs())
    ax.annotate(numpmerr%(pli_lmxb,res_lmxb),xy=(70,0.0095),
                color=colors['lmxb'],ha='left',fontsize=ftsize-3,
                **graphics.load_ann_kwargs())
    ax.annotate(numpmerr%(pli_l36,res_l36),xy=(70,0.0013),
                color=colors['3.6um'],ha='left',fontsize=ftsize-3,
                **graphics.load_ann_kwargs())
    
    graphics.put_feh_on_top(contour_pc,ax)
    
    plt.savefig('res_stellar_pops%s'%EXT)
    plt.close('all')
    print 'wrote res_stellar_pops%s'%EXT
    return

def poisson_error(arr,norm=1.):
    p = arr/norm+np.sqrt(arr)/norm
    m = arr/norm-np.sqrt(arr)/norm
    return p,m

def power_law_slopes_hand(x,y):
    ''' A+B*x '''
    logx = np.log10(x)
    logy = np.log10(y)
    N = float(len(logx))
    Anum = np.sum(logx*logx)*np.sum(logy)-np.sum(logx)*np.sum(logx*logy)
    Bnum = N*np.sum(logx*logy)-np.sum(logx)*np.sum(logy)
    D = N*np.sum(logx*logx)-np.sum(logx)*np.sum(logx) 
    
    A = Anum/D
    B = Bnum/D
    
    sig = np.sqrt(1/(N-2)*np.sum((logy-A-B*logx)**2))
    return B,sig

def cmds_and_tracks(phot,rad_inds,cut_inds,**kwargs):
    from matplotlib.patches import FancyArrow
    dmod = kwargs.get('dmod')
    ftsize = kwargs.get('ftsize')
    disk_fit = kwargs.get('disk_fit')
    mag1_cut = kwargs.get('mag1_cut')
    mag2_cut = kwargs.get('mag2_cut')

    btrack = kwargs.get('btrack')
    btrackmore = kwargs.get('btrackmore')
    dtrack = kwargs.get('dtrack')
    bmstrack = kwargs.get('bmstrack')

    Av_bulge = kwargs.get('Av_bulge')
    Av_disk = kwargs.get('Av_disk')
    af275w = kwargs.get('af275w')
    af336w = kwargs.get('af336w')
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')

    print disk_fit
    d1 = pyfits.open(disk_fit)[1].data
    dmag1 = d1.field('mag1')
    dmag2 = d1.field('mag2')
    dqc_inds = np.nonzero( (dmag1 < mag1_cut) & (dmag2 < mag2_cut))
    dmag = dmag2[dqc_inds]
    dcolor = dmag1[dqc_inds]-dmag
    
    rad_inds_inner2 = np.concatenate((rad_inds[0:2]))
    inds_inner2 = list(set(rad_inds_inner2) & set(cut_inds))
    m2 = phot['mag2'][inds_inner2]
    bcolor = phot['mag1'][inds_inner2]-phot['mag2'][inds_inner2]
    m1 = phot['mag1'][inds_inner2]
    
    ax_ll,ax_lr,ax_ul,ax_ur = graphics.setup_four_panel() 
    for ax in [ax_ll,ax_ul]:
        ax.plot(m1-m2,m2,'o',ms=1.5,mec=None,color='black')
    
    for ax in [ax_lr,ax_ur]:
        ax.plot(dcolor,dmag,'o',ms=1.5,mec=None,color='black')

    print 'Disk Tracks'
    graphics.plot_tracks_on_four_panel_smooth(dtrack,ax_lr,ax_ur,Av_disk,
                                              af275w,af336w)
    print 'HB Bulge Tracks'
    graphics.plot_tracks_on_four_panel_smooth(btrack,ax_ll,ax_ul,Av_bulge,
                                              af275w,af336w)
    #graphics.plot_tracks_on_four_panel(btrack,ax_ll,ax_ul,Av_bulge,**kwargs)

    # for Bulge MS tracks:
    print 'Bulge MS tracks'
    for track in bmstrack:
        tm1,tm2,age,HeB,Hfus = PadovaTracksUtils.get_Hfus_HeB(track,
                                                              Filter1,
                                                              Filter2)
        #Mag 2 mag:
        tm1= tm1+dmod+Av_bulge*af275w
        tm2= tm2+dmod+Av_bulge*af336w
        t_col = tm1-tm2
        #tmag = tm1
        tmag = tm2
        mass = PadovaTracksUtils.info_from_track_filename(track.name)['M']
        if int(mass)/mass == 1.:
            try:
                ax_ul.plot(t_col[Hfus],tmag[Hfus],
                           color='darkblue',zorder=3,lw=2)
                ax_ul.plot(t_col[Hfus],tmag[Hfus],
                           color='blue',zorder=3,lw=1.5,label=r'MS')
                tmp_col,tmp_mag = graphics.age_ticks(age[Hfus],
                                                     ax_ul,
                                                     t_col[Hfus],
                                                     tmag[Hfus])
                if len(Hfus) != 0: 
                    ax_ul.annotate('%.2f'%mass,
                                   xy=(min(t_col[Hfus]),
                                       tmag[np.argmin(t_col[Hfus])]),
                                   color='blue',ha='right',
                                   **graphics.load_ann_kwargs())
            except ValueError:
                pass

    ax_ll.annotate(r'%s'%('0.48'),xy=(-0.36,23.6),
                   color='red',ha='left',
                   **graphics.load_ann_kwargs())
    ax_ul.annotate(r'%s'%('MS'),xy=(0.6,19),
                   color='blue',ha='left',fontsize=15,
                   **graphics.load_ann_kwargs())
    ax_ll.annotate(r'%s'%('HP-HB'),xy=(0.6,22),
                   color='red',ha='left',fontsize=15,
                   **graphics.load_ann_kwargs())
    ax_ll.annotate(r'%s'%('P-AGB'),xy=(0.6,19),
                   color='purple',ha='left',fontsize=15,
                   **graphics.load_ann_kwargs())

    Rv31 = Astronomy.reddening_vector(Filter1,Filter2,'wfc3uvis',
                                     3.1,Teff=15000)

    Rv5 = Astronomy.reddening_vector(Filter1,Filter2,'wfc3uvis',
                                    5,Teff=15000)
    arr5 = Astronomy.reddening_vector_patch(Rv5,
                                           mag2_lim=(19,20),
                                           color_start=0.4,
                                           color='darkred',
                                           linewidth=1,
                                           length_includes_head=False,
                                           head_width=.05)
    
    arr31 = Astronomy.reddening_vector_patch(Rv31,
                                            mag2_lim=(19,20),
                                            color_start=0.6,
                                            color='darkred',
                                            linewidth=1,
                                            length_includes_head=False,
                                            head_width=.05)
    ax_ur.add_patch(arr31)
    ax_ur.add_patch(arr5)
    ax_ur.annotate('%s'%(r'$3.1$'),xy=(0.6,19),
                   color='darkred',ha='center',fontsize=15,
                   **graphics.load_ann_kwargs())
    ax_ur.annotate('%s'%(r'$5$'),xy=(0.4,19),
                   color='darkred',ha='center',fontsize=15,
                   **graphics.load_ann_kwargs())
    ax_ur.annotate('%s'%(r'$R_{V}$'),xy=(0.2,19),
                   color='darkred',ha='center',fontsize=15,
                   **graphics.load_ann_kwargs())
    
    for ax in [ax_ul,ax_ll,ax_ur,ax_lr]:
	    majorLocator = MultipleLocator(0.4)
	    ax.xaxis.set_major_locator(majorLocator)
    plt.savefig('comp_cmds%s'%EXT)
    plt.close('all')

    print 'wrote comp_cmds%s'%EXT
    return

def one_cmd_and_track(phot,rad_inds,cut_inds,**kwargs):
    dmod = kwargs.get('dmod')
    ftsize = kwargs.get('ftsize')
    mag1_cut = kwargs.get('mag1_cut')
    mag2_cut = kwargs.get('mag2_cut')

    btrack = kwargs.get('btrack')
    
    Av_bulge = kwargs.get('Av_bulge')
    af275w = kwargs.get('af275w')
    af336w = kwargs.get('af336w')
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    yfilter = kwargs.get('yfilter',Filter2)
    
    rad_inds_inner2 = np.concatenate((rad_inds[0:2]))
    inds_inner2 = list(set(rad_inds_inner2) & set(cut_inds))
    m2 = phot['mag2'][cut_inds]
    bcolor = phot['mag1'][inds_inner2]-phot['mag2'][inds_inner2]
    m1 = phot['mag1'][cut_inds]
    if yfilter == Filter2:
        m = m2
    else:
        m = m1
    fig = plt.figure(figsize=(5,5))
    ax = plt.axes()
    ax.plot(m1-m2,m,'o',ms=1.5,mec=None,color='black')    
    print 'HB Bulge Tracks'
    graphics.plot_tracks_on_one_panel_smooth(btrack,ax,Av_bulge,**kwargs)
    majorLocator = MultipleLocator(0.4)
    ax.xaxis.set_major_locator(majorLocator)
    ax.set_xlim(-0.8,1.5)
    ax.set_ylim(25,18.5)
    ax.set_xlabel('$%s-%s$'%(Filter1,Filter2))
    ax.set_ylabel('$%s$'%yfilter)
    plt.savefig('comp_cmd%s'%EXT)
    plt.close('all')
    print 'wrote comp_cmd%s'%EXT
    return

def load_tracks(track_location,track_ext):
    '''
    returns a list of PadovaTracksUtils.Track objects
    for all files that match track_ext in directory track_location.
    ''' 
    track_files = GenUtils.get_afile(track_location,'*'+track_ext)
    tracks = [PadovaTracksUtils.get_tracks(t) for t in track_files]
    return tracks

def test_other_tracks(**kwargs):
    track_mixs = ['C09_AF02_I1_Z0.02_Y0.284',
                  'C09_AF02_I1_Z0.03_Y0.302',
                  'C09_AF02_I1_Z0.04_Y0.321',
                  'C09_AF02_I1_Z0.05_Y0.339',
                  'C09_AF02_I1_Z0.07_Y0.375'][::-1]
    from params import TRACKS_BASE
    track_dir = TRACKS_BASE
    track_ext = 'HB*0.01*UVbright'
    for track_mix in track_mixs:
        track_location = os.path.join(track_dir,track_mix)
        tracks = np.array(load_tracks(track_location,track_ext))
        kwargs['btrack'] = tracks
        time_on_cmd = time_spent_on_cmd(**kwargs)
        print ''
        print track_mix,np.mean(time_on_cmd),np.median(time_on_cmd)
        print ''
        # these are all about 1.1 to 1.2 e7!
        
def time_spent_on_cmd(**kwargs):
    Av_bulge = kwargs.get('Av_bulge')
    af275w = kwargs.get('af275w')
    af336w = kwargs.get('af336w')
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    btrack = kwargs.get('btrack')
    mag2_cut = kwargs.get('mag2_cut')
    dmod = kwargs.get('dmod')
    xlim = (-0.75,1.0)
    ylim = (mag2_cut,18.)
    xwide = .1
    ywide = .5
    xlim_wide = (xlim[0]-xwide,xlim[1]+xwide)
    ylim_wide = (ylim[0]+ywide,ylim[1]-ywide)
    time_on_cmd = []
    plt.figure()
    ax = plt.axes()
    for track in btrack:
        #if not 0.54 <= track.mass <= 0.58: continue
        # print track.mass
        LY = PadovaTracksUtils.Tracks.get_stable_col(track,'LY')
        tmag1 = PadovaTracksUtils.Tracks.get_stable_col(track,Filter1+'1')
        tmag2 = PadovaTracksUtils.Tracks.get_stable_col(track,Filter2+'1')
        if type(tmag1) == int:
            tmag1 = PadovaTracksUtils.Tracks.get_stable_col(track,Filter1)
            tmag2 = PadovaTracksUtils.Tracks.get_stable_col(track,Filter2)
        
        age = PadovaTracksUtils.Tracks.get_stable_col(track,'AGE') 
        HeB = np.nonzero(LY > 0.)[0] 
        
        tM1= tmag1+dmod+Av_bulge*af275w
        tM2= tmag2+dmod+Av_bulge*af336w
        tcol = tM1-tM2
        tmag = tM2
        in_cmd = np.nonzero((tcol < xlim[1]) &
                            (tcol > xlim[0]) &
                            (tmag < ylim[0]) &
                            (tmag > ylim[1]))[0]
        if len(in_cmd) < 2: continue
        # print age[in_cmd]
        print '%.2f enters cmd at %.3g'%(track.mass,age[in_cmd[0]])
        print '%.2f leaves cmd at %.3g'%(track.mass,age[in_cmd[-1]])
        print 'total time on cmd, %.3g'%(age[in_cmd[-1]]-age[in_cmd[0]])
        
        time_on_cmd.append(age[in_cmd[-1]]-age[in_cmd[0]])
        
        ax.plot(tcol[HeB],tmag[HeB],zorder=3)
        ax.plot(tcol[in_cmd],tmag[in_cmd],'.',color='black',zorder=3)
        ax.plot(tcol[in_cmd[0]],tmag[in_cmd[0]],'o',color='blue',zorder=3)
        ax.plot(tcol[in_cmd[-1]],tmag[in_cmd[-1]],'o',color='green',zorder=3)
        ax.annotate('%f'%track.mass,xy=(tcol[np.argmin(tmag[HeB])],min(tmag[HeB])),color='red',**graphics.load_ann_kwargs())
        ax.set_xlim(xlim_wide)
        ax.set_ylim(ylim_wide)
        xs = [xlim[0],xlim[1],xlim[1],xlim[0],xlim[0]]
        ys = [ylim[0],ylim[0],ylim[1],ylim[1],ylim[0]]
        ax.plot(xs,ys,'--',lw=3,color='grey')
    return time_on_cmd

'''
0.48 enters cmd at 1.48e+08
0.48 leaves cmd at 1.77e+08
total time on cmd, 2.87e+07

0.58 enters cmd at 1.41e+08
0.58 leaves cmd at 1.41e+08
total time on cmd, 2.19e+05

0.54 enters cmd at 1.45e+08
0.54 leaves cmd at 1.46e+08
total time on cmd, 9.94e+05

0.50 enters cmd at 1.37e+08
0.50 leaves cmd at 1.6e+08
total time on cmd, 2.31e+07

0.56 enters cmd at 1.43e+08
0.56 leaves cmd at 1.43e+08
total time on cmd, 4.6e+05

ts = [2.19e+05,9.94e+05,2,4.6e+05]
median(ts)
994000.0

mean(ts)
10694600.0
'''
