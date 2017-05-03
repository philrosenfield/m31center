import os
import numpy as np
import pyfits
#from uv_plot_setups import *
from scipy import *
#from PHATDataUtils import points_inside_reg, read_fullbrick, read_basefits
import GenUtils
import data
import file_io
import calcs
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import graphics

def CDFPoisson():
    ax = plt.axes()
    A = 600.
    S = 650000./A
    aas = 1/200.,1/100.,1/50.
    for a in aas:
        ntot = (A/a)*np.array([1.-np.exp(-(S*a))*np.sum([(S*a)**i/scipy.misc.factorial(i) for i in range(int(x+S*a))]) for x in range(25)])
        ns = A/a
        ax.semilogy(ntot,'o',label = r'$a=%.3f^{\prime\prime}\ N_{r}=%i$'%(a,ns))
        
    #ax.set_ylabel(r'$N_{tot} \geq N_s+N_b$')
    ax.set_xlabel(r'$N_s$',fontsize=20)
    #ax.set_ylim(0,1)
    plt.legend(loc=0,frameon=False)
    
def get_ave_blend(**kwargs):
    '''
    kwargs:
    n_unblended = 120,000
    area = 600 arcsec^2
    surf_bright = 18 mag/arcsec^2
    m_unresolved = 26 # ~mag EHB F336W
    '''

    n_unblended = kwargs.get('n_unblended',120000.)
    area = kwargs.get('area',600.)
    surf_bright = kwargs.get('surf_bright',18.)
    m_unresolved = kwargs.get('m_unresolved',26.)
    
    flux_unresolved = 10**(-0.4*m_unresolved)
    star_density = 10**(surf_bright/-2.5)/flux_unresolved
    ave_blend = star_density*area/n_unblended
    return ave_blend

def poisson_method1(**kwargs):
    area = kwargs.get('area')
    n_unblended = kwargs.get('n_unblended',200.*area)
    nstars = kwargs.get('nstars',705.)
    int_flx = kwargs.get('int_flx',98.33)
    flx_stars = kwargs.get('flx_stars',1.03)
    m_unresolved = kwargs.get('m_unresolved')
    expected_mag = kwargs.get('expected_mag',m_unresolved-3.)
    loud = kwargs.get('loud',False)
    # amount of UV flux from unresolved stars (table 2)
    completeness_1 = int_flx/flx_stars
    
    # the amount of flux if the resolved stars were from fainter mags
    completeness_2 = 10**(-0.4*(m_unresolved-expected_mag)) 
    
    # the number of stars that are totally blened
    n_blends = nstars*completeness_1/completeness_2
    
    # number density
    background = n_blends/n_unblended
    
    over_density = background+1./completeness_2
    prob = poisson([background,over_density])
    nstars_expected = n_unblended*prob
    
    if loud == True:
        print 'assuming the mag %.1f stars are blends from mag %.1f stars:'%(expected_mag,m_unresolved)
        print 'area',area
        print 'blended stars',n_blends
        print 'stellar density',n_blends/area
        print 'ave bg',background
        print 'density above bg needed', 1./completeness_2
        print 'prob',prob
        print 'n stars expected from blends over whole region',nstars_expected
        print 'n stars expected from blends arcsec^2',nstars_expected/area
        
    return background
    
def poisson_probability(**kwargs):
    '''
    kwargs used:
    n_unblended = 120,000
    m_unresolved = 26 # ~mag EHB F336W
    mags = np.arange(21.,m_unresolved+.01,0.1)
    
    kwargs passed to ave_blends:
    area = 600 arcsec^2
    surf_bright = 18 mag/arcsec^2
    '''
    import scipy.integrate
    area = kwargs.get('area',600.)
    n_unblended = kwargs.get('n_unblended',200.*area)    
    print 'poss prob nunblended',n_unblended
    m_unresolved = kwargs.get('m_unresolved',26.)
    expected_mag = kwargs.get('expected_mag',m_unresolved-3.)
    mags = kwargs.get('mags',np.arange(16.,m_unresolved+.01,0.1))
        
    #ave_blend = get_ave_blend(**kwargs)
    temp_kwargs = kwargs
    ave_blends = []
    for mag in mags:
        temp_kwargs['expected_mag'] = mag
        #ave_blends.append(get_ave_blend(**temp_kwargs))
        ave_blends.append(poisson_method1(**temp_kwargs))
    factors = 1/10**(-0.4*(m_unresolved-mags))
    over_densities = ave_blends+factors
    
    prob = np.array([])
    for ave_blend,over_density in zip(ave_blends,over_densities):
        p = [ave_blend,over_density]
        prob = np.append(prob,poisson(p))
    
    #scipy.integrate.simps(prob,rnge)
    fins = np.nonzero(np.isfinite(prob))
    nstars_expected = n_unblended*prob[fins]
    print 'mag expected from mag blended',expected_mag,m_unresolved
    #print 'number of stars expected in the entire area', np.sum(nstars_expected[:expected_mag])
    print 'normed and integrated probability, number expected in area:'
    ind, = np.nonzero(np.round(mags[fins],1)==expected_mag)[0]
    tot  = scipy.integrate.simps(prob[fins],mags[fins])
    part = scipy.integrate.simps(prob[fins][:ind],mags[fins][:ind])
    print part/tot,part*n_unblended
    return mags,nstars_expected

def poisson(p):
    '''
    poisson probability
    p[0] = mu
    p[1] = k
    returns (mu^k)e^-mu/k!
    '''
    try:
        scipy.misc.factorial(10.)
    except NameError:
        import scipy.misc
    
    return np.exp(-1.*p[0])*(p[0]**p[1])/scipy.misc.factorial(p[1])

def load_saved_lf(filename='lf_data.dat'):
    lf_lines = open(filename,'r').readlines()
    dic = {}
    for line in lf_lines:
        if line.startswith('#'): continue
        line = line.strip()
        filter = line.split(' ')[0]
        region = line.split(' ')[1]
        data = line[7:]
        bins = map(float,data.split(';')[0].split())
        hist = map(float,data.split(';')[1].split())
        if not filter in dic.keys(): dic[filter] = {}
        if not region in dic[filter].keys(): dic[filter][region] = {}
        #if not 'bins' in dic[filter][region].keys(): dic[filter][region]['bins'] = {}
        dic[filter][region] = {'bins': bins, 'hist': hist}
    return dic
    
def lf_bin_comparison(phot, inds_cut1, inds_cut2, rad_inds, **kwargs):
    lf_data = load_saved_lf('/Users/phil/research/PHAT/code/UVBrights/diy_final/lf_data.dat')
    sup_phot = kwargs.get('sup_phot')
    rif_36 = sup_phot['rif_36']
    comp_file = kwargs.get('comp_file')
    comp_data = kwargs.get('comp_data')
    frac1s, frac2s, magbins = file_io.load_completeness_corrections(comp_file)
    comp50mag1, comp50mag2 = file_io.get_comp50mags(comp_data)
    ylab =r'$\rm{Number}/10^7 M_\odot$'
    xlab1 = '$%s$'%kwargs.get('Filter1')
    xlab2 = '$%s$'%kwargs.get('Filter2')
    cols = kwargs.get('cols')
    filter1 = kwargs.get('Filter1')
    filter2 = kwargs.get('Filter2')
    ax1s = graphics.setup_five_panel_plot(fig=2)
    for i in range(len(rad_inds)):
        ax1 = ax1s[i]
        ax1.set_xlim(26,19)
        ax1.xaxis.set_major_locator(MultipleLocator(2))
        if i == 0: ax1.set_ylabel(ylab)
        ax1.set_xlabel(xlab1)
        # M/L ratio 1e7 scaling is for AXIS
        norm = rif_36[i+1]/2.315e-9*0.95*1e-7
        inds1 = list(set(inds_cut1) & set(rad_inds[i]))
        hist1, bins = np.histogram(phot['mag1'][inds1],bins=magbins)
        known = np.nonzero(bins[1:] <= comp50mag1[i])[0]
        poor =  np.nonzero(bins[1:] >= comp50mag1[i])[0]
        vals = hist1[known]*frac1s[i][known]/norm
        errlow, errhigh = np.transpose([cdf_bitches(v) for v in vals])
        ax1.fill_between(bins[1:][known],errlow,errhigh,color=cols[i],alpha=0.2)
        ax1.plot(bins[1:],hist1*frac1s[i]/norm,
                 drawstyle='steps',color=cols[i],lw=1)
        ax1.plot(bins[1:][known],hist1[known]*frac1s[i][known]/norm,
                 drawstyle='steps',color=cols[i],lw=2)
        lf_hist,lf_bins = lf_data[filter1][str(i)].values()
        ax1.plot(lf_bins,lf_hist,drawstyle='steps',color='black',lw=2)
        ax1.set_yscale('log')
        ax1.set_ylim(.01,10)
    plt.savefig('radialLF_cdf_%s%s'%EXT%kwargs.get('Filter1'))
    plt.close()

    ax2s = graphics.setup_five_panel_plot(fig=1)
    #ax1,ax2 = graphics.two_panel_plot(11,5,xlab1,xlab2,ylab)
    for i in range(len(rad_inds)):
        ax2 = ax2s[i]
        ax2.set_xlim(26,19)
        #ax2.set_ylim(0,8)
        ax2.xaxis.set_major_locator(MultipleLocator(2))
        if i == 0:  ax2.set_ylabel(ylab)
        ax2.set_xlabel(xlab2)
        # M/L ratio 1e7 scaling is for AXIS
        norm = rif_36[i+1]/2.315e-9*0.95*1e-7
        inds2 = list(set(inds_cut2) & set(rad_inds[i]))
        hist2, bins = np.histogram(phot['mag2'][inds2],bins=magbins)
        known = np.nonzero(bins[1:] <= comp50mag2[i])[0]
        poor =  np.nonzero(bins[1:] >= comp50mag2[i])[0]
        vals = hist2[known]*frac2s[i][known]/norm
        errlow, errhigh = np.transpose([cdf_bitches(v) for v in vals])
        ax2.fill_between(bins[1:][known],errlow,errhigh,color=cols[i],alpha=0.2)
        ax2.plot(bins[1:],hist2*frac2s[i]/norm,
                 drawstyle='steps',color=cols[i],lw=1)        
        ax2.plot(bins[1:][known],hist2[known]*frac2s[i][known]/norm,
                 drawstyle='steps',color=cols[i],lw=2)    
        lf_hist,lf_bins = lf_data[filter2][str(i)].values()
        ax2.plot(lf_bins,lf_hist,drawstyle='steps',color='black',lw=2)
        ax2.set_yscale('log')
        ax2.set_ylim(.01,10)
    plt.savefig('radialLF_cdf_%s%s'%EXT%kwargs.get('Filter2'))
    plt.close()
    return
    
def cdf_bitches(val,limits=(0.05,0.95)):
    from scipy.stats import poisson
    import scipy.interpolate
    lamb = np.arange(0,150,0.1)
    vect = np.linspace(0,1,1500)
    func = scipy.special.gammainc(val,lamb)/scipy.misc.common.factorial(scipy.floor(val))
    g = scipy.interpolate.interp1d(func/np.max(func),lamb)
    return g(limits)
    
def normalize_poisson(mu,**kwargs):
    '''
    finds poisson prob over range and integrates 
    returns integration result

    kwargs
    rnge = np.arange(0,100)
    '''
    arr = []
    rnge = kwargs.get('rnge',np.arange(0,100))
    for i in rnge:
        p = [mu,i]
        arr.append(poisson(p))
    norm =  scipy.integrate.simps(arr,rnge)
    return norm

def get_fake(**kwargs): 
    fake = '/astro/net/angst2/philrose/PHAT/paper/referree/test.st'
    d = np.loadtxt(fake,unpack=True)
    
    x = d[2]
    y = d[3]
    Chi = d[4]
    SNR = d[5]
    Sharp = d[6]
    Round = d[7]
    Crowd = d[9]
    Type = d[10]
    mag_F275W = d[15]
    chi_F275W = d[18]
    snr_F275W = d[19]
    sharp_F275W = d[20]
    round_F275W = d[21]
    crowd_F275W = d[22]
    flag_F275W = d[23]
    mag_F336W = d[28]
    chi_F336W = d[31]
    snr_F336W = d[32]
    sharp_F336W = d[33]
    round_F336W = d[34]
    crowd_F336W = d[35]
    flag_F336W = d[36]
    
    cut_F336W = data.quality_cut(sharp_F336W,snr_F336W,crowd_F336W,round_F336W,**kwargs)
    cut_F275W = data.quality_cut(sharp_F275W,snr_F275W,crowd_F275W,round_F275W,**kwargs)
    
    cut_both = list( set(cut_F336W) & set(cut_F275W) )
    return mag_F336W[cut_both]

def add_rectangle_patch(xedge):
    try:
        mpatches
    except NameError:
        import matplotlib.patches as mpatches

    art = mpatches.Rectangle((xedge, 0), 2.9, 100,color='black',alpha=0.2)
    return art
		
if __name__ == "__main__":
    default_file = 'defaults_single_crowd.dat'
    kwargs = file_io.load_defaults(default_file)    
    bulge_fit = file_io.get_bulge_file(st=False)
    phot_cat=os.path.split(bulge_fit)[1]
    phot = file_io.read_data_diy(bulge_fit)
    kwargs['mag1_cut'], kwargs['mag2_cut'] = data.get_mag_cuts(**kwargs)
    sd1, sd2, both, comp, good = data.photometry_cuts(phot,**kwargs)
    
    # divide stars into radial bins [0] = inner [-1] = outter
    from params import PHAT_DIR_EXTRA3
    contour = PHAT_DIR_EXTRA3 + 'contours_v3.reg'
    rad_inds = calcs.points_inside_ds9_polygon(contour,phot['ra'],phot['dec'])
    
    reg1 = list(set(good) & set(rad_inds[0]))
    F336W = phot['mag2'][reg1]
    
    fake336W = get_fake(**kwargs)
    
    f336w_range = np.arange(18.9,30,0.1)
    offset = 0.01
    f336w_bins = np.arange(19.,27.,offset)
    
    reg1_hist,reg1_bins = np.histogram(F336W,bins=f336w_range)
    hist,bins = np.histogram(fake336W,bins=f336w_range)
    
    pkwargs = {'area' : 600.,'m_unresolved': 26., 'surf_bright': 18,'expected_mag':23.}
    print 'obs stars brighter than %.1f: %i'%(pkwargs['expected_mag'],
                                            len(np.nonzero(F336W<=pkwargs['expected_mag'])[0]))
    mags, nstars_expected = poisson_probability(**pkwargs)
    
    print ''
    pkwargs1 = {'area' : 600.,'m_unresolved': 26.5, 'surf_bright': 18, 'expected_mag':24.,'nstars':280}
    print 'obs stars brighter than %.1f: %i'%(pkwargs1['expected_mag'],
                                            len(np.nonzero(F336W<=pkwargs1['expected_mag'])[0]))
    mags1, nstars_expected1 = poisson_probability(**pkwargs1)
    
    #poisson_method1(**pkwargs)
    #poisson_method1(**pkwargs1)
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(reg1_bins[1:],reg1_hist/pkwargs['area'],label=r'$\rm{Region\ 1\ Sources}$',
            color='green',lw=2,ls='steps')
    ax.plot(mags,nstars_expected/pkwargs['area'],
            label=r'$\rm{Poisson\ Expectation}$ $m_{F336W}=%.1f$'%pkwargs['m_unresolved'],
            color='black',lw=2)
    ax.plot(mags1,nstars_expected1/pkwargs1['area'],'--',
            label=r'$\rm{Poisson\ Expectation}$ $m_{F336W}=%.1f$'%pkwargs1['m_unresolved'],
            color='black',lw=2)   
    ax.plot(bins[1:],hist/pkwargs['area'],
            label=r'$\rm{Simulated\ Field}$ $m_{F336W}=26.0$',
            color='blue',lw=2,ls='steps')
    ax.add_patch(add_rectangle_patch(24.6))
    ax.set_xlabel(r'$F336W$')
    ax.set_ylabel(r'$\rm{\#\ arcsec}^{-2}$')
    ax.set_xlim(26,19)
    ax.set_ylim(0,.4)
    plt.legend(loc=0,frameon=False)
    plt.savefig('crowdtest2%s'%EXT)    
    '''
    print 'doing region 2'
    # test region 2.
    r2kwargs= {'int_flx': 351.08,
               'flx_stars': 2.41,
               'surf_bright': 18.36,
               'area': 3197.,
               'm_unresolved': 26.,
               'nstars': 1491}
    reg2 = list(set(good) & set(rad_inds[1]))
    mags, nstars_expected = poisson_probability(**r2kwargs)
    F336W = phot['mag2'][reg2]
    #plt.show()
    #ax.set_yscale('log')
    '''
    '''
    m_tests = np.arange(26.3,27,0.1)
    for m_test in m_tests:
        kwargs['m_unresolved'] = m_test
        mags, nstars_expected = poisson_probability(**kwargs)
        ax.plot(mags,nstars_expected/kwargs['area'],
                label=r'm test=%f'%m_test,lw=2)
    '''
