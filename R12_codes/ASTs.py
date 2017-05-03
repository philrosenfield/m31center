import matplotlib.pyplot as plt
import os
import sys
#from kapteyn import wcs
import pyfits
import numpy as np
import matplotlib.nxutils as nxutils
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import NullFormatter, MultipleLocator
nullfmt   = NullFormatter() # no labels
import scipy.interpolate
import graphics
from data import quality_cut
from GenUtils import get_afile,smooth,closest_match
import file_io
import calcs
from params import *

def ast_diff_test(ast_rec,ast_in,axs=None,color='black',spacing = 0.1):
    ast_diff = ast_rec-ast_in
    
    if ax==None:
        axs = graphics.two_panel_plot_vert_bsquat()
        axs[0].set_xlabel('Recovered')
        axs[1].set_ylabel('standard deviation')
        axs[0].set_ylabel('residual')
    bins = np.arange(np.min(ast_rec),np.max(ast_rec)+spacing,spacing)
    inds = np.digitize(ast_rec,bins)
        
    faints = np.nonzero(ast_diff>=0)[0]
    brights = np.nonzero(ast_diff<=0)[0]
    
    fstdevs = np.array([])
    bstdevs = np.array([])
    filled_bins = np.array([])
    for i in range(len(bins)):
        ibin = np.nonzero(inds==i)[0]
        if len(ibin) == 0: continue
        fs = list(set(faints) & set(ibin))
        bs = list(set(brights) & set(ibin))
        fstdevs = np.append(fstdevs,np.std(ast_diff[fs]))
        bstdevs = np.append(bstdevs,np.std(ast_diff[bs]))
        filled_bins = np.append(filled_bins,bins[i])
    
    axs[1].plot(filled_bins,fstdevs,color=color)
    axs[1].plot(filled_bins,bstdevs,color=color,lw=2)
    axs[0].plot(filled_bins,bstdevs-fstdevs,color=color,ls='--')
    return

def how_does_the_lf_vary():
    lf_o = CrowdTest.load_saved_lf(filename='/Users/phil/research/PHAT/code/UVBrights/diy_final/lf_data.dat')    
    filters = ('F275W','F336W')
    vary = np.array([])
    for filter in filters:
        sortkeys = np.sort(lf_o[filter].keys())
        vary = np.array([])
        for region in sortkeys:
            hist,bins = lf_o[filter][region].values()
            hist = np.array(hist)
            bins = np.array(bins)
            maxbin = bins[np.argmax(hist)]
            maxhist = np.max(hist)
            filled = np.nonzero(hist)[0]
            minbin = bins[filled][0]
            minhist =hist[filled][0]
            magdiff =  maxbin - minbin
            print filter, region, (maxhist-minhist)/magdiff
            vary = np.append(vary,maxhist/magdiff)
        print '    ',filter,np.median(vary),np.mean(vary)
def power_law(logx,logy,logyerr):
    ''''
    y = ax**k
    '''
    powerlaw = lambda x, amp, index: amp * (x**index)
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
    #x = np.logspace(-1,6,len(yvals))
    pinit = [1.0, -1.0]
    out = scipy.optimize.leastsq(errfunc, pinit,args=(logx, logy, logyerr), full_output=1)
    pfinal = out[0]
    covar = out[1]
    
    index = pfinal[1]
    amp = 10.0**pfinal[0]
    
    indexErr = sqrt( covar[0][0] )
    ampErr = sqrt( covar[1][1] ) * amp
    print index, amp
    return index, amp

''' tests
            hist_o, bins_o = lf_o[filter][region].values()
            hist_o = np.asfarray(hist_o)
            bins_o = np.asfarray(bins_o)
            filled = np.nonzero(hist_o)[0]
            logx = np.log10(bins_o[filled])
            logy = np.log10(hist_o[filled])
            yerr = 0.2 * hist_o[filled]
            logyerr = yerr / hist_o[filled]
            hist_o = np.asfarray(hist_o)
            index,amp = power_law(logx,logy,logyerr)
            #powerlaw(bins_o, amp, index)
            #axs[k].plot(bins,hist,ls='steps',color=cols[k])
            #plt.plot(bins,hist,ls='steps',color=cols[k])
            #plt.plot(bins_o,np.array(hist_o)*1e3,ls='steps',color=cols[k])
            plt.plot(bins_o,powerlaw(bins_o, amp, index))
            plt.plot(bins_o,hist_o)
            plt.plot(bins,powerlaw(bins, amp, index))
            plt.yscale('log')
'''
    
def ast_LF(ast,**kwargs):
    cont_file = kwargs.get('cont_file')
    one_plot = kwargs.get('one_plot',True)
    powerlaw = lambda x, amp, index: amp * (x**index)
    import CrowdTest
    cols = kwargs.get('cols')
    double_plot = kwargs.get('double_plot',False)
    two_panel = kwargs.get('two_panel',True)
    lf_o = CrowdTest.load_saved_lf(filename='/Users/phil/research/PHAT/code/UVBrights/diy_final/lf_data.dat')
    comp_data = kwargs['comp_data']
    mag1_cut =  file_io.get_comp90mags(comp_data)[0][0]
    mag2_cut = file_io.get_comp90mags(comp_data)[1][0]
    print '# mag1_cut, mag2_cut', mag1_cut, mag2_cut
    if cont_file != None:
        cf = open(cont_file,'w')
        cf.write('# mag1_cut, mag2_cut %.1f %.1f\n'%(mag1_cut, mag2_cut))
    from params import PHAT_DIR_EXTRA3
    contour = PHAT_DIR_EXTRA3 + 'contours_v3.reg'
    inds = calcs.points_inside_ds9_polygon(contour,ast['ra'],ast['dec'])
    filters = ('F275W','F336W')
    
    rec1 = recovered_asts(ast['mag1in'],ast['mag1out'])
    rec2 = recovered_asts(ast['mag2in'],ast['mag2out'])
        
    # quality cut on all asts
    qc1 = quality_cut(ast['sharp1'],
                      ast['snr1'], 
                      ast['crowd1'],
                      ast['round1'],
                      loud=True,
                      **kwargs)
    
    qc2 = quality_cut(ast['sharp2'],
                      ast['snr2'],
                      ast['crowd2'],
                      ast['round2'],
                      loud=True,
                      **kwargs)
    spacing = .05
    if two_panel == True:
        ax1,ax2 = graphics.two_panel_plot(11,5,'x','x','$\#$')
        axs = [ax1,ax2,ax1,ax2,ax1]
    for j,filter in enumerate(filters):
        fcolor = 'blue'
        if double_plot == True:
            axs,baxs =  graphics.setup_ten_panel_plot(fig=j+1,noy=True)        
        elif two_panel == False:
            fig = j+1
            if one_plot == True: fig = 1
            axs = graphics.setup_five_panel_plot(fig=fig,noy=True)
            if j != 0: fcolor = 'green'
        sortkeys = np.sort(lf_o[filter].keys())
        for k,region in enumerate(sortkeys):
            if filter == 'F336W': 
                qc = qc2
                ast_mag_in = ast['mag2in']
                ast_mag_out = ast['mag2out']
                mag = mag2_cut
                if two_panel == True: axs[k] = ax2
                axs[k].set_xlim(26.1,16)
                index = 16.943523372280001
                amp = 1e-17
            if filter == 'F275W': 
                qc = qc1
                ast_mag_in = ast['mag1in']
                ast_mag_out = ast['mag1out']
                mag = mag1_cut
                if two_panel == True: axs[k] = ax1
                axs[k].set_xlim(25.9,16)
                index = 19.246755161139998
                amp = 1e-20
                
            Rinds = list(set(rec1) & set(qc) & set(inds[k]))
            bins = np.arange(np.min(ast_mag_in[Rinds]),np.max(ast_mag_in[Rinds]),spacing)
            hist = powerlaw(bins, 1e-17, 15)
            #hist = powerlaw(bins,amp,index)
            # paste tests here
            axs[k].plot(bins,hist,ls='steps',color='black')
            hist_o, bins_o = lf_o[filter][region].values()
            #axs[k].plot(bins_o,np.array(hist_o)*1e2,ls='steps',color=cols[k])
            # bin up input mag
            in1s = np.digitize(ast_mag_in[Rinds],bins)
            new_mags = np.array([])
            for i in range(len(bins)):
                if hist[i] == 0: continue
                if len(np.array(in1s==i).nonzero()[0]) <= 1: continue
                # output mag corresponding to the input mag bin
                ast_mag = ast_mag_out[Rinds][np.array(in1s==i).nonzero()[0]]
                # make cumulative distribution
                sort_ast = np.sort(ast_mag)
                yvals = np.arange(len(sort_ast))/float(len(sort_ast))
                # interpolate it from 0,1
                f = scipy.interpolate.interp1d(yvals,sort_ast,bounds_error=False)
                # pick randomly from distribution the same number as in the current histogram
                new_mags = np.append(new_mags,f(np.random.random(hist[i])))
            # find new LF
            new_hist,lixo =  np.histogram(new_mags,bins=bins)
            #print sum(new_hist), sum(hist)
            axs[k].plot(bins[1:],new_hist,ls='steps',color=fcolor,lw=2)
            ofst=0.1
            #mag = closest_match(bins[1:],(new_hist-hist[1:])/hist[1:],.2)
            g = scipy.interpolate.interp1d(bins[1:],(new_hist-hist[1:])/hist[1:])
            print  filter,k,g(mag)
            if cont_file != None:
                cf.write('%s %i %f\n'%(filter,k,g(mag)))
            if double_plot == True:
                baxs[k].plot(np.zeros(10)+mag,np.linspace(baxs[k].get_ylim()[0],baxs[k].get_ylim()[1],10),'--',color='black',label='%.2f'%g(mag))
                baxs[k].set_xlim(np.max(ast_mag_in)-ofst,np.min(ast_mag_in)+ofst)
                baxs[k].legend(loc=0,frameon=False)
                baxs[k].set_xlabel('$%s$'%filter)
                baxs[k].plot(bins[1:],(new_hist-hist[1:])/hist[1:],color='black',lw=2,ls='steps')
            else:
                if one_plot == True:
                    axs[k].set_xlabel(r'$\rm{mag}$')
                    if k!=0: axs[k].yaxis.set_major_formatter(nullfmt)
                else:
                    axs[k].set_xlabel('$%s$'%filter)
            axs[k].xaxis.set_major_locator(MultipleLocator(2))
            axs[k].xaxis.set_minor_locator(MultipleLocator(1))
            axs[k].set_ylim(1e1,1e5)
            axs[k].plot(np.zeros(3)+mag,np.logspace(1,5,3),'--',lw=2,color=fcolor,label='%.2f'%g(mag))
            axs[k].set_yscale('log')
        axs[0].set_ylabel('$\#$')
        axs[k].annotate(r'$F275W$',xy=(22,5e4),fontsize=20,color='blue')        
        axs[k].annotate(r'$F336W$',xy=(22,3e4),fontsize=20,color='green')
        [ax.set_ylim(1e2,1e5) for ax in axs]
        [ax.set_xlim(26,19) for ax in axs]
        if double_plot == True: baxs[0].set_ylabel('(input-obs)/input')
        if one_plot !=True:
            plt.figure(j+1)
            plt.savefig('ast_lf_test_%s%s'%EXT%filter)
        else:
            plt.savefig('ast_lf_test%s'%EXT)
    if cont_file != None: cf.close()
    return axs

def closest_match(x,y,val):
    g = scipy.interpolate.interp1d(x,y)
    arr = np.arange(x[0],x[-1],1e-4)
    t = np.array([g(a) for a in arr])
    ind = np.argmin(abs(val-t[np.isfinite(t)]))
    return arr[ind]
    
def new_mag_cut(**kwargs):
    '''
    this can run by a call to ASTs, not just by running main
    set the magnitude cutoff to be the faintest recovered magnitude of a 
    star whose input magnitude was at the 50% completeness limit
    '''
    ast_file = kwargs.get('ast_file')
    comp_data = kwargs.get('comp_data')
    
    ast = load_ast_file(ast_file)
    comp = file_io.load_comp_data(comp_data)
    pm = 0.05
    
    # within dmag range
    rec1 = recovered_asts(ast['mag1in'],ast['mag1out'])
    rec2 = recovered_asts(ast['mag2in'],ast['mag2out'])
        
    # quality cut on all asts
    qc1 = quality_cut(ast['sharp1'],
                      ast['snr1'], 
                      ast['crowd1'],
                      ast['round1'],
                      loud=True,
                      **kwargs)
    
    qc2 = quality_cut(ast['sharp2'],
                      ast['snr2'],
                      ast['crowd2'],
                      ast['round2'],
                      loud=True,
                      **kwargs)
    # each radial bin
    from params import PHAT_DIR_EXTRA3
    contour = PHAT_DIR_EXTRA3 + 'contours_v3.reg'
    inds = points_inside_ds9_polygon(contour,ast['ra'],ast['dec'])
    filters = ('F275W','F336W')
    ax1s =  graphics.setup_five_panel_plot(fig=1)   
    [ax.set_xlabel('$%s$ Recovered'%filters[0]) for ax in ax1s]
    ax1s[0].set_ylabel('$%s$ Recovered-Input'%filters[0])
    cols = kwargs.get('cols')
    axs = graphics.two_panel_plot_vert_bsquat()
    axs[0].set_xlabel('F275W Recovered')
    axs[1].set_ylabel('standard deviation')
    axs[0].set_ylabel('residual')    
    for i in range(len(inds)):
        Rinds1 = list(set(rec1) & set(qc1) & set(inds[i]))
        #Rinds1 = list(set(rec1) & set(inds[i]))
        #Rinds1 = list(set(rec1) & set(qc1) & set(qc2) & set(inds[i]))
        within1 = np.nonzero((ast['mag1in'][Rinds1]<=comp['comp50mag1'][i]+pm) & 
                             (ast['mag1in'][Rinds1]>=comp['comp50mag1'][i]-pm))[0]
        mag1_cut = np.max(ast['mag1out'][Rinds1][within1])
        print mag1_cut
        ast_diff_test(ast['mag1out'][Rinds1],ast['mag1in'][Rinds1],axs=axs,color=cols[i],spacing=0.4)
        axs[0].set_ylim(-0.01,0.1)
        ax1s[i].plot(ast['mag1out'][inds[i]],ast['mag1out'][inds[i]]-ast['mag1in'][inds[i]],'.',color='black')
        ax1s[i].plot(ast['mag1out'][Rinds1],ast['mag1out'][Rinds1]-ast['mag1in'][Rinds1],'.',color='blue')
        #ax1s[i].plot(ast['mag1in'][inds[i]],ast['mag1out'][inds[i]]-ast['mag1in'][inds[i]],'.',color='black')
        #ax1s[i].plot(ast['mag1in'][Rinds1],ast['mag1out'][Rinds1]-ast['mag1in'][Rinds1],'.',color='blue')
        ax1s[i].plot(np.zeros(100)+mag1_cut,np.linspace(-2,10,100),color='green')
        ax1s[i].plot(np.zeros(100)+comp['comp50mag1'][i],np.linspace(-2,10,100),color='red')
        ax1s[i].xaxis.set_major_locator(MultipleLocator(2))
        ax1s[i].set_xlim(19,27) 
        ax1s[i].set_ylim(-2,2)
    plt.savefig('ast_diff_mag1%s'%EXT)
    
    ax2s =  graphics.setup_five_panel_plot(fig=2)
    [ax.set_xlabel('$%s$ Recovered'%filters[1]) for ax in ax2s]
    ax2s[0].set_ylabel('$%s$ Recovered-Input'%filters[1])
    axs = graphics.two_panel_plot_vert_bsquat()
    axs[0].set_xlabel('F336W Recovered')
    axs[1].set_ylabel('standard deviation')
    axs[0].set_ylabel('residual')
    for i in range(len(inds)):
        Rinds2 = list(set(rec2) & set(qc2) & set(inds[i]))
        #Rinds2 = list(set(rec2) & set(qc1) & set(qc2) & set(inds[i]))
        within2 = np.nonzero((ast['mag2in'][Rinds2]<=comp['comp50mag2'][i]+pm) & 
                             (ast['mag2in'][Rinds2]>=comp['comp50mag2'][i]-pm))[0]
        mag2_cut = np.max(ast['mag2out'][Rinds2][within2])
        print mag2_cut
        ast_diff_test(ast['mag2out'][Rinds2],ast['mag2in'][Rinds2],axs=axs,color=cols[i],spacing=0.4)
        axs[0].set_ylim(-0.01,0.1)
        ax2s[i].plot(ast['mag2out'][inds[i]],ast['mag2out'][inds[i]]-ast['mag2in'][inds[i]],'.',color='black')
        ax2s[i].plot(ast['mag2out'][Rinds2],ast['mag2out'][Rinds2]-ast['mag2in'][Rinds2],'.',color='blue')
        ax2s[i].plot(np.zeros(100)+mag2_cut,np.linspace(-2,10,100),color='green')
        ax2s[i].plot(np.zeros(100)+comp['comp50mag2'][i],np.linspace(-2,10,100),color='red')
        ax2s[i].xaxis.set_major_locator(MultipleLocator(2))
        ax2s[i].set_xlim(19,27)
        ax2s[i].set_ylim(-2,2)
    plt.savefig('ast_diff_mag2%s'%EXT)
    return
    
def hist_bin_optimization(arr,limit = None):
    '''
    I got this algorythm from 
    http://toyoizumilab.brain.riken.jp/hideaki/res/histogram.html
    '''
    costs = np.array([])
    tol = 1000.
    limits = (arr.min(),arr.max())
    delta_max = np.diff(limits)/5.
    if limit != None: delta_max = limit
    delta_min = delta_max/tol
    deltas = np.arange(delta_min,delta_max,delta_min)
    for delta in deltas:
        hist_range = np.arange(arr.min(),arr.max(),delta)
        hist,bin_edges = np.histogram(arr,bins=hist_range)
        k = np.mean(hist)
        v = np.var(hist)
        cost = (2*k-v)/delta**2
        costs = np.append(costs,cost)
        
    cs = smooth(costs)
    ds = smooth(deltas)
    return costs,deltas,ds[np.argmin(cs)]
    
def test_hist_bins(arr,limit=None):
    ''' This is meant to be done in ipython, it's just for copy and pasting ease'''
    cs,ds,sp = hist_bin_optimization(arr,limit=limit)
    plt.plot(ds,cs)
    
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

def points_inside_ds9_polygon(reg_name,ra_points,dec_points):
    '''
    uses read_reg
    '''
    radec = np.column_stack((ra_points,dec_points))
    ras,decs = read_reg(reg_name,shape='polygon')
    verts = [np.column_stack((ras[i],decs[i])) for i in range(len(ras))]
    masks = [nxutils.points_inside_poly(radec, vert) for vert in verts]
    
    inds = [np.nonzero(masks[i]-masks[i+1])[0] for i in range(len(masks)-1)]
    return inds

def read_all_asts(verts):
    '''
    Reads dolphot ast output. 20 columns of photometry data + 170 columns of dolphot out.
    verts is the vertices of some region where to kick out the rest of the ASTs.
    '''    
    ra = dec = mag1in = mag2in = objtype = np.array([])
    mag1out = mag1err = snr1 = round1 = crowd1 = sharp1 = qflag1 = np.array([])
    mag2out = mag2err = snr2 = round2 = crowd2 = sharp2 = qflag2 = np.array([])
    
    names = ['X','Y','mag1in','mag2in','objtype','mag1out','mag1err','snr1','sharp1','round1','crowd1','qflag1','mag2out','mag2err','snr2','sharp2','round2','crowd2','qflag2']
    
    allfakes = get_afile('/astro/net/angst2/philrose/PHAT/simulated_data/','*fake')
    for fake in allfakes:
        #print '# '+fake
        drz = fake.replace('all.fake','F336W_drz.chip1.fits')
        proj = wcs.Projection(pyfits.open(drz)[0].header)
        
        t = np.genfromtxt(fake,unpack=True,
                          usecols = (2,3,5,15,30,33,35,37,38,39,40,46,49,51,53,54,55,56,62),
                          names=names)

        tra,tdec = proj.toworld((t['X'],t['Y']))
        points = np.column_stack((tra,tdec))
        inds = np.nonzero(nxutils.points_inside_poly(points, verts))[0]
        ra = np.append(ra, tra[inds])
        dec = np.append(dec, tdec[inds])
        mag1in = np.append(mag1in,t['mag1in'][inds])
        mag2in = np.append(mag2in,t['mag2in'][inds])
        objtype = np.append(objtype,t['objtype'][inds])
        mag1out = np.append(mag1out,t['mag1out'][inds])
        mag1err = np.append(mag1err,t['mag1err'][inds])
        snr1 = np.append(snr1,t['snr1'][inds])
        sharp1 = np.append(sharp1,t['sharp1'][inds])
        round1 = np.append(round1,t['round1'][inds])
        crowd1 = np.append(crowd1,t['crowd1'][inds])
        qflag1 = np.append(qflag1,t['qflag1'][inds])
        mag2out = np.append(mag2out,t['mag2out'][inds])
        mag2err = np.append(mag2err,t['mag2err'][inds])
        snr2 = np.append(snr2,t['snr2'][inds])
        sharp2 = np.append(sharp2,t['sharp2'][inds])
        round2 = np.append(round2,t['round2'][inds])
        crowd2 = np.append(crowd2,t['crowd2'][inds])
        qflag2 = np.append(qflag2,t['qflag2'][inds])
        
    ast={'ra': ra,
         'dec': dec,
         'mag1in': mag1in,
         'mag2in': mag2in,
         'objtype': objtype,
         'mag1out': mag1out,
         'mag1err': mag1err,
         'snr1': snr1,
         'sharp1': sharp1,
         'round1': round1,
         'crowd1': crowd1,
         'qflag1': qflag1,
         'mag2out': mag2out,
         'mag2err': mag2err,
         'snr2': snr2,
         'sharp2': sharp2,
         'round2': round2,
         'crowd2': crowd2,
         'qflag2': qflag2}
    return ast

def write_asts(ast_dict,filename,**kwargs):
    o = open(filename,'w')
    sharp_cut = kwargs.get('sharp_cut',0.075)
    snr_cut = kwargs.get('snr_cut',4.)
    crowd_cut = kwargs.get('crowd_cut',None)
    round_cut = kwargs.get('round_cut',1.5)
    o.write('# abs(sharp)<%f snr>%f abs(round)<%f crowd<%s\n'%(sharp_cut,
                                                       snr_cut,
                                                       round_cut,
                                                       str(crowd_cut)))
    o.write('# ra, dec, mag1in, mag2in, objtype, mag1out, mag1err, snr1, sharp1, round1, crowd1, qflag1, mag2out, mag2err, snr2, sharp2, round2, crowd2, qflag2 \n')
    fmt = '%.8f %.8f %.3f %.3f %i %.3f %.3f %.1f %.3f %.3f %.3f %.1f %.3f %.3f %.1f %.3f %.3f %.3f %.1f\n'
    for i in range(len(ast['ra'])):
        o.write(fmt%(ast_dict['ra'][i],ast_dict['dec'][i],ast_dict['mag1in'][i],
        ast_dict['mag2in'][i],ast_dict['objtype'][i],
        ast_dict['mag1out'][i],ast_dict['mag1err'][i],ast_dict['snr1'][i],ast_dict['sharp1'][i],
        ast_dict['round1'][i],ast_dict['crowd1'][i],ast_dict['qflag1'][i],
        ast_dict['mag2out'][i],ast_dict['mag2err'][i],ast_dict['snr2'][i],
        ast_dict['sharp2'][i],ast_dict['round2'][i],ast_dict['crowd2'][i],ast_dict['qflag2'][i]))
    return filename
    
def recovered_asts(magin,magout,dmaglim=0.7):
    return np.nonzero(abs(magin-magout) < dmaglim)[0]

def load_ast_file(filename):
    ra, dec, mag1in, mag2in, objtype, mag1out, mag1err, snr1, sharp1, round1, crowd1, qflag1, mag2out, mag2err, snr2, sharp2, round2, crowd2, qflag2 = np.loadtxt(filename,unpack=True)
    ast={'ra': ra,
         'dec': dec,
         'mag1in': mag1in,
         'mag2in': mag2in,
         'objtype': objtype,
         'mag1out': mag1out,
         'mag1err': mag1err,
         'snr1': snr1,
         'sharp1': sharp1,
         'round1': round1,
         'crowd1': crowd1,
         'qflag1': qflag1,
         'mag2out': mag2out,
         'mag2err': mag2err,
         'snr2': snr2,
         'sharp2': sharp2,
         'round2': round2,
         'crowd2': crowd2,
         'qflag2': qflag2}
    return ast

def ast_diff_plot(ast,Rinds1s,Rinds2s,inds,**kwargs):
    ast_file_ext = kwargs.get('ast_file_ext','')
    filters = ('F275W','F336W')
    ax1s =  graphics.setup_five_panel_plot(fig=1)
    [ax.set_xlabel('$%s$ Recovered'%filters[0]) for ax in ax1s]
    ax1s[0].set_ylabel('$%s$ Recovered-Input'%filters[0])

    mag1in = ast['mag1in']
    mag1out = ast['mag1out']
    mag2in = ast['mag2in']
    mag2out = ast['mag2out']
    for j in range(len(Rinds1s)):
        ax1s[j].plot(mag1out[inds[j]],mag1out[inds[j]]-mag1in[inds[j]],'.',color='black')
        ax1s[j].plot(mag1out[Rinds1s[j]],
                     mag1out[Rinds1s[j]]-mag1in[Rinds1s[j]],'.',color='blue')
        ax1s[j].xaxis.set_major_locator(MultipleLocator(2))
        ax1s[j].set_xlim(19,27)
        ax1s[j].set_ylim(-2,2)
    plt.savefig('ast_diff_mag1%s%s'%EXT%(ast_file_ext))
    plt.close()
    
    ax2s =  graphics.setup_five_panel_plot(fig=2)
    [ax.set_xlabel('$%s$ Recovered'%filters[1]) for ax in ax2s]
    ax2s[0].set_ylabel('$%s$ Recovered-Input'%filters[1])
    for j in range(len(Rinds2s)):
        ax2s[j].plot(mag2out[inds[j]],mag2out[inds[j]]-mag2in[inds[j]],'.',color='black')
        ax2s[j].plot(mag2out[Rinds2s[j]],
                     mag2out[Rinds2s[j]]-mag2in[Rinds2s[j]],'.',color='blue')
        ax2s[j].xaxis.set_major_locator(MultipleLocator(2))
        ax2s[j].set_xlim(19,27)
        ax2s[j].set_ylim(-2,2)
    
    plt.savefig('ast_diff_mag2%s%s'%EXT%(ast_file_ext))
    plt.close()
    return

def ast_completeness_plot(Bins,Frac1s,Frac2s,**kwargs):
    ast_file_ext = kwargs.get('ast_file_ext','')
    cols = kwargs.get('cols')
    ax1, ax2 = graphics.two_panel_plot(10,5,'$F275W$','$F336W$','completeness',fig=6)
    for i in range(len(Frac1s)):
        ax1.plot(Bins[i][1:],1/Frac1s[i],ls='steps',label='%i'%i,color=cols[i])
        ax2.plot(Bins[i][1:],1/Frac2s[i],ls='steps',label='%i'%i,color=cols[i])

    for ax in [ax1,ax2]:
        ax.grid()
        ax.set_ylim(0,1)
        ax.plot(Bins[0][1:],Bins[0][1:]*0+0.5,lw=2,color='gray') 
        ax.plot(Bins[0][1:],Bins[0][1:]*0+0.9,lw=2,color='gray')
        ax.set_xlim(19,27)
        ax.legend(loc=0,frameon=False)
    
    plt.savefig('completeness%s%s'%EXT%ast_file_ext)
    #plt.close()
    return ax1,ax2

def ast_hist_plot(Bins,Rhist1s,Qhist1s,Rhist2s,Qhist2s,**kwargs):
    cols = kwargs.get('cols')
    ast_file_ext = kwargs.get('ast_file_ext','')
    ax3, ax4 = graphics.two_panel_plot(10,5,'$F275W$','$F336W$','ASTs',fig=7)
    for i in range(len(Rhist1s)):
        ax3.plot(Bins[i][1:],Rhist1s[i],ls='steps',
                 label='%i Recovered'%i,color=cols[i])
        ax3.plot(Bins[i][1:],Qhist1s[i],ls='steps--',
                 color=cols[i])
        ax4.plot(Bins[i][1:],Rhist2s[i],ls='steps',
                 label='%i Recovered'%i,color=cols[i])
        ax4.plot(Bins[i][1:],Qhist2s[i],ls='steps--',
                 color=cols[i])
    
    for ax in [ax3,ax4]:
        ax.set_xlim(19,27)
        ax.legend(loc=0,frameon=False)

    plt.savefig('ast_hist%s%s'%EXT%ast_file_ext)
    plt.close()
    return ax3,ax4

def ast_test_lasttime():
    write_comps = kwargs.get('write_comps',False)
    ast_file = kwargs.get('ast_file')
    ast = load_ast_file(ast_file)
    from params import PHAT_DIR_EXTRA3
    contour = PHAT_DIR_EXTRA3 + 'contours_v3.reg'
    cras, cdecs = file_io.read_reg(contour)
    outverts = np.column_stack((cras[-1],cdecs[-1]))
    
    rec1 = recovered_asts(ast['mag1in'],ast['mag1out'])
    rec2 = recovered_asts(ast['mag2in'],ast['mag2out'])
    
    # quality cut on all asts
    qc1 = quality_cut(ast['sharp1'],
                      ast['snr1'], 
                      ast['crowd1'],
                      ast['round1'],
                      **kwargs)
    
    qc2 = quality_cut(ast['sharp2'],
                      ast['snr2'],
                      ast['crowd2'],
                      ast['round2'],
                      **kwargs)
    
    qcd = list(set(qc1) & set(qc2))
    inds = points_inside_ds9_polygon(contour,ast['ra'],ast['dec'])
    if write_comps == True:
        comp_corr_file = open(kwargs.get('comp_file'),'w') 
        comp_corr_file.write('# region (0=in) mags comp_frac \n')
        comp_data_file = open(kwargs.get('comp_data'),'w') 
        print 'comp correction file:',kwargs.get('comp_file')
        print 'comp data file:',comp_data_file
    single = kwargs.get('single')
    kwargs['ast_file_ext'] = '_'+kwargs.get('comp_file').replace('.dat','').split('_')[-1]
    if single == False:
        qc1 = qcd
        qc2 = qcd
        kwargs['ast_file_ext'] += '_dd'
        color_range = (-0.7,0.5)
        ast_color = ast['mag1in']-ast['mag2in']
        color_inds = np.nonzero( (ast_color>color_range[0]) &
                                 (ast_color<color_range[1]))[0]
        
    Rinds1s , Rinds2s , Bins , Frac1s , Frac2s = [], [], [], [], []
    Rhist1s , Rhist2s , Qhist1s , Qhist2s = [], [], [], []
    if write_comps == True: comp_data_file.write('# Region comp50mag1 comp50mag2 comp90mag1 comp90mag2\n')
    print '# Region comp50mag1 comp50mag2 comp90mag1 comp90mag2'
    for i in range(len(inds)):
        # in each annulus, passed qc, recovered
        if single == False:
            Rinds1 = list(set(rec1) & set(qc1) & set(inds[i]) &
                          set(color_inds))
            Rinds2 = list(set(rec2) & set(qc2) & set(inds[i]) &
                          set(color_inds))
            # REWRITING INDS OF I FOR THIS COLOR RANGE.
            inds[i] = list(set(inds[i]) & set(color_inds))
        else:
            Rinds1 = list(set(rec1) & set(qc1) & set(inds[i]))
            Rinds2 = list(set(rec2) & set(qc2) & set(inds[i]))
        
        hist_range = (18.34,28.)   
        # use hist_bin_optimization to get this value
        #spacing = 0.04 
        spacing = 0.2 
        bins = np.arange(22,hist_range[1],spacing)
        Rhist1,bins = np.histogram(ast['mag1in'][Rinds1],bins=bins)
        Rhist2,bins = np.histogram(ast['mag2in'][Rinds2],bins=bins)
        
        Qhist1,bins = np.histogram(ast['mag1in'][inds[i]],bins=bins)
        Qhist2,bins = np.histogram(ast['mag2in'][inds[i]],bins=bins)
        
        # need to be float...
        arr = np.arange(bins[1],bins[-1],1e-4)
        R1 = scipy.interpolate.interp1d(bins[1:],np.asfarray(Rhist1),bounds_error=False)
        R2 = scipy.interpolate.interp1d(bins[1:],np.asfarray(Rhist2),bounds_error=False)
        Q1 = scipy.interpolate.interp1d(bins[1:],np.asfarray(Qhist1),bounds_error=False)
        Q2 = scipy.interpolate.interp1d(bins[1:],np.asfarray(Qhist2),bounds_error=False)
        
        frac1 = R1(arr)/Q1(arr)
        frac2 = R2(arr)/Q2(arr)
        
        c_lim1 = np.argmin(abs(0.5-frac1[np.isfinite(frac1)]))
        c_lim2 = np.argmin(abs(0.5-frac2[np.isfinite(frac2)]))
        c_lim1a = np.argmin(abs(0.9-frac1[np.isfinite(frac1)]))
        c_lim2a = np.argmin(abs(0.9-frac2[np.isfinite(frac2)]))
        
        
        print i,arr[c_lim1],arr[c_lim2],arr[c_lim1a],arr[c_lim2a]
        plt.plot(arr,frac1,color=cols[i])
        plt.plot(arr,frac2,color=cols[i],lw=2)
        if i == 0: 
            inmag1 = c_lim1
            inmag2 = c_lim2
            if write_comps == True:
                [comp_corr_file.write('%.2f '%bin) for bin in Rbins1]
                comp_corr_file.write('\n')
        else:    
            print 1./frac1[inmag1],1./frac2[inmag2]
            if write_comps == True:
                comp_data_file.write('# inner mag limit corresponds to 50 percent completeness \n# 275: %.3f \n# 336: %.3f \n'%(1./frac1[inmag1],1./frac2[inmag2]))

        if write_comps == True:
            comp_corr_file.write('F336W %i '%i)
            [comp_corr_file.write('%.8f '%frac) for frac in frac1]
            comp_corr_file.write('\nF275W %i '%i)
            [comp_corr_file.write('%.8f '%frac) for frac in frac2]
            comp_corr_file.write('\n')
        
        # for Table1 need comp50,98 index 0 to be inner ring.
        comp50mag1 = Rbins1[c_lim1]
        comp50mag2 = Rbins1[c_lim2]
        comp92mag1 = Rbins1[c_lim1a]
        comp92mag2 = Rbins1[c_lim2a]
        if write_comps == True:
            comp_data_file.write('%i %.1f %.1f %.1f %.1f\n'%(i,
                                                             comp50mag1,
                                                             comp50mag2,
                                                             comp92mag1,
                                                             comp92mag2))
        print '%i %.1f %.1f %.1f %.1f'%(i,comp50mag1,comp50mag2,comp92mag1,comp92mag2)
        
        Rinds1s.append(Rinds1)
        Rinds2s.append(Rinds2)
        Bins.append(Rbins1)
        Frac1s.append(frac1)
        Frac2s.append(frac2)
        Rhist1s.append(Rhist1)
        Rhist2s.append(Rhist2)
        Qhist1s.append(Qhist1)
        Qhist2s.append(Qhist2)
        
    ax1,ax2 = ast_completeness_plot(Bins,Frac1s,Frac2s,**kwargs)
    ast_diff_plot(ast,Rinds1s,Rinds2s,inds,**kwargs)
    ast_hist_plot(Bins,Rhist1s,Qhist1s,Rhist2s,Qhist2s,**kwargs)
    
    comp_corr_file.close()
    comp_data_file.close()

def main(**kwargs):
    ast_file = kwargs.get('ast_file')
    ast = load_ast_file(ast_file)
    from params import PHAT_DIR_EXTRA3
    contour = PHAT_DIR_EXTRA3 + 'contours_v3.reg'
    cras, cdecs = file_io.read_reg(contour)
    outverts = np.column_stack((cras[-1],cdecs[-1]))
    
    rec1 = recovered_asts(ast['mag1in'],ast['mag1out'])
    rec2 = recovered_asts(ast['mag2in'],ast['mag2out'])
    
    # quality cut on all asts
    qc1 = quality_cut(ast['sharp1'],
                      ast['snr1'], 
                      ast['crowd1'],
                      ast['round1'],
                      **kwargs)
    
    qc2 = quality_cut(ast['sharp2'],
                      ast['snr2'],
                      ast['crowd2'],
                      ast['round2'],
                      **kwargs)
    
    qcd = list(set(qc1) & set(qc2))
    inds = points_inside_ds9_polygon(contour,ast['ra'],ast['dec'])
    comp_corr_file = open(kwargs.get('comp_file'),'w') 
    comp_corr_file.write('# region (0=in) mags comp_frac \n')
    comp_data_file = open(kwargs.get('comp_data'),'w') 
    print 'comp correction file:',kwargs.get('comp_file')
    print 'comp data file:',comp_data_file
    single = kwargs.get('single')
    kwargs['ast_file_ext'] = '_'+kwargs.get('comp_file').replace('.dat','').split('_')[-1]
    if single == False:
        qc1 = qcd
        qc2 = qcd
        kwargs['ast_file_ext'] += '_dd'
        color_range = (-0.7,0.5)
        ast_color = ast['mag1in']-ast['mag2in']
        color_inds = np.nonzero( (ast_color>color_range[0]) &
                                 (ast_color<color_range[1]))[0]
        
    Rinds1s , Rinds2s , Bins , Frac1s , Frac2s = [], [], [], [], []
    Rhist1s , Rhist2s , Qhist1s , Qhist2s = [], [], [], []
    comp_data_file.write('# Region comp50mag1 comp50mag2 comp90mag1 comp90mag2\n')
    print '# Region comp50mag1 comp50mag2 comp90mag1 comp90mag2'
    for i in range(len(inds)):
        # in each annulus, passed qc, recovered
        if single == False:
            Rinds1 = list(set(rec1) & set(qc1) & set(inds[i]) &
                          set(color_inds))
            Rinds2 = list(set(rec2) & set(qc2) & set(inds[i]) &
                          set(color_inds))
            # REWRITING INDS OF I FOR THIS COLOR RANGE.
            inds[i] = list(set(inds[i]) & set(color_inds))
        else:
            Rinds1 = list(set(rec1) & set(qc1) & set(inds[i]))
            Rinds2 = list(set(rec2) & set(qc2) & set(inds[i]))
        
        hist_range = (18.34,28.)   
        # use hist_bin_optimization to get this value
        #spacing = 0.04 
        spacing = 0.2 
        Rhist1,Rbins1 = np.histogram(ast['mag1in'][Rinds1],bins=np.arange(hist_range[0],hist_range[1],spacing))
        Rhist2,Rbins2 = np.histogram(ast['mag2in'][Rinds2],bins=np.arange(hist_range[0],hist_range[1],spacing))
        Qhist1,Qbins1 = np.histogram(ast['mag1in'][inds[i]],bins=np.arange(hist_range[0],hist_range[1],spacing))
        Qhist2,Qbins2 = np.histogram(ast['mag2in'][inds[i]],bins=np.arange(hist_range[0],hist_range[1],spacing))
        
        # need to be float...
        frac1 = np.array(map(float,Qhist1))/np.array(map(float,Rhist1))
        frac2 = np.array(map(float,Qhist2))/np.array(map(float,Rhist2))
        
        c_lim1 = GenUtils.closest_match(0.5,1/frac1) 
        c_lim2 = GenUtils.closest_match(0.5,1/frac2)
        c_lim1a = GenUtils.closest_match(0.90,1/frac1)
        c_lim2a = GenUtils.closest_match(0.90,1/frac2)
        
        if i == 0: 
            inmag1 = c_lim1
            inmag2 = c_lim2
            [comp_corr_file.write('%.2f '%bin) for bin in Rbins1]
            comp_corr_file.write('\n')
        else:    
            print 1./frac1[inmag1],1./frac2[inmag2]
            comp_data_file.write('# inner mag limit corresponds to 50 percent completeness \n# 275: %.3f \n# 336: %.3f \n'%(1./frac1[inmag1],1./frac2[inmag2]))

        comp_corr_file.write('F336W %i '%i)
        [comp_corr_file.write('%.8f '%frac) for frac in frac1]
        comp_corr_file.write('\nF275W %i '%i)
        [comp_corr_file.write('%.8f '%frac) for frac in frac2]
        comp_corr_file.write('\n')
        
        # for Table1 need comp50,98 index 0 to be inner ring.
        comp50mag1 = Rbins1[c_lim1]
        comp50mag2 = Rbins1[c_lim2]
        comp92mag1 = Rbins1[c_lim1a]
        comp92mag2 = Rbins1[c_lim2a]
        
        comp_data_file.write('%i %.1f %.1f %.1f %.1f\n'%(i,
                                                         comp50mag1,
                                                         comp50mag2,
                                                         comp92mag1,
                                                         comp92mag2))
        print '%i %.1f %.1f %.1f %.1f'%(i,comp50mag1,comp50mag2,comp92mag1,comp92mag2)
        
        Rinds1s.append(Rinds1)
        Rinds2s.append(Rinds2)
        Bins.append(Rbins1)
        Frac1s.append(frac1)
        Frac2s.append(frac2)
        Rhist1s.append(Rhist1)
        Rhist2s.append(Rhist2)
        Qhist1s.append(Qhist1)
        Qhist2s.append(Qhist2)
        
    ax1,ax2 = ast_completeness_plot(Bins,Frac1s,Frac2s,**kwargs)
    ast_diff_plot(ast,Rinds1s,Rinds2s,inds,**kwargs)
    ast_hist_plot(Bins,Rhist1s,Qhist1s,Rhist2s,Qhist2s,**kwargs)
    
    comp_corr_file.close()
    comp_data_file.close()

if __name__ == "__main__":
    home = os.getcwd()
    default_file = sys.argv[1]
    #default_file = 'defaults.dat'
    kwargs = file_io.load_defaults(default_file)
    kwargs['cols'] = ['#007929', '#0D56A6', '#4c063e', '#FF2800', 'black']
    
    os.chdir(kwargs.get('ast_diag'))
    ast_file = kwargs.get('ast_file')
    if ast_file == None:
        ast = read_all_asts(verts=outverts)
        write_asts(ast,'B01_ASTs.dat',**kwargs)
        kwargs['ast_file'] = ast_file
    
    #main(**kwargs)
    #new_mag_cut(**kwargs)
    kwargs['two_panel'] = False
    kwargs['double_plot'] = False
    ast = load_ast_file(ast_file)
    ast_LF(ast,**kwargs)
    os.chdir(home)
