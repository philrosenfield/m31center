import pyfits
import Astronomy
import graphics
import matplotlib.nxutils as nxutils
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter, MultipleLocator
nullfmt   = NullFormatter() # no labels

# Check for Disk contamination.
def disk_cont(phot,rad_inds,cut_inds,**kwargs):
    ftsize = kwargs.get('ftsize')
    disk_fit = kwargs.get('disk_fit')
    mag1_cut = kwargs.get('mag1_cut')
    mag2_cut = kwargs.get('mag2_cut')
    # exculding innermost region
    areas = kwargs.get('areas')[1:]
    distances = kwargs.get('distances')[1:]
    
    d1 = pyfits.open(disk_fit)[1].data
    dmag1 = d1.field('mag1')
    dmag2 = d1.field('mag2')
    dqc_inds = np.nonzero( (dmag1 < mag1_cut) & (dmag2 < mag2_cut))
    dmag = dmag2[dqc_inds]
    dcolor = dmag1[dqc_inds]-dmag
    
    mag_bot = 23.
    mag_top = 22.
    # outline of MS between 22-23 - done roughly by eye/plot
    col_bl = 0.013
    col_tl = -0.13
    col_br = col_bl+0.15
    col_tr = col_tl+0.15

    test_col = [col_tl,col_tr,col_tr,col_tl,col_tl]
    test_mag = [mag_top,mag_top,mag_bot,mag_bot,mag_top]
    test_verts = np.column_stack((test_col,test_mag))
    # outline of NOT MS between 22-23
    col_r = -0.35
    col_l = -0.50
    control_col = [col_l,col_r,col_r,col_l,col_l]
    control_mag = test_mag # just to be clear, same mag range...
    control_verts = np.column_stack((control_col,control_mag))
    axs = graphics.setup_five_panel_plot()
    Ncontrol,Ntest = [],[]
    control,test = [],[]
    
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[i]) & set(cut_inds))
        mag1 = phot['mag1'][inds]
        mag2 = phot['mag2'][inds]
        color = mag1-mag2
        colmag = np.column_stack((color,mag2))
        control=np.nonzero(nxutils.points_inside_poly(colmag,control_verts))[0]
        test = np.nonzero(nxutils.points_inside_poly(colmag,test_verts))[0]
        Ncontrol.append(float(len(control)))
        Ntest.append(float(len(test)))
        ax = axs[i]
        ax.plot(dcolor,dmag,',',color='grey',mec='grey')
        ax.plot(color[control],mag2[control],'.',color='black')
        ax.plot(control_verts[:,0],control_verts[:,1],lw=2,color='blue')
        ax.plot(test_verts[:,0],test_verts[:,1],lw=2,color='white')
        ax.plot(color[test],mag2[test],'.',color='red')
    
    Ncfit = np.polyfit(distances,Ncontrol,1)
    Ntfit = np.polyfit(distances,Ntest,1)
    plt.savefig('disk_cont.png')
    print 'wrote disk_cont.png'
    plt.close()    

    rel_errorNt = np.sqrt(Ntest)/Ntest
    rel_errorNc = np.sqrt(Ncontrol)/Ncontrol
    pct_err = [rt**2+rc**2+2*rt*rc for rt,rc in zip(rel_errorNt,rel_errorNc)]
    ratio = np.array(Ntest)/np.array(Ncontrol)
    ax1,ax2,ax3 = graphics.set_up_three_panel_plot()
    ax3.errorbar(distances,ratio,yerr=ratio*pct_err,color='black',lw=2)
    ax2.plot(dcolor,dmag,',',color='grey')
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[i]) & set(cut_inds))
        mag1 = phot['mag1'][inds]
        mag2 = phot['mag2'][inds]
        color = mag1-mag2
        ax1.plot(color,mag2,',',color='black')

    for ax in [ax1,ax2]:
        ax.plot(control_verts[:,0],control_verts[:,1],'--',lw=2,color='red')
        ax.plot(test_verts[:,0],test_verts[:,1],lw=2,color='darkblue')
        ax.set_ylim(24.1,19)
        ax.set_xlim(-0.75,2)
        ax.set_xlabel(r'$F275W-F336W$',fontsize=ftsize)

    ax2.yaxis.set_major_formatter(nullfmt)
    ax1.set_ylabel(r'$F336W$',fontsize=ftsize)
    ax3.set_xlabel(r'$\rm{Radius\ (pc)}$',fontsize=ftsize)
    ax3.set_ylabel(r'$\rm{Number\ Ratio Bluer/Redder}$',fontsize=ftsize)
    plt.savefig('nratio_cont.png')
    print 'wrote nratio_cont.png'
    plt.close()

    plt.errorbar(distances,Ntest/areas,yerr=np.sqrt(np.array(Ntest/areas)),
                 color='red',lw=2,label='Possible cont')
    plt.errorbar(distances,Ncontrol/areas,np.sqrt(np.array(Ncontrol/areas)),
                 color='black',lw=2,label='No cont')
    plt.xlabel('Region pc')
    plt.ylabel('Number in box/arcsec$^2$')
    plt.legend(loc=2)
    plt.savefig('density_cont.png')
    print 'wrote density_cont.png'
    plt.close()

    plt.figure()
    plt.plot(distances,Ntest,color='red',lw=2,label='Possible cont')
    plt.plot(distances,Ncontrol,color='black',lw=2,label='No cont')
    plt.xlabel('Region Out -- to -- In')
    plt.ylabel('Number in box')
    plt.legend(loc=0)
    plt.savefig('number_cont.png')
    print 'wrote number_cont.png'
    
    return

def red_star(phot,rad_inds,cut_inds,**kwargs):
    axs = graphics.setup_five_panel_plot()
    flux1,flux2 = [],[]
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[i]) & set(cut_inds))
        mag1 = phot['mag1'][inds]
        mag2 = phot['mag2'][inds]
        color = mag1-mag2
        red = np.nonzero((color > 0.8) & (mag2 < 19))[0]
        flux2.append(
            np.sum(Astronomy.vegamag2flux(mag2[red],Filter2,'wfc3uvis')))
        flux1.append(
            np.sum(Astronomy.vegamag2flux(mag1[red],Filter1,'wfc3uvis')))
        ax = axs[i]
        ax.plot(color,mag2,',',color='black')
        ax.plot(mag1[red]-mag2[red],mag2[red],',',color='red')
    return flux1,flux2
