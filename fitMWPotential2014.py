import os, os.path
import numpy
from scipy import optimize, integrate
from optparse import OptionParser
from galpy import potential
from galpy.util import save_pickles, bovy_plot, bovy_conversion
from matplotlib import pyplot
import readTerminalData
_REFR0= 8.
_REFV0= 220.
def fitMWPotential2014(options,args):
    #First read the surface densities
    surffile= 'mwpot14data/bovyrix13kzdata.csv'
    if not surffile is None and os.path.exists(surffile):
        surf= numpy.loadtxt(surffile,delimiter=',')
        surfrs= surf[:,2]
        kzs= surf[:,6]
        kzerrs= surf[:,7]
    else:
        raise IOError("%s does not exist" % surffile)
    #Read the terminal velocity data if necessary
    cl_glon, cl_vterm, cl_corr= readTerminalData.readClemens(dsinl=options.termdsinl)
    mc_glon, mc_vterm, mc_corr= readTerminalData.readMcClureGriffiths(dsinl=options.termdsinl)
    termdata= (cl_glon,cl_vterm/_REFV0,cl_corr,
               mc_glon,mc_vterm/_REFV0,mc_corr)
    if options.twodisks:
        init_params= [0.25,0.25,0.45,
                      numpy.log(3.5/_REFR0),numpy.log(0.4/_REFR0),
                      numpy.log(2./_REFR0),numpy.log(0.4/_REFR0),
                      numpy.log(20./_REFR0),
                      0.,0.]
    else:
        init_params= [0.5,0.45,
                      numpy.log(2.5/_REFR0),numpy.log(0.4/_REFR0),
                      numpy.log(20./_REFR0),
                      0.,0.]
    funcargs= (surfrs,kzs,kzerrs,termdata,options)
    #Optimize likelihood
    if True:
        params= optimize.fmin_powell(like_func,init_params,
                                     args=funcargs)
    else:
        params= init_params
    #Round the parameters, such that they are simpler numbers
    params[0]= numpy.round(20.*params[0])/20.
    params[1]= numpy.round(20.*params[1])/20.
    if options.twodisks:
        params[2]= numpy.round(20.*params[2])/20.
    params[2+options.twodisks]= numpy.log(numpy.round(_REFR0*numpy.exp(params[2+options.twodisks])*10.)/10./_REFR0)
    params[3+options.twodisks]= numpy.log(numpy.round(_REFR0*numpy.exp(params[3+options.twodisks])*100.)/100./_REFR0)
    if options.twodisks:
        params[5]= numpy.log(numpy.round(_REFR0*numpy.exp(params[5])*10.)/10./_REFR0)
        params[6]= numpy.log(numpy.round(_REFR0*numpy.exp(params[6])*100.)/100./_REFR0)
    params[4+3*options.twodisks]= numpy.log(numpy.round(_REFR0*numpy.exp(params[4+3*options.twodisks]))/_REFR0)   
    print "Best-fit parameters:"
    print "Disk contribution to F_R(R_0)", params[0]
    print "Halo contribution to F_R(R_0)", params[1]
    print "Disk scale length / kpc:", numpy.exp(params[2])*_REFR0
    print "Disk scale height / pc:", numpy.exp(params[3])*1000*_REFR0
    print "Halo scale radius / kpc:", numpy.exp(params[4])*_REFR0
    save_pickles(args[0],params)
    #Make plots comparing the fit to the data
    if options.twodisks:
        pot= [potential.PowerSphericalPotentialwCutoff(normalize=1.-params[0]-params[1]-params[2],
                                                       alpha=1.8,rc=1.9/_REFR0),
              potential.MiyamotoNagaiPotential(normalize=params[0],
                                               a=numpy.exp(params[3]),
                                               b=numpy.exp(params[4])),
              potential.MiyamotoNagaiPotential(normalize=params[1],
                                               a=numpy.exp(params[5]),
                                               b=numpy.exp(params[6])),
              potential.NFWPotential(normalize=params[2],
                                     a=numpy.exp(params[7]))]
    else:
        pot= [potential.PowerSphericalPotentialwCutoff(normalize=1.-params[0]-params[1],
                                                       alpha=1.8,rc=1.9/_REFR0),
              potential.MiyamotoNagaiPotential(normalize=params[0],
                                               a=numpy.exp(params[2]),
                                               b=numpy.exp(params[3])),
              potential.NFWPotential(normalize=params[1],a=numpy.exp(params[4]))]
    plotRotcurve(pot,options.rotcurvename,options)
    fzrd= plotKz(pot,options.kzcurvename,surfrs,kzs,kzerrs)
    plotTerm(pot,options.termcurvename,termdata)
    plotPot(pot,options.potname)
    plotDens(pot,options.densname)
    writeTable(pot,params,options.tablename,options,fzrd)
    return None

def like_func(params,surfrs,kzs,kzerrs,
              termdata,options):
    #Check ranges
    if params[0] < 0. or params[0] > 1.: return numpy.finfo(numpy.dtype(numpy.float64)).max
    if params[1] < 0. or params[1] > 1.: return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.twodisks and (params[2] < 0. or params[2] > 1.): return numpy.finfo(numpy.dtype(numpy.float64)).max
    if (1.-params[0]-params[1]-options.twodisks*params[2]) < 0. or (1.-params[0]-params[1]-options.twodisks*params[2]) > 1.: return numpy.finfo(numpy.dtype(numpy.float64)).max
    if params[2+options.twodisks] < numpy.log(1./_REFR0) or params[2+options.twodisks] > numpy.log(8./_REFR0):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if params[3+options.twodisks] < numpy.log(0.05/_REFR0) or params[3+options.twodisks] > numpy.log(1./_REFR0):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.twodisks:
        if params[5] < numpy.log(1./_REFR0) or params[5] > numpy.log(8./_REFR0):
            return numpy.finfo(numpy.dtype(numpy.float64)).max
        if params[6] < numpy.log(0.05/_REFR0) or params[6] > numpy.log(1./_REFR0):
            return numpy.finfo(numpy.dtype(numpy.float64)).max
    #Setup potential
    if options.twodisks:
        pot= [potential.PowerSphericalPotentialwCutoff(normalize=1.-params[0]-params[1]-params[2],
                                                       alpha=1.8,rc=1.9/_REFR0),
              potential.MiyamotoNagaiPotential(normalize=params[0],
                                               a=numpy.exp(params[3]),
                                               b=numpy.exp(params[4])),
              potential.MiyamotoNagaiPotential(normalize=params[1],
                                               a=numpy.exp(params[5]),
                                               b=numpy.exp(params[6])),
              potential.NFWPotential(normalize=params[2],
                                     a=numpy.exp(params[7]))]
    else:
        pot= [potential.PowerSphericalPotentialwCutoff(normalize=1.-params[0]-params[1],
                                                       alpha=1.8,rc=1.9/_REFR0),
              potential.MiyamotoNagaiPotential(normalize=params[0],
                                               a=numpy.exp(params[2]),
                                               b=numpy.exp(params[3])),
              potential.NFWPotential(normalize=params[1],a=numpy.exp(params[4]))]
    #Calculate model surface density at surfrs
    modelkzs= numpy.empty_like(surfrs)
    for ii in range(len(surfrs)):
        modelkzs[ii]= -potential.evaluatezforces(surfrs[ii]/_REFR0,1.1/_REFR0,
                                                 pot)*bovy_conversion.force_in_2piGmsolpc2(_REFV0,_REFR0)
    out= 0.5*numpy.sum((kzs-modelkzs)**2./kzerrs**2.)
    #Add terminal velocities
    vrsun= params[5+3*options.twodisks]
    vtsun= params[6+3*options.twodisks]
    cl_glon, cl_vterm, cl_corr, mc_glon, mc_vterm, mc_corr= termdata
    #Calculate terminal velocities at data glon
    cl_vterm_model= numpy.zeros_like(cl_vterm)
    for ii in range(len(cl_glon)):
        cl_vterm_model[ii]= potential.vterm(pot,cl_glon[ii])
    cl_vterm_model+= vrsun*numpy.cos(cl_glon/180.*numpy.pi)\
        -vtsun*numpy.sin(cl_glon/180.*numpy.pi)
    mc_vterm_model= numpy.zeros_like(mc_vterm)
    for ii in range(len(mc_glon)):
        mc_vterm_model[ii]= potential.vterm(pot,mc_glon[ii])
    mc_vterm_model+= vrsun*numpy.cos(mc_glon/180.*numpy.pi)\
        -vtsun*numpy.sin(mc_glon/180.*numpy.pi)
    cl_dvterm= (cl_vterm-cl_vterm_model)/options.termsigma*_REFV0
    mc_dvterm= (mc_vterm-mc_vterm_model)/options.termsigma*_REFV0
    out+= 0.5*numpy.sum(cl_dvterm*numpy.dot(cl_corr,cl_dvterm))
    out+= 0.5*numpy.sum(mc_dvterm*numpy.dot(mc_corr,mc_dvterm))
    #Rotation curve constraint
    out-= logprior_dlnvcdlnr(potential.dvcircdR(pot,1.))
    #K dwarfs, Kz
    out+= 0.5*(-potential.evaluatezforces(1.,1.1/_REFR0,pot)*bovy_conversion.force_in_2piGmsolpc2(_REFV0,_REFR0)-67.)**2./36.
    #K dwarfs, visible
    out+= 0.5*(visible_dens(pot,options)-55.)**2./25.
    #Local density prior
    localdens= potential.evaluateDensities(1.,0.,pot)*bovy_conversion.dens_in_msolpc3(_REFV0,_REFR0)
    out+= 0.5*(localdens-0.102)**2./0.01**2.
    #Bulge velocity dispersion
    out+= 0.5*(bulge_dispersion(pot)-117.)**2./225.
    #Mass at 60 kpc
    out+= 0.5*(mass60(pot,options)-4.)**2./0.7**2.
    #Concentration prior?
    return out

def pdf_func(params,*args):
    return -like_func(params,*args)

def mass60(pot,options):
    """The mass at 60 kpc in 10^11 msolar"""
    tR= 60./_REFR0
    #For the MN potential, we just assume that all of its mass is enclosed within 60 kpc, even though this isn't technically correct (but it's WRONG)
    if options.twodisks:
        return (pot[0].mass(tR)+pot[3].mass(tR)+pot[1]._amp+pot[2]._amp)\
            *bovy_conversion.mass_in_1010msol(_REFV0,_REFR0)/10.
    else:
        return (pot[0].mass(tR)+pot[2].mass(tR)+pot[1]._amp)\
            *bovy_conversion.mass_in_1010msol(_REFV0,_REFR0)/10.

def bulge_dispersion(pot):
    """The expected dispersion in Baade's window, in km/s"""
    bar, baz= 0.0175, 0.068
    return numpy.sqrt(1./pot[0].dens(bar,baz)*integrate.quad(lambda x: -potential.evaluatezforces(bar,x,pot)*pot[0].dens(bar,x),baz,numpy.inf)[0])*_REFV0

def visible_dens(pot,options,r=1.):
    """The visible surface density at 8 kpc from the center"""
    if options.twodisks:
        return 2.*integrate.quad((lambda zz: potential.evaluateDensities(r,zz,pot[1:3])),0.,2.)[0]*bovy_conversion.surfdens_in_msolpc2(_REFV0,_REFR0)
    else:
        return 2.*integrate.quad((lambda zz: potential.evaluateDensities(r,zz,pot[1])),0.,2.)[0]*bovy_conversion.surfdens_in_msolpc2(_REFV0,_REFR0)

def logprior_dlnvcdlnr(dlnvcdlnr):
    sb= 0.04
    if dlnvcdlnr > sb or dlnvcdlnr < -0.5:
        return -numpy.finfo(numpy.dtype(numpy.float64)).max
    return numpy.log((sb-dlnvcdlnr)/sb)-(sb-dlnvcdlnr)/sb

#########################################PLOTS#################################
def plotRotcurve(pot,plotfilename,options):
    potential.plotRotcurve(pot,xrange=[0.,4.],color='k',lw=2.,yrange=[0.,1.4])
    #Constituents
    line1= potential.plotRotcurve(pot[0],overplot=True,color='k',ls='-.',lw=2.)
    if options.twodisks:
        line2= potential.plotRotcurve(pot[1:3],overplot=True,color='k',ls='--',lw=2.)
        line3= potential.plotRotcurve(pot[3],overplot=True,color='k',ls=':',lw=2.)
    else:
        line2= potential.plotRotcurve(pot[1],overplot=True,color='k',ls='--',lw=2.)
        line3= potential.plotRotcurve(pot[2],overplot=True,color='k',ls=':',lw=2.)
    #Add legend
    pyplot.legend((line1[0],line2[0],line3[0]),
                  (r'$\mathrm{Bulge}$',
                   r'$\mathrm{Disk}$',
                   r'$\mathrm{Halo}$'),
                  loc='upper right',#bbox_to_anchor=(.91,.375),
                  numpoints=8,
                  prop={'size':16},
                  frameon=False)
    bovy_plot.bovy_end_print(plotfilename)
    return None

def plotKz(pot,plotfilename,surfrs,kzs,kzerrs):
    krs= numpy.linspace(4./_REFR0,10./_REFR0,1001)
    modelkz= numpy.array([-potential.evaluatezforces(kr,1.1/_REFR0,pot)\
                               *bovy_conversion.force_in_2piGmsolpc2(_REFV0,_REFR0) for kr in krs])
    bovy_plot.bovy_print(fig_height=3.5)
    bovy_plot.bovy_plot(krs*_REFR0,modelkz,'-',color='0.6',lw=2.,
                        xlabel=r'$R\ (\mathrm{kpc})$',
                        ylabel=r'$F_{Z}(R,|Z| = 1.1\,\mathrm{kpc})\ (2\pi G\,M_\odot\,\mathrm{pc}^{-2})$',
                        semilogy=True,
                        yrange=[10.,1000.],
                        xrange=[4.,10.],
                        zorder=0)
    pyplot.errorbar(surfrs,
                    kzs,
                    yerr=kzerrs,
                    marker='o',
                    elinewidth=1.,capsize=3,zorder=1,
                    color='k',linestyle='none')  
    pyplot.errorbar([_REFR0],[69.],yerr=[6.],marker='d',ms=10.,
                    elinewidth=1.,capsize=3,zorder=10,
                    color='0.4',linestyle='none')
    bovy_plot.bovy_end_print(plotfilename)
    #Do an exponential fit to the model Kz and return the scale length
    indx= krs < 9./_REFR0
    p= numpy.polyfit(krs[indx],numpy.log(modelkz[indx]),1)
    return -1./p[0]

def plotTerm(pot,plotfilename,termdata):
    mglons= numpy.linspace(-90.,-20.,1001)
    pglons= numpy.linspace(20.,90.,1001)
    mterms= numpy.array([potential.vterm(pot,mgl)*_REFV0 for mgl in mglons])
    pterms= numpy.array([potential.vterm(pot,pgl)*_REFV0 for pgl in pglons])
    bovy_plot.bovy_print(fig_height=3.5)
    bovy_plot.bovy_plot(mglons,mterms,'-',color='0.6',lw=2.,zorder=0,
                        xlabel=r'$\mathrm{Galactic\ longitude\, (deg)}$',
                        ylabel=r'$\mathrm{Terminal\ velocity}\, (\mathrm{km\,s}^{-1})$',
                        xrange=[-100.,100.],
                        yrange=[-150.,150.])
    bovy_plot.bovy_plot(pglons,pterms,'-',color='0.6',lw=2.,zorder=0,
                        overplot=True)
    cl_glon,cl_vterm,cl_corr,mc_glon,mc_vterm,mc_corr= termdata
    bovy_plot.bovy_plot(cl_glon,cl_vterm*_REFV0,'ko',overplot=True)
    bovy_plot.bovy_plot(mc_glon-360.,mc_vterm*_REFV0,'ko',overplot=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None

def plotPot(pot,plotfilename):
    potential.plotPotentials(pot,rmin=0.,rmax=1.5,nrs=201,
                             zmin=-0.5,zmax=0.5,nzs=201,ncontours=21,
                             justcontours=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None
    
def plotDens(pot,plotfilename):
    potential.plotDensities(pot,rmin=0.01,rmax=1.5,nrs=201,
                            zmin=-0.5,zmax=0.5,nzs=201,ncontours=21,
                            log=True,justcontours=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None
    
def writeTable(pot,params,tablename,options,fzrd):
    outfile= open(tablename,'w')
    delimiter= ' & '
    #First list the explicit parameters
    printline= '$R_0\,(\kpc)$%s$8$%sfixed\\\\\n' % (delimiter,delimiter)
    outfile.write(printline)
    printline= '$v_c(R_0)\,(\kms)$%s$220$%sfixed\\\\\n' % (delimiter,delimiter)
    outfile.write(printline)
    if options.twodisks:
        printline= '$f_{b}$%s$%.2f$%s\ldots\\\\\n' % (delimiter,1.-params[0]-params[1]-params[2],delimiter)        
        outfile.write(printline)
        printline= '$f_{d,1}$%s$%.2f$%s\ldots\\\\\n' % (delimiter,params[0],delimiter)
        outfile.write(printline)
        printline= '$f_{d,2}$%s$%.2f$%s\ldots\\\\\n' % (delimiter,params[1],delimiter)        
        outfile.write(printline)
        printline= '$f_{h}$%s$%.2f$%s\ldots\\\\\n' % (delimiter,params[2],delimiter)        
        outfile.write(printline)
    else:
        printline= '$f_{b}$%s$%.2f$%s\ldots\\\\\n' % (delimiter,1.-params[0]-params[1],delimiter)        
        outfile.write(printline)
        printline= '$f_{d}$%s%.2f%s\ldots\\\\\n' % (delimiter,params[0],delimiter)
        outfile.write(printline)
        printline= '$f_{h}$%s$%.2f$%s\ldots\\\\\n' % (delimiter,params[1],delimiter)        
        outfile.write(printline)
    #Bulge power-law and rc
    printline= '$\mathrm{Bulge\ power}$%s$-1.8$%sfixed\\\\\n' % (delimiter,delimiter)        
    outfile.write(printline)
    printline= '$\mathrm{Bulge\ cut}\,(\kpc)$%s$1.9$%sfixed\\\\\n' % (delimiter,delimiter)        
    outfile.write(printline)
    #Disk scale length and height
    printline= '$a\,(\kpc)$%s$%.1f$%s\ldots\\\\\n' % (delimiter,numpy.exp(params[2+options.twodisks])*_REFR0,delimiter)        
    outfile.write(printline)
    printline= '$b\,(\pc)$%s$%.0f$%s\ldots\\\\\n' % (delimiter,1000.*numpy.exp(params[3+options.twodisks])*_REFR0,delimiter)        
    outfile.write(printline)
    printline= '$\mathrm{Halo}\ r_s\,(\kpc)$%s$%.0f$%s\ldots\\\\\n' % (delimiter,numpy.exp(params[4+options.twodisks])*_REFR0,delimiter)        
    outfile.write(printline)
    outfile.write('%s%s\\\\\n' % (delimiter,delimiter))
    #Now the constraints
    printline= '$\sigma_b\,(\kms)$%s$%.0f$%s$117\pm15$\\\\\n' % (delimiter,bulge_dispersion(pot),delimiter)
    outfile.write(printline)
    printline= '$F_Z(R_0,1.1\kpc)\,(2\pi G\,M_\odot\pc^{-2})$%s$%.0f$%s$67\pm6$\\\\\n' % (delimiter,-potential.evaluatezforces(1.,1.1/_REFR0,pot)*bovy_conversion.force_in_2piGmsolpc2(_REFV0,_REFR0),delimiter)
    outfile.write(printline)
    printline= '$\Sigma_{\mathrm{vis}}(R_0)\,(M_\odot\pc^{-2})$%s$%.0f$%s$55\pm5$\\\\\n' % (delimiter,visible_dens(pot,options),delimiter)
    outfile.write(printline)
    printline= '$F_Z\ \mathrm{scale\ length}\,(\kpc)$%s%.1f%s$2.7\pm0.1$\\\\\n' % (delimiter,fzrd*_REFR0,delimiter)
    outfile.write(printline)
    printline= '$\\rho(R_0,z=0)\,(M_\odot\pc^{-3})$%s$%.2f$%s$0.10\pm0.01$\\\\\n' % (delimiter,potential.evaluateDensities(1.,0.,pot)*bovy_conversion.dens_in_msolpc3(_REFV0,_REFR0),delimiter)
    outfile.write(printline)
    printline= '$(\dd \ln v_c/\dd \ln R)|_{R_0}$%s$%.2f$%s$-0.2\ \mathrm{to}\ 0$\\\\\n' % (delimiter,potential.dvcircdR(pot,1.),delimiter)
    outfile.write(printline)
    printline= '$M(r<60\kpc)\,(10^{11}\,M_\odot)$%s$%.1f$%s$4.0\pm0.7$\\\\\n' % (delimiter,mass60(pot,options),delimiter)
    outfile.write(printline)
    outfile.write('%s%s\\\\\n' % (delimiter,delimiter))
    #Now some derived properties
    printline= '$M_b\,(10^{10}\,M_\odot)$%s$%.1f$%s\ldots\\\\\n' % (delimiter,pot[0].mass(10.)*bovy_conversion.mass_in_1010msol(_REFV0,_REFR0),delimiter)
    outfile.write(printline)
    if options.twodisks:
        printline= '$M_d\,(10^{10}\,M_\odot)$%s$%.1f$%s\ldots\\\\\n' % (delimiter,(pot[1]._amp+pot[2]._amp)*bovy_conversion.mass_in_1010msol(_REFV0,_REFR0),delimiter)
    else:
        printline= '$M_d\,(10^{10}\,M_\odot)$%s$%.1f$%s\ldots\\\\\n' % (delimiter,pot[1]._amp*bovy_conversion.mass_in_1010msol(_REFV0,_REFR0),delimiter)
    outfile.write(printline)
    krs= numpy.linspace(4./_REFR0,9./_REFR0,101)
    disksurf= numpy.array([visible_dens(pot,options,r=kr) for kr in krs])
    p= numpy.polyfit(krs,numpy.log(disksurf),1)
    rd= -1./p[0]*_REFR0
    printline= '$R_{d}\,(\kpc)$%s$%.1f$%s\ldots\\\\\n' % (delimiter,rd,delimiter)
    outfile.write(printline)
    rvir= pot[2+options.twodisks].rvir(_REFV0,_REFR0,wrtcrit=True,overdens=96.7)
    printline= '$\\rho_{\mathrm{DM}}(R_0)\,(M_\odot\pc^{-3})$%s$%.3f$%s\\\\\n' % (delimiter,potential.evaluateDensities(1.,0.,pot[2+options.twodisks])*bovy_conversion.dens_in_msolpc3(_REFV0,_REFR0),delimiter)
    outfile.write(printline)
    printline= '$M_{\mathrm{vir}}\,(10^{12}\,M_\odot)$%s$%.1f$%s\ldots\\\\\n' % (delimiter,pot[2+options.twodisks].mass(rvir)*bovy_conversion.mass_in_1010msol(_REFV0,_REFR0)/100.,delimiter)
    outfile.write(printline)
    printline= '$r_{\mathrm{vir}}\,(\kpc)$%s$%.0f$%s\ldots\\\\\n' % (delimiter,rvir*_REFR0,delimiter)
    outfile.write(printline)
    printline= '$\mathrm{Concentration}$%s$%.1f$%s\ldots\\\\\n' % (delimiter,rvir/pot[2+options.twodisks].a,delimiter)
    outfile.write(printline)
    printline= '$v_{\mathrm{esc}}(R_0)\,(\kms)$%s$%.0f$%s\ldots\n' % (delimiter,potential.vesc(pot,1.)*_REFV0,delimiter)
    outfile.write(printline)
    #printline= '$r_{\mathrm{vir}}\,(\kpc)$%s$%.0f$%s\ldots\\\\\n' % (delimiter,delimiter)
    #outfile.write(printline)
    
    outfile.write('\\enddata\n')
    pass

def get_options():
    usage = "usage: %prog [options] <savefile>\n\nsavefile= name of the file that the fits/samples will be saved to"
    parser = OptionParser(usage=usage)
    #Data options
    parser.add_option("--termdsinl",dest='termdsinl',default=0.125,type='float',
                      help="Correlation length for terminal velocity residuals")
    parser.add_option("--termsigma",dest='termsigma',default=7.,type='float',
                      help="sigma for terminal velocity residuals")
    #Fit options
    parser.add_option("--twodisks",action="store_true", 
                      dest="twodisks",
                      default=False,
                      help="If set, fit two disk components")
    #Output options
    parser.add_option("--init",dest='initfile',default=None,
                      help="Name of the file that has the best-fits")
    parser.add_option("--rotcurve",dest='rotcurvename',default=None,
                      help="Name of the file that holds the plot of the rotation curve")
    parser.add_option("--kzcurve",dest='kzcurvename',default=None,
                      help="Name of the file that holds the plot of the Kz curve")
    parser.add_option("--termcurve",dest='termcurvename',default=None,
                      help="Name of the file that holds the plot of the Kz curve")
    parser.add_option("--potname",dest='potname',default=None,
                      help="Name of the file that holds the plot of the potential")
    parser.add_option("--densname",dest='densname',default=None,
                      help="Name of the file that holds the plot of the density")
    parser.add_option("--tablename",dest='tablename',default=None,
                      help="Name of the file that holds the table with properties")
    #seed
    parser.add_option("--seed",dest='seed',default=1,type='int',
                      help="seed for random number generator")
    return parser

if __name__ == '__main__':
    parser= get_options()
    options,args= parser.parse_args()
    numpy.random.seed(options.seed)
    fitMWPotential2014(options,args)
