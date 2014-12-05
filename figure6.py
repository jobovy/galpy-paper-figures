import sys
from galpy import potential
from galpy.util import bovy_plot
from matplotlib import pyplot
def plot_rotation_curves(plotfilename):
    #Setup all potentials
    dp= potential.DoubleExponentialDiskPotential(normalize=1.,
                                                 hr=0.4,hz=0.04)
    jp= potential.JaffePotential(normalize=1.,a=1.)
    hp= potential.HernquistPotential(normalize=1.,a=1.)
    kp= potential.KeplerPotential(normalize=1.)
    lp= potential.LogarithmicHaloPotential(normalize=1.)
    ip= potential.IsochronePotential(normalize=1.,b=1.)
    mp= potential.MiyamotoNagaiPotential(normalize=1.,a=0.4,b=0.04)
    np= potential.NFWPotential(normalize=1.,a=1.)
    pp= potential.PowerSphericalPotentialwCutoff(normalize=1.,alpha=1.,rc=1.)
    rp= potential.RazorThinExponentialDiskPotential(normalize=1.,hr=0.4)
    bp= potential.BurkertPotential(normalize=1.,a=1.)
    #Plot rotation curves
    Rrange= [0.001,5.]
    grid= 201
    bovy_plot.bovy_print(fig_width=6.)
    dpline= dp.plotRotcurve(Rrange=Rrange,grid=grid,color='r',ls='-',
                            yrange=[0.,1.9],xrange=[0.,5.],lw=2.)
    rpline= rp.plotRotcurve(Rrange=Rrange,grid=grid,ls=':',
                            color='r',lw=2.,
                            overplot=True)
    mpline= mp.plotRotcurve(Rrange=Rrange,grid=grid,ls='--',
                            color='r',lw=2.,
                            overplot=True)
    hpline= hp.plotRotcurve(Rrange=Rrange,grid=grid,ls='-',
                            color='y',lw=2.,
                            overplot=True)
    jpline= jp.plotRotcurve(Rrange=Rrange,grid=grid,ls=':',
                            color='y',lw=2.,
                            overplot=True)
    npline= np.plotRotcurve(Rrange=Rrange,grid=grid,ls='--',
                            color='y',lw=2.,
                            overplot=True)
    bpline= bp.plotRotcurve(Rrange=Rrange,grid=grid,ls='-.',
                            color='y',lw=2.,
                            overplot=True)
    lpline= lp.plotRotcurve(Rrange=Rrange,grid=grid,ls='-',
                            color='k',lw=2.,
                            overplot=True)
    kpline= kp.plotRotcurve(Rrange=Rrange,grid=grid,ls=':',
                            color='k',lw=2.,
                            overplot=True)
    ipline= ip.plotRotcurve(Rrange=Rrange,grid=grid,ls='--',
                            color='k',lw=2.,
                            overplot=True)
    ppline= pp.plotRotcurve(Rrange=Rrange,grid=grid,ls='-.',
                            color='k',lw=2.,
                            overplot=True)
    #Add legend
    legend1= pyplot.legend((dpline[0],rpline[0],mpline[0]),
                           (r'$\mathrm{Double\!-\!Exp}: h_R = 0.4, h_z = 0.04$',
                            r'$\mathrm{Razor\!-\!thin\ Exp}: h_R = 0.4$',
                            r'$\mathrm{Miya\!-\!Nagai}: h_R = 0.4, h_z = 0.04$'),
                           loc='lower right',#bbox_to_anchor=(.91,.375),
                           numpoints=2,
                           prop={'size':14},
                  frameon=False)
    legend2= pyplot.legend((hpline[0],jpline[0],npline[0],bpline[0]),
                           (r'$\mathrm{Hernquist}$',
                            r'$\mathrm{Jaffe}$',
                            r'$\mathrm{NFW}$',
                            r'$\mathrm{Burkert}$'),
                           loc='upper right',#bbox_to_anchor=(.91,.375),
                           numpoints=2,
                           prop={'size':14},
                  frameon=False)
    legend3= pyplot.legend((lpline[0],ipline[0],ppline[0],kpline[0]),
                           (r'$\mathrm{Logarithmic}$',
                            r'$\mathrm{Isochrone}$',
                            r'$\rho \propto r^{-1}\,e^{-r^2}$',
                            r'$\mathrm{Kepler}$'),
                           bbox_to_anchor=(.50,.985),
                           numpoints=2,
                           prop={'size':14},
                  frameon=False)
    pyplot.gca().add_artist(legend1)
    pyplot.gca().add_artist(legend2)
    pyplot.gca().add_artist(legend3)
    bovy_plot.bovy_end_print(plotfilename)

if __name__ == '__main__':
    plot_rotation_curves(sys.argv[1])
