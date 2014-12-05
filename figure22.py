import sys
import numpy
from galpy.df import dehnendf
from galpy.util import bovy_plot
from matplotlib import pyplot, cm
def plot_dfcorrections(plotfilename):
    niters= [1,2,3,4,5,10,15,20,25]
    bovy_plot.bovy_print(fig_height=7.,fig_width=8.)
    ii= 0
    # Load DF
    pyplot.subplot(2,1,1)
    dfc= dehnendf(beta=0.,correct=True,niter=niters[ii])
    bovy_plot.bovy_plot(dfc._corr._rs,
                        numpy.log(dfc._corr._corrections[:,0]),
                        '-',gcf=True,color=cm.jet(1.),lw=2.,zorder=1,
                        xrange=[0.,5.],
                        yrange=[-0.25,0.25],
                        ylabel=r'$\ln \Sigma_{\mathrm{out}}(R)-\ln\Sigma_{\mathrm{DF}}(R)$')
    linthresh= 0.0001
    pyplot.yscale('symlog',linthreshy=linthresh)
    for ii,niter in enumerate(niters[1:]):
        dfcn= dehnendf(beta=0.,correct=True,niter=niter)
        dfcp= dehnendf(beta=0.,correct=True,niter=niter-1)
        bovy_plot.bovy_plot(dfc._corr._rs,
                            numpy.log(dfcn._corr._corrections[:,0])-numpy.log(dfcp._corr._corrections[:,0]),
                            '-',overplot=True,
                            color=cm.jet(1.-(ii+1)/float(len(niters))),lw=2.,
                            zorder=ii+2)
    pyplot.fill_between(numpy.linspace(0.,5.,2.),
                        -linthresh*numpy.ones(2),
                         linthresh*numpy.ones(2),color='0.9',
                         zorder=0)
    bovy_plot.bovy_text(4.,-0.00008,r'$\mathrm{linear\ scale}$',
                        backgroundcolor='w',size=16.)
    pyplot.subplot(2,1,2)
    bovy_plot.bovy_plot(dfc._corr._rs,
                        0.5*numpy.log(dfc._corr._corrections[:,1]),
                        '-',gcf=True,color=cm.jet(1.),lw=2.,zorder=1,
                        xrange=[0.,5.],
                        yrange=[-0.25,0.25],
                        xlabel=r'$R/R_0$',
                        ylabel=r'$\ln \sigma_{R,\mathrm{out}}(R)-\ln\sigma_{R,\mathrm{DF}}(R)$')
    pyplot.yscale('symlog',linthreshy=linthresh)
    for ii,niter in enumerate(niters[1:]):
        dfcn= dehnendf(beta=0.,correct=True,niter=niter)
        dfcp= dehnendf(beta=0.,correct=True,niter=niter-1)
        bovy_plot.bovy_plot(dfc._corr._rs,
                            numpy.log(dfcn._corr._corrections[:,1])-numpy.log(dfcp._corr._corrections[:,1]),
                            '-',overplot=True,
                            color=cm.jet(1.-(ii+1)/float(len(niters))),lw=2.,
                            zorder=ii+2)
    pyplot.fill_between(numpy.linspace(0.,5.,2.),
                        -linthresh*numpy.ones(2),
                         linthresh*numpy.ones(2),color='0.9',
                         zorder=0)
    bovy_plot.bovy_text(4.,-0.00008,r'$\mathrm{linear\ scale}$',
                        backgroundcolor='w',size=16.)
    pyplot.tight_layout()
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_dfcorrections(sys.argv[1])
