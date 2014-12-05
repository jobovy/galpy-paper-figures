import sys
import time
import numpy
from galpy.potential import MWPotential2014
from galpy.potential import evaluatePotentials as evalPot
from galpy.actionAngle import actionAngleStaeckel
from galpy.orbit import Orbit
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import matplotlib.ticker as ticker
def plot_energy_conservation(plotfilename):
    #First orbit
    E, Lz= -1.25, 0.6
    aAS= actionAngleStaeckel(pot=MWPotential2014,c=True,delta=0.5)
    o= Orbit([0.8,0.3,Lz/0.8,0.,numpy.sqrt(2.*(E-evalPot(0.8,0.,MWPotential2014)-(Lz/0.8)**2./2.-0.3**2./2.)),0.])
    orbt= 2.*numpy.pi/aAS.actionsFreqs(o)[4]
    norb= 2000.
    nt= 20001
    ts= numpy.linspace(0.,norb*orbt,nt)
    start= time.time()
    o.integrate(ts,MWPotential2014,method='dopr54_c')
    print 'dopr54_c', time.time()-start
    Es= o.E(ts)
    dEs= numpy.fabs((Es-Es[0])/Es[0])
    breakt= 80.
    pts= list(ts[ts < breakt])
    pts.extend(list((ts[ts >= breakt])[::10]))
    pts= numpy.array(pts)
    pdEs= list(dEs[ts < breakt])
    pdEs.extend(list((dEs[ts >= breakt])[::10]))
    pdEs= numpy.array(pdEs)
    print 'dopr54_c', numpy.mean(dEs)
    bovy_plot.bovy_print(fig_width=3.25,fig_height=9.)
    pyplot.subplot(8,1,8)
    bovy_plot.bovy_plot(pts/orbt,pdEs,color='k',
                        loglog=True,gcf=True,
                        xrange=[0.5,2000.],
                        yrange=[10.**-12.,1.],
                        xlabel=r'$\mathrm{Number\ of\ orbital\ periods}$')
    bovy_plot.bovy_text(r'$\texttt{dopr54\_c}$',
                        top_left=True,size=14.)
    ax= pyplot.gca()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter(r'$%0.f$'))
    ax.yaxis.set_ticks([10.**-12.,10.**-8.,10.**-4.,1.])
    nullfmt   = NullFormatter()         # no labels
    other_integrators= ['odeint','leapfrog','leapfrog_c',
                        'symplec4_c','symplec6_c',
                        'rk4_c','rk6_c']
    for ii,integrator in enumerate(other_integrators):
        start= time.time()
        o.integrate(ts,MWPotential2014,method=integrator)
        print integrator, time.time()-start
        Es= o.E(ts)
        dEs= numpy.fabs((Es-Es[0])/Es[0])
        pdEs= list(dEs[ts < breakt])
        pdEs.extend(list((dEs[ts >= breakt])[::10]))
        pdEs= numpy.array(pdEs)
        print integrator, numpy.mean(dEs)
        pyplot.subplot(8,1,ii+1)
        if ii == 3: ylabel= r'$\left|\Delta E/ E\right|$'
        else: ylabel= None
        bovy_plot.bovy_plot(pts/orbt,pdEs,color='k',
                            loglog=True,gcf=True,
                            xrange=[0.5,2000.],
                            yrange=[10.**-12.,1.],
                            ylabel=ylabel)
        thisax= pyplot.gca()
        thisax.xaxis.set_major_formatter(nullfmt)
        thisax.yaxis.set_ticks([10.**-12.,10.**-8.,10.**-4.,1.])
        bovy_plot.bovy_text(r'$\texttt{%s}$' % (integrator.replace('_','\_')),
                            top_left=True,size=14.)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_energy_conservation(sys.argv[1])
