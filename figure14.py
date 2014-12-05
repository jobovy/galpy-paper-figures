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
def plot_liouville(plotfilename):
    E, Lz= -1.25, 0.6
    aAS= actionAngleStaeckel(pot=MWPotential2014,c=True,delta=0.5)
    #planarOrbit w/ the same energy
    o= Orbit([0.8,numpy.sqrt(2.*(E-evalPot(0.8,0.,MWPotential2014)-(Lz/0.8)**2./2.)),Lz/0.8,0.])
    #For the orbital period
    fo= Orbit([0.8,numpy.sqrt(2.*(E-evalPot(0.8,0.,MWPotential2014)-(Lz/0.8)**2./2.)),Lz/0.8,0.,0.001,0.])
    orbt= 2.*numpy.pi/aAS.actionsFreqs(fo)[4]
    norb= 200.
    nt= 20001
    ts= numpy.linspace(0.,norb*orbt,nt)
    start= time.time()
    integrator= 'dopr54_c'
    o.integrate_dxdv([1.,0.,0.,0.],ts,MWPotential2014,
                     method=integrator,
                     rectIn=True,rectOut=True)
    dx= o.getOrbit_dxdv()[:,:]
    o.integrate_dxdv([0.,1.,0.,0.],ts,MWPotential2014,
                     method=integrator,
                     rectIn=True,rectOut=True)
    dy= o.getOrbit_dxdv()[:,:]
    o.integrate_dxdv([0.,0.,1.,0.],ts,MWPotential2014,
                     method=integrator,
                     rectIn=True,rectOut=True)
    dvx= o.getOrbit_dxdv()[:,:]
    o.integrate_dxdv([0.,0.,0.,1.],ts,MWPotential2014,
                     method=integrator,
                     rectIn=True,rectOut=True)
    dvy= o.getOrbit_dxdv()[:,:]
    print integrator, time.time()-start
    #Calculate Jacobian
    jacs= numpy.array([numpy.linalg.det(numpy.array([dx[ii],dy[ii],
                                                     dvx[ii],dvy[ii]])) for ii in range(nt)])
    breakt= 8.
    pts= list(ts[ts < breakt])
    pts.extend(list((ts[ts >= breakt])[::10]))
    pts= numpy.array(pts)
    pjacs= list(jacs[ts < breakt])
    pjacs.extend(list((jacs[ts >= breakt])[::10]))
    pjacs= numpy.array(pjacs)
    print integrator, numpy.mean(jacs)-1.
    bovy_plot.bovy_print(fig_width=3.25,fig_height=4.5)
    pyplot.subplot(4,1,4)
    bovy_plot.bovy_plot(pts/orbt,numpy.fabs(pjacs-1.),color='k',
                        loglog=True,gcf=True,
                        xrange=[0.5,norb],
                        yrange=[10.**-12.,10.**0.],
                        xlabel=r'$\mathrm{Number\ of\ orbital\ periods}$')
    bovy_plot.bovy_text(r'$\texttt{dopr54\_c}$',
                        top_left=True,size=14.)
    bovy_plot.bovy_text(0.1,10.**44.5,
                        r'$|\mathrm{Determinant\ of\ volume\ transformation}-1|$',
                        fontsize=16.,
                        rotation='vertical')
    ax= pyplot.gca()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter(r'$%0.f$'))
    ax.yaxis.set_ticks([10.**-12.,10.**-8.,10.**-4.,1.])
    nullfmt   = NullFormatter()         # no labels
    other_integrators= ['odeint',
                        'rk4_c','rk6_c']
    for ii,integrator in enumerate(other_integrators):
        start= time.time()
        o.integrate_dxdv([1.,0.,0.,0.],ts,MWPotential2014,
                         method=integrator,
                         rectIn=True,rectOut=True)
        dx= o.getOrbit_dxdv()[:,:]
        o.integrate_dxdv([0.,1.,0.,0.],ts,MWPotential2014,
                         method=integrator,
                         rectIn=True,rectOut=True)
        dy= o.getOrbit_dxdv()[:,:]
        o.integrate_dxdv([0.,0.,1.,0.],ts,MWPotential2014,
                         method=integrator,
                         rectIn=True,rectOut=True)
        dvx= o.getOrbit_dxdv()[:,:]
        o.integrate_dxdv([0.,0.,0.,1.],ts,MWPotential2014,
                         method=integrator,
                         rectIn=True,rectOut=True)
        dvy= o.getOrbit_dxdv()[:,:]
        print integrator, time.time()-start
        jacs= numpy.array([numpy.linalg.det(numpy.array([dx[jj],dy[jj],
                                                         dvx[jj],dvy[jj]])) for jj in range(nt)])
        pts= list(ts[ts < breakt])
        pts.extend(list((ts[ts >= breakt])[::10]))
        pts= numpy.array(pts)
        pjacs= list(jacs[ts < breakt])
        pjacs.extend(list((jacs[ts >= breakt])[::10]))
        pjacs= numpy.array(pjacs)
        print integrator, numpy.mean(jacs)-1.
        pyplot.subplot(4,1,ii+1)
        if integrator == 'odeint':
            yrange=[10.**-8.,10.**4.]
            yticks= [10.**-8.,10.**-4.,1.,10.**4.]
        else:
            yrange=[10.**-12.,10.**0.]
            yticks= [10.**-12.,10.**-8.,10.**-4.,1.]
        bovy_plot.bovy_plot(pts/orbt,numpy.fabs(pjacs-1.),color='k',
                            loglog=True,gcf=True,
                            xrange=[0.5,norb],
                            yrange=yrange)
        thisax= pyplot.gca()
        thisax.xaxis.set_major_formatter(nullfmt)
        thisax.yaxis.set_ticks(yticks)
        bovy_plot.bovy_text(r'$\texttt{%s}$' % (integrator.replace('_','\_')),
                            top_left=True,size=14.)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_liouville(sys.argv[1])
