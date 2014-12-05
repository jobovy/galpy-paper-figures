import sys
import numpy
from galpy.potential import MWPotential2014
from galpy.potential import evaluatePotentials as evalPot
from galpy.orbit import Orbit
from galpy.util import bovy_plot
def plot_twoorbits(plotfilename1,plotfilename2):
    #First orbit
    E, Lz= -1.25, 0.6
    o1= Orbit([0.8,0.,Lz/0.8,0.,numpy.sqrt(2.*(E-evalPot(0.8,0.,MWPotential2014)-(Lz/0.8)**2./2.)),0.])
    ts= numpy.linspace(0.,100.,2001)
    o1.integrate(ts,MWPotential2014)
    print "First orbit: E, L = %f,%f" % (o1.E(),o1.L()[:,2])
    o2= Orbit([0.8,0.3,Lz/0.8,0.,numpy.sqrt(2.*(E-evalPot(0.8,0.,MWPotential2014)-(Lz/0.8)**2./2.-0.3**2./2.)),0.])
    o2.integrate(ts,MWPotential2014)
    print "Second orbit: E, L = %f,%f" % (o2.E(),o2.L()[:,2])
    print "First orbit: zmax = %f" % (o1.zmax())
    print "Second orbit: zmax = %f" % (o2.zmax())
    o1.plot(xrange=[0.3,1.],yrange=[-0.2,0.2],color='k')
    bovy_plot.bovy_end_print(plotfilename1)
    o2.plot(xrange=[0.3,1.],yrange=[-0.2,0.2],color='k')
    bovy_plot.bovy_end_print(plotfilename2)
    return (o1,o2)

def plot_poincare(o1,o2,plotfilename):
    ts= numpy.linspace(0.,1000.,20001)
    o1.integrate(ts,MWPotential2014)
    o2.integrate(ts,MWPotential2014)
    #First orbit
    Rs= o1.R(ts); zs= o1.z(ts); vRs= o1.vR(ts)
    sect1Rs, sect1vRs= surface_section(Rs,zs,vRs)
    #Second orbit
    Rs= o2.R(ts); zs= o2.z(ts); vRs= o2.vR(ts)
    sect2Rs, sect2vRs= surface_section(Rs,zs,vRs)
    bovy_plot.bovy_print(fig_height=3.5)
    bovy_plot.bovy_plot(sect1Rs,sect1vRs,'bo',mec='none',
                        xlabel=r'$R$',ylabel=r'$v_R$',
                        xrange=[0.3,1.],
                        yrange=[-0.69,0.69])
    bovy_plot.bovy_plot(sect2Rs,sect2vRs,'yo',mec='none',overplot=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None

def surface_section(Rs,zs,vRs):
    #Find points where the orbit crosses z from negative to positive
    shiftzs= numpy.roll(zs,-1)
    indx= (zs[:-1] < 0.)*(shiftzs[:-1] > 0.)
    return (Rs[:-1][indx],vRs[:-1][indx])

if __name__ == '__main__':
    o1,o2= plot_twoorbits(sys.argv[1],sys.argv[2])
    plot_poincare(o1,o2,sys.argv[3])
