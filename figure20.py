import sys
import numpy
from galpy.potential import MWPotential2014
from galpy.potential import evaluatePotentials as evalPot
from galpy.orbit import Orbit
from galpy.actionAngle import actionAngleAdiabatic, actionAngleAdiabaticGrid, \
    actionAngleStaeckel, actionAngleStaeckelGrid
from galpy.util import bovy_plot
from matplotlib import pyplot
def plot_aagrid(plotfilename1,plotfilename2):
    #Setup orbit
    E, Lz= -1.25, 0.6
    o= Orbit([0.8,0.3,Lz/0.8,0.,numpy.sqrt(2.*(E-evalPot(0.8,0.,MWPotential2014)-(Lz/0.8)**2./2.-0.3**2./2.)),0.])
    delta= 0.434
    #Integrate the orbit, setup Staeckel already to calculate the period
    aAS= actionAngleStaeckel(pot=MWPotential2014,delta=delta,c=True)
    orbt= 2.*numpy.pi/aAS.actionsFreqs(o)[4]  
    norb= 5.
    nt= 501
    ts= numpy.linspace(0.,norb*orbt,nt)
    o.integrate(ts,MWPotential2014,method='symplec4_c')
    #First do adiabatic
    aAA= actionAngleAdiabatic(pot=MWPotential2014,gamma=1.,c=True)
    aAAG= actionAngleAdiabaticGrid(pot=MWPotential2014,gamma=1.,c=True,
                                   nR=31,nEz=31,nEr=51,nLz=51)
    jfa= aAA(o.R(ts),o.vR(ts),o.vT(ts),o.z(ts),o.vz(ts),o.phi(ts))
    jfag= aAAG(o.R(ts),o.vR(ts),o.vT(ts),o.z(ts),o.vz(ts),o.phi(ts))
    #First do adiabatic
    #aAS already setup
    aASG= actionAngleStaeckelGrid(pot=MWPotential2014,delta=delta,c=True,
                                  nE=51,npsi=51,nLz=51)
    jfs= aAS(o.R(ts),o.vR(ts),o.vT(ts),o.z(ts),o.vz(ts),o.phi(ts))
    jfsg= aASG(o.R(ts),o.vR(ts),o.vT(ts),o.z(ts),o.vz(ts),o.phi(ts))
    bovy_plot.bovy_print()
    line1= bovy_plot.bovy_plot(jfa[0],jfa[2],'r.',
                               xrange=[0.045,0.055],
                               yrange=[0.0075,0.011],
                               xlabel=r'$J_R$',ylabel=r'$J_z$',zorder=2)
    line2= bovy_plot.bovy_plot(jfag[0],jfag[2],'rx',overplot=True,zorder=1)
    bovy_plot.bovy_plot(jfs[0],jfs[2],'k,',overplot=True)
    pyplot.legend((line1[0],line2[0]),
                  (r'$\mathrm{\texttt{actionAngleAdiabatic}}$',
                   r'$\mathrm{\texttt{actionAngleAdiabaticGrid}}$',),
                  loc='upper right',#bbox_to_anchor=(.91,.375),
                  numpoints=1,
                  prop={'size':14},
                  frameon=False)
    bovy_plot.bovy_end_print(plotfilename1)
    #Zoom of Staeckel
    line1= bovy_plot.bovy_plot(jfs[0],jfs[2],'k.',
                               xrange=[0.05025,0.05145],
                               yrange=[0.0086,0.00933],
                               xlabel=r'$J_R$',ylabel=r'$J_z$')
    line2= bovy_plot.bovy_plot(jfsg[0],jfsg[2],'kx',overplot=True)
    pyplot.legend((line1[0],line2[0]),
                  (r'$\mathrm{\texttt{actionAngleStaeckel}}$',
                   r'$\mathrm{\texttt{actionAngleStaeckelGrid}}$',),
                  loc='upper right',#bbox_to_anchor=(.91,.375),
                  numpoints=1,
                  prop={'size':14},
                  frameon=False)
    bovy_plot.bovy_end_print(plotfilename2)
    return None

if __name__ == '__main__':
    plot_aagrid(sys.argv[1],sys.argv[2])
