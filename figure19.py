import os, os.path
import sys
import pickle
import numpy
from galpy.potential import MWPotential2014, rl
from galpy.potential import evaluatePotentials as evalPot
from galpy.orbit import Orbit
from galpy.actionAngle import estimateDeltaStaeckel
from galpy.util import bovy_plot, save_pickles
def calc_delta_MWPotential2014(savefilename,plotfilename):
    Lmin, Lmax= 0.01, 10.
    if not os.path.exists(savefilename):
        #Setup grid
        nL, nE= 101,101
        Ls= 10.**numpy.linspace(numpy.log10(Lmin),numpy.log10(Lmax),nL)
        #Integration times
        ts= numpy.linspace(0.,20.,1001)
        deltas= numpy.empty((nL,nE))
        Einf= evalPot(10.**12.,0.,MWPotential2014)
        print Einf
        for ii in range(nL):
            #Calculate Ec
            Rc= rl(MWPotential2014,Ls[ii])
            print ii, "Rc = ", Rc*8.
            Ec= evalPot(Rc,0.,MWPotential2014)+0.5*Ls[ii]**2./Rc**2.
            Es= Ec+(Einf-Ec)*10.**numpy.linspace(-2.,0.,nE)
            for jj in range(nE):
                #Setup an orbit with this energy and angular momentum
                Er= 2.*(Es[jj]-Ec) #Random energy times 2 = vR^2 + vz^2
                vR= numpy.sqrt(4./5.*Er)
                vz= numpy.sqrt(1./5.*Er)
                o= Orbit([Rc,vR,Ls[ii]/Rc,0.,vz,0.])
                o.integrate(ts,MWPotential2014,method='symplec4_c')
                deltas[ii,jj]= estimateDeltaStaeckel(o.R(ts),o.z(ts),
                                                     pot=MWPotential2014)
        #Save
        save_pickles(savefilename,deltas)
    else:
        savefile= open(savefilename,'rb')
        deltas= pickle.load(savefile)
        savefile.close()
    #Plot
    print numpy.nanmax(deltas)
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(deltas.T,origin='lower',cmap='jet',
                          xrange=[numpy.log10(Lmin),numpy.log10(Lmax)],
                          yrange=[-2,0.],
                          xlabel=r'$\log_{10} L$',
                          ylabel=r'$\log_{10}\left(\frac{E-E_c(L)}{E(\infty)-E_c(L)}\right)$',
                          colorbar=True,shrink=0.78,
                          zlabel=r'$\mathrm{Approximate\ focal\ length}$',
#                          interpolation='nearest',
                          vmin=0.,vmax=0.6)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    calc_delta_MWPotential2014(sys.argv[1],sys.argv[2])
