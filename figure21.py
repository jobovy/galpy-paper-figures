import sys
import numpy
from galpy.potential import IsochronePotential
from galpy.orbit import Orbit
from galpy.actionAngle import actionAngleIsochrone
from galpy.util import bovy_plot
from galpy.potential import Potential
def smoothInterp(t,dt,tform):
    """Smooth interpolation in time, following Dehnen (2000)"""
    if t < tform: smooth= 0.
    elif t > (tform+dt): smooth= 1. 
    else:
        xi= 2.*(t-tform)/dt-1.
        smooth= (3./16.*xi**5.-5./8*xi**3.+15./16.*xi+.5)
    return smooth
class TimeInterpPotential(Potential):
    """Potential that smoothly interpolates in time between two static potentials"""
    def __init__(self,pot1,pot2,dt=100.,tform=50.):
        """pot1= potential for t < tform, pot2= potential for t > tform+dt, dt: time over which to turn on pot2, tform: time at which the interpolation is switched on"""
        Potential.__init__(self,amp=1.)
        self._pot1= pot1
        self._pot2= pot2
        self._tform= tform
        self._dt= dt
        return None
    def _Rforce(self,R,z,phi=0.,t=0.):
        smooth= smoothInterp(t,self._dt,self._tform)
        return (1.-smooth)*self._pot1.Rforce(R,z)+smooth*self._pot2.Rforce(R,z)
    def _zforce(self,R,z,phi=0.,t=0.):
        smooth= smoothInterp(t,self._dt,self._tform)
        return (1.-smooth)*self._pot1.zforce(R,z)+smooth*self._pot2.zforce(R,z)

def illustrate_adiabatic_invariance(plotfilename1,plotfilename2):
    # Initialize two different IsochronePotentials
    ip1= IsochronePotential(normalize=1.,b=1.)
    ip2= IsochronePotential(normalize=0.5,b=1.)
    # Use TimeInterpPotential to interpolate smoothly between the two
    tip= TimeInterpPotential(ip1,ip2,dt=100.,tform=50.)
    # Integrate the orbit, in three parts
    # 1) Orbit in the first isochrone potential
    o1= Orbit([1.,0.1,1.1,0.0,0.1,0.])
    ts= numpy.linspace(0.,50.,1001)
    o1.integrate(ts,tip)
    bovy_plot.bovy_print()
    o1.plot(d1='x',d2='y',xrange=[-1.6,1.6],yrange=[-1.6,1.6],color='b',
            gcf=True)
    # 2) Orbit in the transition
    o2= o1(ts[-1]) # Last time step = initial time step of the next integration
    ts2= numpy.linspace(50.,150.,1001)
    o2.integrate(ts2,tip)
    o2.plot(d1='x',d2='y',overplot=True,color='g')
    # 3) Orbit in the second isochrone potential
    o3= o2(ts2[-1])
    ts3= numpy.linspace(150.,200.,1001)
    o3.integrate(ts3,ip2)
    o3.plot(d1='x',d2='y',overplot=True,color='r')
    bovy_plot.bovy_end_print(plotfilename1)
    # Also plot the R,z projection
    bovy_plot.bovy_print(fig_height=2.3333)
    o1.plot(d1='R',d2='z',xrange=[0.9,1.65],yrange=[-.175,.175],color='b',
            gcf=True)
    o2.plot(d1='R',d2='z',overplot=True,color='g')
    o3.plot(d1='R',d2='z',overplot=True,color='r')
    bovy_plot.bovy_end_print(plotfilename2)   
    # Now we calculate the energy, eccentricity, mean radius, and maximum height
    print o1.E(pot=ip1), o1.e(), 0.5*(o1.rperi()+o1.rap()), o1.zmax()
    print o3.E(pot=ip2), o3.e(), 0.5*(o3.rperi()+o3.rap()), o3.zmax()
   # The orbit has clearly moved to larger radii, the actions are however conserved
    aAI1= actionAngleIsochrone(ip=ip1)
    aAI2= actionAngleIsochrone(ip=ip2)
    print aAI1(o1)
    print aAI2(o3)
    return None

if __name__ == '__main__':
    illustrate_adiabatic_invariance(sys.argv[1],sys.argv[2])
