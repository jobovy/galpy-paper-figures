import sys
import numpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014, IsochronePotential
from galpy.potential import evaluatePotentials as evalPot
from galpy.actionAngle import actionAngleIsochrone, actionAngleSpherical, \
    estimateBIsochrone
from galpy.actionAngle import dePeriod
from galpy.util import bovy_plot
from matplotlib import pyplot
import matplotlib.ticker as ticker
from matplotlib.ticker import NullFormatter
def plot_aaspher_conservation(plotfilename1,plotfilename2):
    #Setup orbit
    E, Lz= -1.25, 0.6
    o= Orbit([0.8,0.3,Lz/0.8,0.,numpy.sqrt(2.*(E-evalPot(0.8,0.,MWPotential2014)-(Lz/0.8)**2./2.-0.3**2./2.)),0.])
    #Integrate the orbit to estimate an equivalent b
    nt= 1001
    ts= numpy.linspace(0.,20.,nt)
    o.integrate(ts,MWPotential2014,method='symplec4_c')
    b= estimateBIsochrone(o.R(ts),o.z(ts),pot=MWPotential2014)
    print b
    b= 0.3
    #Now integrate the orbit in the isochronePotential
    ip= IsochronePotential(normalize=1.,b=b)
    aAI= actionAngleIsochrone(ip=ip)
    orbt= 2.*numpy.pi/aAI.actionsFreqs(o)[4]    
    norb= 200.
    nt= 20001
    ts= numpy.linspace(0.,norb*orbt,nt)
    o.integrate(ts,ip,method='symplec4_c')
    #Calculate actions, frequencies, and angles
    jfa= aAI.actionsFreqsAngles(o.R(ts),o.vR(ts),o.vT(ts),
                                o.z(ts),o.vz(ts),o.phi(ts))
    dJs= numpy.fabs((jfa[0]-numpy.mean(jfa[0]))/numpy.mean(jfa[0]))
    dOrs= numpy.fabs((jfa[3]-numpy.mean(jfa[3]))/numpy.mean(jfa[3]))
    dOzs= numpy.fabs((jfa[5]-numpy.mean(jfa[5]))/numpy.mean(jfa[5]))
    print "frequencies", numpy.mean(dOrs), numpy.mean(dOzs)
    ar= dePeriod(numpy.reshape(jfa[6],(1,len(ts)))).flatten()
    az= dePeriod(numpy.reshape(jfa[8],(1,len(ts)))).flatten()
    danglers= numpy.fabs(ar-numpy.mean(jfa[3])*ts-jfa[6][0])/2./numpy.pi
    danglezs= numpy.fabs(az-numpy.mean(jfa[5])*ts-jfa[8][0])/2./numpy.pi
    #Break up
    breakt= 50.
    pts= parse_break(ts,ts < breakt)
    pdJs= parse_break(dJs,ts < breakt)
    pdanglers= parse_break(danglers,ts < breakt)
    pdanglezs= parse_break(danglezs,ts < breakt)
    #dAngles
    bovy_plot.bovy_print()
    pyplot.subplot(2,1,1)
    bovy_plot.bovy_plot(pts/orbt,
                        pdJs,
                        color='k',loglog=True,gcf=True,
                        xrange=[0.5,norb],
                        yrange=[10.**-12.,1.])
    bovy_plot.bovy_text(r'$\texttt{actionAngleIsochrone}$',
                        top_left=True,size=14.)
    ax= pyplot.gca()
    ax.yaxis.set_ticks([10.**-12.,10.**-8.,10.**-4.,1.])
    nullfmt   = NullFormatter()         # no labels
    ax.xaxis.set_major_formatter(nullfmt)
    #Same for actionAngleSpherical
    aAS= actionAngleSpherical(pot=ip)
    tts= ts[::1]
    jfa= aAS.actionsFreqsAngles(o.R(tts),o.vR(tts),o.vT(tts),
                                o.z(tts),o.vz(tts),o.phi(tts),
                                fixed_quad=True)
    #dJr
    dJs= numpy.fabs((jfa[0]-numpy.mean(jfa[0]))/numpy.mean(jfa[0]))
    dOrs= numpy.fabs((jfa[3]-numpy.mean(jfa[3]))/numpy.mean(jfa[3]))
    dOzs= numpy.fabs((jfa[5]-numpy.mean(jfa[5]))/numpy.mean(jfa[5]))
    print "frequencies", numpy.mean(dOrs), numpy.mean(dOzs)
    #dAngles
    ar= dePeriod(numpy.reshape(jfa[6],(1,len(tts)))).flatten()
    az= dePeriod(numpy.reshape(jfa[8],(1,len(tts)))).flatten()
    danglers= numpy.fabs(ar-numpy.mean(jfa[3])*tts-jfa[6][0])/2./numpy.pi
    danglezs= numpy.fabs(az-numpy.mean(jfa[5])*tts-jfa[8][0])/2./numpy.pi
    print numpy.mean(danglers)
    print numpy.mean(danglezs)
    ptts= parse_break(tts,tts < breakt)
    pdJs= parse_break(dJs,tts < breakt)
    pyplot.subplot(2,1,2)
    bovy_plot.bovy_plot(ptts/orbt,
                        pdJs,
                        color='k',loglog=True,gcf=True,
                        xrange=[0.5,norb],
                        yrange=[10.**-12.,1.],
                        xlabel=r'$\mathrm{Number\ of\ orbital\ periods}$')
    bovy_plot.bovy_text(r'$\texttt{actionAngleSpherical}$',
                        top_left=True,size=14.)
    bovy_plot.bovy_text(0.175,10.**2.,r'$\left|\Delta J_R / J_R\right|$',
                        fontsize=16.,
                        rotation='vertical')
    ax= pyplot.gca()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter(r'$%0.f$'))
    ax.yaxis.set_ticks([10.**-12.,10.**-8.,10.**-4.,1.])
    bovy_plot.bovy_end_print(plotfilename1)    
    #Now plot the deviations in the angles
    bovy_plot.bovy_print()
    pyplot.subplot(2,1,1)
    liner= bovy_plot.bovy_plot(pts/orbt,
                               pdanglers,
                               color='k',ls='-',loglog=True,gcf=True,
                               xrange=[0.5,norb],
                               yrange=[10.**-12.,1.])
    linez= bovy_plot.bovy_plot(pts/orbt,
                               pdanglezs,
                               color='k',ls='--',overplot=True)
    legend1= pyplot.legend((liner[0],linez[0]),
                           (r'$\theta_R$',
                            r'$\theta_z$'),
                           loc='lower right',#bbox_to_anchor=(.91,.375),
                           numpoints=2,
                           prop={'size':14},
                  frameon=False)
    bovy_plot.bovy_text(r'$\texttt{actionAngleIsochrone}$',
                        top_left=True,size=14.)
    ax= pyplot.gca()
    ax.yaxis.set_ticks([10.**-12.,10.**-8.,10.**-4.,1.])
    nullfmt   = NullFormatter()         # no labels
    ax.xaxis.set_major_formatter(nullfmt)
    #Same for Spherical
    pdanglers= parse_break(danglers,tts < breakt)
    pdanglezs= parse_break(danglezs,tts < breakt)
    pyplot.subplot(2,1,2)
    bovy_plot.bovy_plot(ptts/orbt,
                        pdanglers,
                        color='k',ls='-',loglog=True,gcf=True,
                        xrange=[0.5,norb],
                        yrange=[10.**-12.,1.],
                        xlabel=r'$\mathrm{Number\ of\ orbital\ periods}$')
    bovy_plot.bovy_plot(ptts/orbt,
                        pdanglezs,
                        color='k',ls='--',overplot=True)
    bovy_plot.bovy_text(r'$\texttt{actionAngleSpherical}$',
                        top_left=True,size=14.)
    bovy_plot.bovy_text(0.175,10.**4.,r'$\left|\Delta \theta_{R,z} / 2\,\pi\right|$',
                        fontsize=16.,
                        rotation='vertical')
    ax= pyplot.gca()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter(r'$%0.f$'))
    ax.yaxis.set_ticks([10.**-12.,10.**-8.,10.**-4.,1.])
    bovy_plot.bovy_end_print(plotfilename2)    
    return None

def parse_break(quant,b4indx):
    pquant= list(quant[b4indx])
    pquant.extend(list((quant[True-b4indx])[::10]))
    pquant= numpy.array(pquant)
    return pquant

if __name__ == '__main__':
    plot_aaspher_conservation(sys.argv[1],
                              sys.argv[2])
