import sys
import numpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from galpy.potential import evaluatePotentials as evalPot
from galpy.actionAngle import actionAngleAdiabatic, actionAngleStaeckel, \
    estimateDeltaStaeckel, actionAngleIsochroneApprox
from galpy.actionAngle import dePeriod
from galpy.util import bovy_plot
from matplotlib import pyplot
import matplotlib.ticker as ticker
from matplotlib.ticker import NullFormatter
def plot_aaaxi_conservation(plotfilename1,plotfilename2):
    #Setup orbit
    E, Lz= -1.25, 0.6
    o= Orbit([0.8,0.3,Lz/0.8,0.,numpy.sqrt(2.*(E-evalPot(0.8,0.,MWPotential2014)-(Lz/0.8)**2./2.-0.3**2./2.)),0.])
    #Integrate the orbit to estimate an equivalent delta
    nt= 1001
    ts= numpy.linspace(0.,20.,nt)
    o.integrate(ts,MWPotential2014,method='symplec4_c')
    delta= estimateDeltaStaeckel(o.R(ts),o.z(ts),pot=MWPotential2014)
    print delta
    delta= 0.434
    #Now integrate the orbit
    aAS= actionAngleStaeckel(pot=MWPotential2014,delta=delta,c=True)
    orbt= 2.*numpy.pi/aAS.actionsFreqs(o)[4]  
    norb= 200.
    nt= 20001
    ts= numpy.linspace(0.,norb*orbt,nt)
    o.integrate(ts,MWPotential2014,method='symplec4_c')
    #Calculate actions, frequencies, and angles in the adiabatic approximation
    aAA= actionAngleAdiabatic(pot=MWPotential2014,c=True)
    jfa= aAA(o.R(ts),o.vR(ts),o.vT(ts),
             o.z(ts),o.vz(ts),o.phi(ts))
    print "actions", numpy.mean(jfa[0][True-numpy.isnan(jfa[0])]), numpy.mean(jfa[2])
    mjr= numpy.mean(jfa[0][True-numpy.isnan(jfa[0])])
    mjz= numpy.mean(jfa[2][True-numpy.isnan(jfa[2])])
    dJs= numpy.fabs((jfa[0]-mjr)/mjr)
    dJzs= numpy.fabs((jfa[2]-mjz)/mjz)
    print numpy.mean(dJs[True-numpy.isnan(dJs)]), numpy.mean(dJzs[True-numpy.isnan(dJzs)])
    #Break up
    breakt= 50.
    pts= parse_break(ts,ts < breakt)
    pdJs= parse_break(dJs,ts < breakt)
    pdJzs= parse_break(dJzs,ts < breakt)
    bovy_plot.bovy_print(fig_height=7.)
    pyplot.subplot(3,1,1)
    liner= bovy_plot.bovy_plot(pts/orbt,
                               pdJs,
                               color='k',loglog=True,gcf=True,
                               xrange=[0.5,norb],
                               yrange=[10.**-6.,1.])
    linez= bovy_plot.bovy_plot(pts/orbt,
                               pdJzs,
                               color='k',ls='--',overplot=True)
    bovy_plot.bovy_text(r'$\texttt{actionAngleAdiabatic}$',
                        top_left=True,size=14.)
    pyplot.legend((liner[0],linez[0]),
                  (r'$J_R$',
                   r'$J_z$'),
                  loc='lower left',#bbox_to_anchor=(.91,.375),
                  numpoints=2,
                  prop={'size':13},
                  frameon=False)
    ax= pyplot.gca()
    ax.yaxis.set_ticks([10.**-6.,10.**-4.,10.**-2.,1.])
    nullfmt   = NullFormatter()         # no labels
    ax.xaxis.set_major_formatter(nullfmt)
    #Same for actionAngleStaeckel
    jfa= aAS.actionsFreqsAngles(o.R(ts),o.vR(ts),o.vT(ts),
                                o.z(ts),o.vz(ts),o.phi(ts))
    print "actions", numpy.mean(jfa[0]), numpy.mean(jfa[2])
    #dJr
    dJs= numpy.fabs((jfa[0]-numpy.mean(jfa[0]))/numpy.mean(jfa[0]))
    dJzs= numpy.fabs((jfa[2]-numpy.mean(jfa[2]))/numpy.mean(jfa[2]))
    print numpy.mean(dJs), numpy.mean(dJzs)
    #dO
    dOrs= numpy.fabs((jfa[3]-numpy.mean(jfa[3]))/numpy.mean(jfa[3]))
    dOps= numpy.fabs((jfa[4]-numpy.mean(jfa[4]))/numpy.mean(jfa[4]))
    dOzs= numpy.fabs((jfa[5]-numpy.mean(jfa[5]))/numpy.mean(jfa[5]))
    print numpy.mean(dOrs), numpy.mean(dOps), numpy.mean(dOzs)
    #dAngles
    ar= dePeriod(numpy.reshape(jfa[6],(1,len(ts)))).flatten()
    ap= dePeriod(numpy.reshape(jfa[7],(1,len(ts)))).flatten()
    az= dePeriod(numpy.reshape(jfa[8],(1,len(ts)))).flatten()
    danglers= numpy.fabs(ar-numpy.mean(jfa[3])*ts-jfa[6][0])/2./numpy.pi
    dangleps= numpy.fabs(ap-numpy.mean(jfa[4])*ts-jfa[7][0])/2./numpy.pi
    danglezs= numpy.fabs(az-numpy.mean(jfa[5])*ts-jfa[8][0])/2./numpy.pi
    print numpy.mean(danglers)
    print numpy.mean(dangleps)
    print numpy.mean(danglezs)
    pts= parse_break(ts,ts < breakt)
    pdJs= parse_break(dJs,ts < breakt)
    pdJzs= parse_break(dJzs,ts < breakt)
    pdanglers= parse_break(danglers,ts < breakt)
    pdangleps= parse_break(dangleps,ts < breakt)
    pdanglezs= parse_break(danglezs,ts < breakt)
    pyplot.subplot(3,1,2)
    bovy_plot.bovy_plot(pts/orbt,
                        pdJs,
                        color='k',loglog=True,gcf=True,
                        xrange=[0.5,norb],
                        yrange=[10.**-6.,1.],
                        ylabel=r'$\left|\Delta J_{R,z} / J_{R,z}\right|$')
    bovy_plot.bovy_plot(pts/orbt,
                        pdJzs,
                        color='k',ls='--',overplot=True)
    bovy_plot.bovy_text(r'$\texttt{actionAngleStaeckel}$',
                        top_left=True,size=14.)
    ax= pyplot.gca()
    ax.yaxis.set_ticks([10.**-6.,10.**-4.,10.**-2.,1.])
    nullfmt   = NullFormatter()         # no labels
    ax.xaxis.set_major_formatter(nullfmt)
    #Same for actionAngleIsochroneApprox
    tts= ts[::1]
    aAIA= actionAngleIsochroneApprox(pot=MWPotential2014,b=0.3)
    #Process separately to avoid MemoryError
    aiajr, aiajp, aiajz, aiaor, aiaop, aiaoz, aiaar, aiaap, aiaaz= [],[],[],[],[],[],[],[],[]
    for ii in range(len(tts)):
        tjfa= aAIA.actionsFreqsAngles(o.R(tts[ii]),o.vR(tts[ii]),o.vT(tts[ii]),
                                      o.z(tts[ii]),o.vz(tts[ii]),o.phi(tts[ii]))
        aiajr.append(tjfa[0])
        aiajp.append(tjfa[1])
        aiajz.append(tjfa[2])
        aiaor.append(tjfa[3])
        aiaop.append(tjfa[4])
        aiaoz.append(tjfa[5])
        aiaar.append(tjfa[6])
        aiaap.append(tjfa[7])
        aiaaz.append(tjfa[8])
    jfa= (numpy.array(aiajr),numpy.array(aiajp),numpy.array(aiajz),
          numpy.array(aiaor),numpy.array(aiaop),numpy.array(aiaoz),
          numpy.array(aiaar).flatten(),numpy.array(aiaap).flatten(),
          numpy.array(aiaaz).flatten())
    print "actions", numpy.mean(jfa[0]), numpy.mean(jfa[2])
    #dJr
    dJs= numpy.fabs((jfa[0]-numpy.mean(jfa[0]))/numpy.mean(jfa[0]))
    dJzs= numpy.fabs((jfa[2]-numpy.mean(jfa[2]))/numpy.mean(jfa[2]))
    print numpy.mean(dJs), numpy.mean(dJzs)
    #dO
    dOrs= numpy.fabs((jfa[3]-numpy.mean(jfa[3]))/numpy.mean(jfa[3]))
    dOps= numpy.fabs((jfa[4]-numpy.mean(jfa[4]))/numpy.mean(jfa[4]))
    dOzs= numpy.fabs((jfa[5]-numpy.mean(jfa[5]))/numpy.mean(jfa[5]))
    print numpy.mean(dOrs), numpy.mean(dOps), numpy.mean(dOzs)
    if len(tts) > 5000:
        ar= dePeriod(numpy.reshape(jfa[6],(1,len(tts)))).flatten()
        ap= dePeriod(numpy.reshape(jfa[7],(1,len(tts)))).flatten()
        az= dePeriod(numpy.reshape(jfa[8],(1,len(tts)))).flatten()
        danglers= numpy.fabs(ar-numpy.mean(jfa[3])*tts-jfa[6][0])/2./numpy.pi
        dangleps= numpy.fabs(ap-numpy.mean(jfa[4])*tts-jfa[7][0])/2./numpy.pi
        danglezs= numpy.fabs(az-numpy.mean(jfa[5])*tts-jfa[8][0])/2./numpy.pi
    else:
        danglers= numpy.fabs(jfa[6]-numpy.mod(numpy.mean(jfa[3])*tts+jfa[6][0],2.*numpy.pi))/2./numpy.pi
        dangleps= numpy.fabs(jfa[7]-numpy.mod(numpy.mean(jfa[4])*tts+jfa[7][0],2.*numpy.pi))/2./numpy.pi
        danglezs= numpy.fabs(jfa[8]-numpy.mod(numpy.mean(jfa[5])*tts+jfa[8][0],2.*numpy.pi))/2./numpy.pi
    print numpy.mean(danglers)
    print numpy.mean(dangleps)
    print numpy.mean(danglezs)
    ptts= parse_break(tts,tts < breakt)
    pdJs= parse_break(dJs,tts < breakt)
    pdJzs= parse_break(dJzs,tts < breakt)
    pyplot.subplot(3,1,3)
    bovy_plot.bovy_plot(ptts/orbt,
                        pdJs,
                        color='k',loglog=True,gcf=True,
                        xrange=[0.5,norb],
                        yrange=[10.**-6.,1.],
                        xlabel=r'$\mathrm{Number\ of\ orbital\ periods}$')
    bovy_plot.bovy_plot(ptts/orbt,
                        pdJzs,
                        color='k',ls='--',overplot=True)
    bovy_plot.bovy_text(r'$\texttt{actionAngleIsochroneApprox}$',
                        top_left=True,size=14.)
    ax= pyplot.gca()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter(r'$%0.f$'))
    ax.yaxis.set_ticks([10.**-6.,10.**-4.,10.**-2.,1.])
    bovy_plot.bovy_end_print(plotfilename1)    
    #Now plot the deviations in the angles
    bovy_plot.bovy_print()
    pyplot.subplot(2,1,1)
    liner= bovy_plot.bovy_plot(pts/orbt,
                               pdanglers,
                               color='k',ls='-',loglog=True,gcf=True,
                               xrange=[0.5,norb],
                               yrange=[10.**-4.,1.])
    linep= bovy_plot.bovy_plot(pts/orbt,
                               pdangleps,
                               color='k',ls='-.',overplot=True)
    linez= bovy_plot.bovy_plot(pts/orbt,
                               pdanglezs,
                               color='k',ls='--',overplot=True)
    pyplot.legend((liner[0],linep[0],linez[0]),
                  (r'$\theta_R$',
                   r'$\theta_\phi$',
                   r'$\theta_z$'),
                  loc='lower right',#bbox_to_anchor=(.91,.375),
                  numpoints=2,
                  prop={'size':14},
                  frameon=False)
    bovy_plot.bovy_text(r'$\texttt{actionAngleStaeckel}$',
                        top_left=True,size=14.)
    ax= pyplot.gca()
    ax.yaxis.set_ticks([10.**-4.,10.**-2.,1.])
    nullfmt   = NullFormatter()         # no labels
    ax.xaxis.set_major_formatter(nullfmt)
    #Same for Spherical
    pdanglers= parse_break(danglers,tts < breakt)
    pdangleps= parse_break(dangleps,tts < breakt)
    pdanglezs= parse_break(danglezs,tts < breakt)
    pyplot.subplot(2,1,2)
    bovy_plot.bovy_plot(ptts/orbt,
                        pdanglers,
                        color='k',ls='-',loglog=True,gcf=True,
                        xrange=[0.5,norb],
                        yrange=[10.**-8.,1.],
                        xlabel=r'$\mathrm{Number\ of\ orbital\ periods}$')
    bovy_plot.bovy_plot(ptts/orbt,
                        pdangleps,
                        color='k',ls='-.',overplot=True)
    bovy_plot.bovy_plot(ptts/orbt,
                        pdanglezs,
                        color='k',ls='--',overplot=True)
    bovy_plot.bovy_text(r'$\texttt{actionAngleIsochroneApprox}$',
                        top_left=True,size=14.)
    bovy_plot.bovy_text(0.175,10.**2.,r'$\left|\Delta \theta_{R,\phi,z} / 2\,\pi\right|$',
                        fontsize=16.,
                        rotation='vertical')
    ax= pyplot.gca()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter(r'$%0.f$'))
    ax.yaxis.set_ticks([10.**-8.,10.**-6.,10.**-4.,10.**-2.,1.])
    bovy_plot.bovy_end_print(plotfilename2)    
    return None

def parse_break(quant,b4indx):
    pquant= list(quant[b4indx])
    pquant.extend(list((quant[True-b4indx])[::10]))
    pquant= numpy.array(pquant)
    return pquant

if __name__ == '__main__':
    plot_aaaxi_conservation(sys.argv[1],
                            sys.argv[2])
