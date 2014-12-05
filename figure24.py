###############################################################################
# Axisymmetric kinematics plots
###############################################################################
import sys
import os, os.path
import cPickle as pickle
import math as m
import numpy
from matplotlib import pyplot
from matplotlib.ticker import LogFormatter
from galpy.df import dehnendf, shudf
from galpy.util import bovy_plot

class myLogFormatter(LogFormatter):
	def __call__(self,x,pos=None):
	    'Return the format for tick val *x* at position *pos*'
	    vmin, vmax = self.axis.get_view_interval()
	    d = abs(vmax - vmin)
	    b = self._base
	    if x == 0.0:
		    return '0'
	    sign = numpy.sign(x)
	    # only label the decades
	    fx = m.log(abs(x))/m.log(b)
	    isDecade = is_close_to_int(fx)
	    if not isDecade and self.labelOnlyBase: s = ''
	    elif x>10000: s= '%1.e'%x
	    elif x<1: s =  '%.2f'%x
	    else        : s =  self.pprint_val(x, d)
	    s= r'$'+s+'$'
	    if sign == -1:
		    s =  '-%s' % s
	    return self.fix_minus(s)
def vaOort(plotfilename):
    justCalc= False
    noTheory= True
    #Load all the dfs
    so= numpy.exp(numpy.linspace(numpy.log(0.0125),numpy.log(0.2),21))
    ndfs= 9
    dftypes= ['dehnen' for ii in range(ndfs)]
    betas= [0 for ii in range(ndfs)]
    hrs= [1./3. for ii in range(ndfs)]
    hss= [1. for ii in range(ndfs)]
    symbols=['o','^','s','p','D','1','2','3','4']
    colors= ['k' for ii in range(ndfs)]
    labels= [r'$\mathrm{Fiducial}\ (' for ii in range(ndfs)]
    #set up different types
    dftypes[1]= 'shu'
    labels[1]= r'$\mathrm{Shu\ DF}\ ('
    hrs[2]= 1./4.
    labels[2]= r'$h_R = R_0 / 4\ ('
    hrs[3]= 1./2.
    labels[3]= r'$h_R = R_0 / 2\ ('
    hss[4]= 2./3.
    labels[4]= r'$h_\sigma = 2 R_0 / 3\ ('
    betas[5]= -0.1
    labels[5]= r'$\beta = -0.1\ ('
    betas[6]= 0.1
    labels[6]= r'$\beta = \phantom{-}0.1\ ('
    betas[7]= -0.2
    labels[7]= r'$\beta = -0.2\ ('
    betas[8]= 0.2
    labels[8]= r'$\beta = \phantom{-}0.2\ ('
    
    dfs= []
    for ii in range(ndfs):
	    print "Working on DF %i ..." % ii
	    thesedfs= []
	    for jj in range(len(so)):
		    sys.stdout.write('\r'+"Working on %i out of %i" % \
				     (jj+1,len(so)))
		    sys.stdout.flush()
		    if dftypes[ii] == 'dehnen':
			    thesedfs.append(dehnendf(beta=betas[ii],
						     correct=True,niter=20,
						     profileParams=(hrs[ii],hss[ii],so[jj])))
		    else:
			    thesedfs.append(shudf(beta=betas[ii],
						  correct=True,niter=20,
						  profileParams=(hrs[ii],hss[ii],so[jj])))
	    dfs.append(thesedfs)
	    sys.stdout.write('\n')
	    sys.stdout.flush()

    #Make plots for all
    #1) asymmetric drift
    print "Working on asymmetric drift ..."
    vafilename= 'axi_va.sav'
    if os.path.exists(vafilename):
	    vafile= open(vafilename,'rb')
	    vas= pickle.load(vafile)
	    vafile.close()
	    precalc=True
    else:
	    precalc=False
    if not justCalc:
	    bovy_plot.bovy_print()
    ii= 0
    if not precalc:
	    print "Working on DF %i ..." % ii
	    vas= []
	    va= numpy.zeros(len(so))
	    for jj in range(len(so)):
		    sys.stdout.write('\r'+"Working on %i out of %i" % \
				     (jj+1,len(so)))
		    sys.stdout.flush()
		    va[jj]= 1.-dfs[ii][jj].meanvT(1.)
	    vas.append(va)
	    sys.stdout.write('\n')
	    sys.stdout.flush()
    ks= [numpy.mean(so**2./vas[ii])]
    for ii in range(1,ndfs):
	    print "Working on DF %i ..." % ii
	    if not precalc:
		    va= numpy.zeros(len(so))
		    for jj in range(len(so)):
			    sys.stdout.write('\r'+"Working on %i out of %i" % \
					     (jj+1,len(so)))
			    sys.stdout.flush()
			    va[jj]= 1.-dfs[ii][jj].meanvT(1.)
		    vas.append(va)
		    sys.stdout.write('\n')
		    sys.stdout.flush()
	    ks.append(numpy.mean(so**2./vas[ii]))
    if not precalc:
	    vafile= open(vafilename,'wb')
	    pickle.dump(vas,vafile)
	    vafile.close()

    #2) oortA
    print "Working on Oort ..."
    oortfilename= 'axi_oort.sav'
    if os.path.exists(oortfilename):
	    oortfile= open(oortfilename,'rb')
	    oortAs= pickle.load(oortfile)
	    oortBs= pickle.load(oortfile)
	    oortfile.close()
	    precalc= True
    else:
	    precalc= False
	    oortAs= []
	    oortBs= []
    if not justCalc:
	    bovy_plot.bovy_print()
    ii= 0
    if not precalc:
	    print "Working on DF %i ..." % ii
	    oortA= numpy.zeros(len(so))
	    oortB= numpy.zeros(len(so))
	    for jj in range(len(so)):
		    sys.stdout.write('\r'+"Working on %i out of %i" % \
				     (jj+1,len(so)))
		    sys.stdout.flush()
		    oortA[jj]= dfs[ii][jj].oortA(1.)-.5
		    oortB[jj]= .5+dfs[ii][jj].oortB(1.)
	    oortAs.append(oortA)
	    oortBs.append(oortB)
	    sys.stdout.write('\n')
	    sys.stdout.flush()
    va= vas[ii]
    xrange=[0.0001,0.2]
    yrange=[0.00001,.2]
    lines= []
    if not justCalc:
	    Apos= numpy.all((-oortAs[ii]+betas[ii]/2.) > 0.)
	    Aneg= numpy.all((-oortAs[ii]+betas[ii]/2.) < 0.)
	    Bpos= numpy.all((oortBs[ii]+betas[ii]/2.) > 0.)
	    Bneg= numpy.all((oortBs[ii]+betas[ii]/2.) < 0.)
	    if Apos: labels[ii]+= '-,'
	    elif Aneg: labels[ii]+= '+,'
	    else: labels[ii]+= '?,'
	    if Bpos: labels[ii]+= '+)$'
	    elif Bneg: labels[ii]+= '-)$'
	    else: labels[ii]+= '?)'
	    lines.append(bovy_plot.bovy_plot(va,
					     do_fabs(-oortAs[ii]-betas[ii]/2.),
					     loglog=True,
					     xrange=xrange,
					     yrange=yrange,
					     xlabel=r'$v_a / v_0$',
					     ylabel=r'$|\Delta \mathrm{Oort}\ A\,,B|$',
					     ls='none',marker=symbols[ii],color=colors[ii],mfc='none')[0])
	    bovy_plot.bovy_plot(va,
				do_fabs(oortBs[ii]+betas[ii]/2.),
				loglog=True,overplot=True,
				ls='none',marker=symbols[ii],mec='k',mfc='none')
	    ax= pyplot.gca()
	    ax.xaxis.set_major_formatter(myLogFormatter(base=10.,labelOnlyBase=False))
	    ax.yaxis.set_major_formatter(myLogFormatter(base=10.,labelOnlyBase=False))
	    ax.yaxis.set_minor_formatter(myLogFormatter(base=10.,labelOnlyBase=False))
    for ii in range(1,ndfs):
	    print "Working on DF %i ..." % ii
	    if not precalc:
		    oortA= numpy.zeros(len(so))
		    oortB= numpy.zeros(len(so))
		    for jj in range(len(so)):
			    sys.stdout.write('\r'+"Working on %i out of %i" % \
					     (jj+1,len(so)))
			    sys.stdout.flush()
			    oortA[jj]= .5-dfs[ii][jj].oortA(1.)
			    oortB[jj]= .5+dfs[ii][jj].oortB(1.)
		    oortAs.append(oortA)
		    oortBs.append(oortB)
		    sys.stdout.write('\n')
		    sys.stdout.flush()
	    va= vas[ii]
	    if not justCalc:
		    Apos= numpy.all((oortAs[ii]-betas[ii]/2.) > 0.)
		    Aneg= numpy.all((oortAs[ii]-betas[ii]/2.) < 0.)
		    Bpos= numpy.all((oortBs[ii]+betas[ii]/2.) > 0.)
		    Bneg= numpy.all((oortBs[ii]+betas[ii]/2.) < 0.)
		    if Apos: labels[ii]+= '-,'
		    elif Aneg: labels[ii]+= '+,'
		    else: labels[ii]+= '?,'
		    if Bpos: labels[ii]+= '+)$'
		    elif Bneg: labels[ii]+= '-)$'
		    else: labels[ii]+= '?)$'
		    lines.append(bovy_plot.bovy_plot(va,
						     do_fabs(oortAs[ii]-betas[ii]/2.),
						     loglog=True,
						     overplot=True,
						     ls='none',marker=symbols[ii],color=colors[ii],mfc='none')[0])
		    if betas[ii] == 0.2: print do_fabs(oortBs[ii]+betas[ii]/2.)
		    if betas[ii] == 0.1: print do_fabs(oortBs[ii]+betas[ii]/2.)
		    bovy_plot.bovy_plot(va,
					do_fabs(oortBs[ii]+betas[ii]/2.),
					loglog=True,overplot=True,
					ls='none',marker=symbols[ii],mec='k',mfc='none')
		    if not noTheory:
			    bovy_plot.bovy_plot(xrange,
						(1.+2./hss[ii]-ks[ii]*(1./hrs[ii]+2./hss[ii])/2.-betas[ii])/2.*numpy.array(xrange),
						'-',color='0.5',overplot=True)
			    bovy_plot.bovy_plot(xrange,
						-(-1.+2./hss[ii]-ks[ii]*(1./hrs[ii]+2./hss[ii])/2.-betas[ii])/2.*numpy.array(xrange),
					'-',color='0.5',overplot=True)
		    size= 15
		    bovy_plot.bovy_text(0.00011,0.0003,r'$|\Delta A|$',fontsize=size)
		    bovy_plot.bovy_text(0.00011,0.00003,r'$|\Delta B|$',fontsize=size)
		    #Add legend
		    l1= pyplot.legend(lines[0:5],labels[0:5],loc=2,
				      frameon=False,numpoints=1)
		    l2= pyplot.legend(lines[5:ndfs],labels[5:ndfs],loc=4,
				      frameon=False,numpoints=1)
		    pyplot.gca().add_artist(l1)
    if not justCalc:
	    bovy_plot.bovy_end_print(plotfilename)
    if not precalc:
	    oortfile= open(oortfilename,'wb')
	    pickle.dump(oortAs,oortfile)
	    pickle.dump(oortBs,oortfile)
	    oortfile.close()

def do_fabs(x):
	return numpy.fabs(x)

if __name__ == '__main__':
	vaOort(sys.argv[1])
