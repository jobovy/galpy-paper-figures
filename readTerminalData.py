import numpy
def readClemens(dsinl=0.5/8.):
    data= numpy.loadtxt('mwpot14data/clemens1985_table2.dat',delimiter='|',
                        comments='#')
    glon= data[:,0]
    vterm= data[:,1]
    #Remove l < 30 and l > 80
    indx= (glon > 40.)*(glon < 80.)
    glon= glon[indx]
    vterm= vterm[indx]
    if bin:
        #Bin in l=1 bins
        glon, vterm= binlbins(glon,vterm,dl=1.)
        #Remove nan, because 1 bin is empty
        indx= True-numpy.isnan(glon)
        glon= glon[indx]
        vterm= vterm[indx]
    #Calculate correlation matrix
    singlon= numpy.sin(glon/180.*numpy.pi)
    corr= calc_corr(singlon,dsinl)
    return (glon,vterm,numpy.linalg.inv(corr))

def readMcClureGriffiths(dsinl=0.5/8.,bin=True):
    data= numpy.loadtxt('mwpot14data/McClureGriffiths2007.dat',
                        comments='#')
    glon= data[:,0]
    vterm= data[:,1]
    #Remove l > 330 and l > 80
    indx= (glon < 320.)*(glon > 280.)
    glon= glon[indx]
    vterm= vterm[indx]
    if bin:
        #Bin in l=1 bins
        glon, vterm= binlbins(glon,vterm,dl=1.)
    #Calculate correlation matrix
    singlon= numpy.sin(glon/180.*numpy.pi)
    corr= calc_corr(singlon,dsinl)
    return (glon,vterm,numpy.linalg.inv(corr))
    
def calc_corr(singlon,dsinl):
    #Calculate correlation matrix
    corr= numpy.zeros((len(singlon),len(singlon)))
    for ii in range(len(singlon)):
        for jj in range(len(singlon)):
            corr[ii,jj]= numpy.exp(-numpy.fabs(singlon[ii]-singlon[jj])/dsinl)
    corr= 0.5*(corr+corr.T)
    return corr+10.**-10.*numpy.eye(len(singlon)) #for stability

def binlbins(glon,vterm,dl=1.):
    minglon, maxglon= numpy.floor(numpy.amin(glon)), numpy.floor(numpy.amax(glon))
    minglon, maxglon= int(minglon), int(maxglon)
    nout= maxglon-minglon+1
    glon_out= numpy.zeros(nout)
    vterm_out= numpy.zeros(nout)
    for ii in range(nout):
        indx= (glon > minglon+ii)*(glon < minglon+ii+1)
        glon_out[ii]= numpy.mean(glon[indx])
        vterm_out[ii]= numpy.mean(vterm[indx])
    return (glon_out,vterm_out)
