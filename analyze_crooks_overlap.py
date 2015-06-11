#!/usr/bin/env python
# pymacs  Copyright Notice
# ============================
#
# The pymacs source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pymacs is Copyright (C) 2006-2011 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

import sys, os, time
from copy import deepcopy
from pymacs_beta import *
from pymacs_beta.parser import *
from pylab import *
from scipy.integrate import simps
from scipy.optimize import fmin, leastsq
from scipy.special import erf
from random import gauss, randint, uniform
import numpy
import warnings

debug = True

params = {#'backend': 'ps',
#          'axes.labelsize': 10,
#          'text.fontsize': 10,
          'legend.fontsize': 12,
#          'xtick.labelsize': 8,
#          'ytick.labelsize': 8,
#          'text.usetex': True,
}#          'figure.figsize': fig_size}
rcParams.update(params)


def tee( fp, s ):
    print >>fp, s
    print s
    
def cgi_error_from_mean(nruns, mu1, sig1, n1, mu2, sig2, n2):
    iseq = []

    for k in range(nruns):
        g1 = []
        g2 = []
        for i in range(n1):
            g1.append( gauss(mu1, sig1))
        for i in range(n2):
            g2.append( gauss(mu2, sig2))
        m1 = average(g1)
        s1 = std(g1)
        m2 = average(g2)
        s2 = std(g2)
        p1 = 1./(s1*sqrt(2*pi))
        p2 = 1./(s2*sqrt(2*pi))
        iq = (m1+m2)/2.
        iseq.append(iq)
    mean = average(iseq)
    err = std(iseq)
    return err

def cgi_error(nruns, mu1, sig1, n1, mu2, sig2, n2):
    iseq = []
    for k in range(nruns):
        g1 = []
        g2 = []
        for i in range(n1):
            g1.append( gauss(mu1, sig1))
        for i in range(n2):
            g2.append( gauss(mu2, sig2))
        m1 = average(g1)
        s1 = std(g1)
        m2 = average(g2)
        s2 = std(g2)
        p1 = 1./(s1*sqrt(2*pi))
        p2 = 1./(s2*sqrt(2*pi))
        iq = gauss_intersection([p1,m1,s1],[p2,m2,s2])
        iseq.append(iq)
    mean = average(iseq)
    err = std(iseq)
    return err


def sort_file_list( lst ):

    # we assume that the directory is numbered
    # guess directory base name first
    dir_name = lst[0].split('/')[-2]
    base_name = ''
    for i, x in enumerate(dir_name):
        if x.isdigit():
            check = True
            for k in range(i, len(dir_name)):
                if not dir_name[k].isdigit():
                    check = False
            if check:
                base_name = dir_name[:i]
                break
    if base_name:
        get_num = lambda s: int(s.split('/')[-2].split(base_name)[1])
        lst.sort( lambda a, b: cmp(get_num(a), get_num(b)) )
        return lst
    else:
        return lst

def process_dgdl( fn, ndata = -1, lambda0 = 0 ):
    sys.stdout.write('\r------>  %s' % fn)
    sys.stdout.flush()
    l = open(fn).readlines()
    if not l: return None, None
    r = []
    for line in l:
        if line[0] not in '#@&':
            try:
                r.append( [ float(x) for x in line.split() ] )
            except:
                print ' !! Skipping %s ' % (fn )
                return None, None
                
    if ndata != -1 and len(r) != ndata:
        try:
            print ' !! Skipping %s ( read %d data points, should be %d )' % (fn, len(r), ndata )
        except:
            print ' !! Skipping %s ' % (fn )
        return None, None
    # convert time to lambda
    ndata = len( r )
    dlambda = 1./ float( ndata )
    if lambda0 == 1: dlambda*=-1
#    if debug:
#        print 'dlambda = ', dlambda
    data = []

    for i, (time, dgdl) in enumerate(r):
        data.append( [ lambda0+i*dlambda, dgdl] )
    x = map( lambda a: a[0], data )
    y = map( lambda a: a[1], data )
    if lambda0 == 1:
        x.reverse()
        y.reverse()
    return simps( y, x ), ndata

def check_first_dgdl( fn, lambda0 ):

    l = open(fn).readlines()
    if not l: return None
    r = []
    for line in l:
        if line[0] not in '#@&':
            r.append( [ float(x) for x in line.split() ] )
    ndata = len( r )
    dlambda = 1./ float( ndata )
    if lambda0 == 1: dlambda*=-1
    print '---------------------------------------------'
    print '\t\t Checking simulation data.....'
    print '\t\t File: %s' % fn
    print '\t\t # data points: %d' % ndata
    print '\t\t Length of trajectory: %8.3f ps' % r[-1][0]
    print '\t\t Delta lambda: %8.5f' % dlambda
    print '---------------------------------------------'
    
def work_from_crooks( lst, lambda0 ):
    print '\nProcessing simulation data......'
    output_data = []
    check_first_dgdl( lst[0], lambda0 )
    first_res, ndata = process_dgdl( lst[0], lambda0 = lambda0 )
    output_data.append( [ lst[0], first_res] )
    results = [ first_res ]
    for f in lst[1:]:
        res, tmp = process_dgdl( f, ndata = ndata, lambda0 = lambda0 )
        if res is not None:
            results.append( res )
            output_data.append( [ f, res] )
    print
    return results, output_data

def data_to_gauss( data ):
    m = average( data )
    dev = std( data )
    A = 1./(dev*sqrt(2*pi))
    return m, dev, A

def gauss_intersection( g1, g2 ):
    A1, m1, s1 = g1
    A2, m2, s2 = g2
    p1 = m1/s1**2-m2/s2**2
    p2 = sqrt(1/(s1**2*s2**2)*(m1-m2)**2+2*(1/s1**2-1/s2**2)*log(s2/s1))
    p3 = 1/s1**2-1/s2**2
    x1 = (p1+p2)/p3
    x2 = (p1-p2)/p3
    # determine which solution to take
    if x1 > m1 and x1 < m2 or \
       x1 > m2 and x1 < m1:
        return x1
    elif x2 > m1 and x2 < m2 or \
       x2 > m2 and x2 < m1:
        return x2
    else:
        return False # we do not take the intersection

def ksref():
    
    f = 1
    potent = 10000
    lamb = arange(0.25,2.5,.001)
    q=array(zeros(len(lamb),float))
    res = []
    for k in range(-potent,potent):
        q=q+f*exp(-2.0*(k**2)*(lamb**2))
        f=-f
    for i in range(len(lamb)):
        res.append((lamb[i],q[i]))
    return res

def ksfunc(lamb):
    f = 1
    potent = 10000
    q=0
    for k in range(-potent,potent):
        q=q+f*exp(-2.0*(k**2)*(lamb**2))
        f*=-1
    return q


def ks(data, alpha=.05, refks = None):
    N = len(data)
    nd, ed = edf(data)
    cd = cdf(data)
    siglev = 1-alpha
    dval=[]
    for i, val in enumerate(ed):
        d = abs(val-cd[i])
        dval.append(d)
        if i:
            d = abs(ed[i-1]-cd[i])
            dval.append(d)
    dmax=max(dval)
    check = math.sqrt(N)*dmax
    if not refks:
        refks = ksref()
    lst = filter(lambda x: x[1] > siglev, refks)
    lam0 = lst[0][0]
    if check >= lam0:
        bOk = False
    else:
        bOk = True

    q = ksfunc(check)
    return (1-q), lam0, check, bOk

def edf( dg_data ):
    edf_=[]
    ndata=[]
    data = deepcopy( dg_data )
    data.sort()
    N=float(len(data))
    cnt=0
    for item in data:
        cnt+=1
        edf_.append(cnt/N)
        ndata.append(item)
    ndata=array(ndata)
    edf_=array(edf_)
    return ndata,edf_

def cdf( dg_data ):
    data = deepcopy( dg_data )
    data.sort()
    mean = average(data)
    sig = std(data)
    cdf=0.5*(1+erf((data-mean)/float(sig*sqrt(2))))
    return cdf


def data_from_file( fn ):
    data =  read_and_format( fn ,'sf')
    return map( lambda a: a[1], data)

def dump_integ_file( fn, data):
    fp = open(fn,'w')
    for fn, w in data:
        print >>fp, fn, w
    fp.close()
    

def BAR(res_ab, res_ba, T = 298):
    kb=0.00831447215
    beta = 1./(kb*T)
    
    nf = float(len(res_ab))
    nr = float(len(res_ba))
    M = kb*T*log(nf/nr)
    
    res_ab = array(res_ab)
    res_ba = array(res_ba)
    
    def func(x, res_ab, res_ba):
        sf = 0
        for v in res_ab:
            sf+=1./(1+exp(beta*(M+v - x)))
        sr = 0
        for v in res_ba:
            sr+=1./(1+exp(-beta*(M+v - x)))
        sf/=nf
        sr/=nr
        r = sf-sr
        return r**2

    avA = average(res_ab)
    avB = average(res_ba)
    x0 = (avA+avB)/2.
    result=fmin(func,x0 = x0, args = (res_ab, res_ba))
    return result

def BAR_err(dG, res_ab, res_ba, T = 298):
    kb=0.00831447215
    beta = 1./(kb*T)
    res_ab = array(res_ab)
    res_ba = array(res_ba)
    nf = float(len(res_ab))
    nr = float(len(res_ba))
    M = kb*T*log(nf/nr)
    err = 0
    for v in res_ab:
       try:
          err+=  1./(2+2*cosh(beta*(M+v-dG)))
       except RuntimeWarning, e:
          print e        
    for v in res_ba:
       try:
          err+=  1./(2+2*cosh(beta*(M+v-dG)))
       except RuntimeWarning, e:
          print e
    N = nf+nr
    err/=float(N)
    tot = 1/(beta**2*N)*(1./err-(N/nf+N/nr))
    return sqrt(tot)

def gauss_func( A, mean, dev, x):
    x = array(x)
    y = A*exp(-(((x-mean)**2)/(2.0*(dev**2))))
    return y

def dG_from_overlap(data1, data2, T=298, fname=None, dpi=300, bins=50):
    kb=0.00831447215
    slope = 1./(kb*T)
    warnings.simplefilter("error")

    # make histograms 
    maxi = max( data1+data2 )
    mini = min( data1+data2 )    
    w1=numpy.ones(len(data1))
    w1=w1*(1./len(data1))
    nf, bins1 = numpy.histogram(data1, bins=bins, range = (mini, maxi), weights=w1)
    w2=numpy.ones(len(data2))
    w2=w2*(1./len(data2))
    nb, bins2 = numpy.histogram(data2, bins=bins, range = (mini, maxi), weights=w2)  

    #create x,y data y = log(Pf(x) / Pb(x))    
    x = []
    y = []
    #print numpy.sum(nf * numpy.diff(bins1)) 
    for i in range(len(nf)):
       if ( (nf[i] > 0) and (nb[i] > 0) ):
	       y.append( log(nf[i]/nb[i]) )
	       x.append( bins1[i] )
          #print nf[i], nb[i]
    
    #linear fit with imposed known slope
    xarr = numpy.array(x)
    yarr = numpy.array(y)
    if (len(x) > 1):
       try:
          # Fit to y=ax+b when "a" is known. 
          fitfunc = lambda p, xs : slope*xs+p[0]              #Target function
          errfunc = lambda p, xs, ys: fitfunc(p, xs) - ys     # Distance to the target function
          #initial guess list, only one param to guess here
          p0 = [1.0]         
          #actual fit
          p1, success = leastsq(errfunc, p0[:], args=(xarr, yarr))
          
          #fit = array([a,b]) (ax+b), without slope constraint
          fit = polyfit(x,y,1)  
                    
          if (fname is not None):
             #time = linspace(xarr.min(), xarr.max(), len(x))                                         
             fit_fn = poly1d(fit)  #for plotting
             figure( figsize = (8, 6) )
             plot(x,y, 'bo', label = "W vs. log(Pf/Pb)")
             plot(x, fit_fn(x), '--r',  label ="y=%2.3f*x + (%4.1f)" %(fit[0], fit[1]))      
             plot(x, fitfunc(p1, xarr), "r-", label="slope imposed: %3.3f" %slope) 
             xlabel('W [kJ/mol]', fontsize=20)
             ylabel('log(Pf/Pb)', fontsize=20)
             title('Work vs. log of the overlapping histograms')
             legend(loc='best', fancybox = True)
             grid(lw = 1)      
             savefig( fname, dpi= dpi )
          #return  (#imposed-beta dG, beta, overlapping points, unimposed-beta dG)	
          # -fit[1]/fit[0] = dG with unimposed beta. -p1[0]/slope = dG with imposed beta
          try:
             return -p1[0]/slope, fit[0], len(x), -fit[1]/fit[0]
          except RuntimeWarning, e:
             print e
             return 0,0,0,0
          # y=ax+b, dG is intersection with x (so y=0, x= -b/a )
       except RankWarning, e:
          print e
          return 0,0,0,0
         
    else:
       print "No overlap, cannot calculate dG with linear regression"
       return 0,0,0,0


#calculates error from the overlap analysis, assuming fwd,back distribute like gaussians
def overlap_error_g(nruns, mu1, sig1, n1, mu2, sig2, n2,T):
    list_of_dG = []
    list_of_beta = []
    for k in range(nruns):
        g1 = []
        g2 = []
        for i in range(n1):
            g1.append( gauss(mu1, sig1))
        for i in range(n2):
            g2.append( gauss(mu2, sig2))
        iq, beta, N, unimpQ = dG_from_overlap(g1,g2,T)        
        if (iq != 0 and beta != 0 and N!=0):
           list_of_dG.append(iq)
           list_of_beta.append(beta)
    mean = average(list_of_dG)
    meanb = average(list_of_beta)    
    err = std(list_of_dG)
#    print mean, err, max(list_of_dG), min(list_of_dG)
    errb = std(list_of_beta)
    return err, errb

#calculating error from overlap analysis, no gaussians assumption
def overlap_error_hists(nruns, data1, data2, T, bins=50):
   kb=0.00831447215
   slope = 1./(kb*T)
   
   warnings.simplefilter("error")
   list_of_dG = []
   list_of_beta = []
   maxi = max(data1+data2)
   mini = min(data1+data2)
   # get the histograms from the data, weighted so that cumulative sum = 1
   w1=numpy.ones(len(data1))
   w1=w1*(1./len(data1))
   w2=numpy.ones(len(data2))
   w2=w2*(1./len(data2))
   P1, E1 = numpy.histogram(data1, bins=bins, range=(mini, maxi), weights=w1)
   P2, E2 = numpy.histogram(data2, bins=bins, range=(mini, maxi), weights=w2)
   #print P1, len(P1)
   #print E1
   #print P2, len(P2)
   figure( figsize = (8, 6) )
   xlabel('W [kJ/mol]', fontsize=20)
   ylabel('log(Pf/Pb)', fontsize=20)
   title('Work vs. log - Bootstrapping histograms')
   grid(lw = 1)
   
   #-Generate new histograms for $nruns times. 
   # 1. zero out historams of the same size (bins)
   # 2. generate random numbers according to P1 and P2 distributions for the Forward and Backward histograms
   # 3. this is done according to the length of data1 (F) and data2(B). 
   #-Compute linear fit ($nruns times).
   #-Compute mean and std of linear fits distribution
   for n in range(nruns):  
      newF = numpy.zeros(bins)
      newB = numpy.zeros(bins)
      for j in range(len(data1)):
         newF[random_from_weighted_dist(P1)] += 1
      for j in range(len(data2)):
         newB[random_from_weighted_dist(P2)] += 1      
      #calc dG from overlapping part of new histograms
      x=[]
      y=[]
      count=0
      for i in range(len(newF)):
         if (newF[i] > 0) and (newB[i] > 0):
            lg = log(newF[i]/newB[i])
            count += 1
            if not (lg != lg): #check for NaN, though it will not happen here.
               x.append(E1[i])
               y.append(lg)            
      #print x
      #print y
      xarr = numpy.array(x)
      yarr = numpy.array(y)
      try:
         if (len(x) > 1):		#more than 1 point for linear regression
            fit = polyfit(x,y,1)   #old fit without slope constraint
            #fit_fn = poly1d(fit)  #for plotting
            # Fit to y=ax+b when "a" is known. 
            fitfunc = lambda p, xs : slope*xs+p[0]              #Target function
            errfunc = lambda p, xs, ys: fitfunc(p, xs) - ys     #Dist from target func
        
            #initial guess list, only one here
            p0 = [1.0]          
            #actual fit
            p1, success = leastsq(errfunc, p0[:], args=(xarr, yarr))
            
            tmpdG =  -p1[0]/slope
            if not (tmpdG != tmpdG):         
               list_of_dG.append(tmpdG )
               plot(x, fitfunc(p1, xarr), "r-") 
               #this list of betas is from unimposed slope fits (old calculation)
               #list_of_beta.append( fit[0] )
               #plot(x, fit_fn(x), 'r' )   old fit
               
    
      except Exception, e:
         print e
         pass
         
   #print list_of_dG
     
   savefig( "bootstrap_reg.png", dpi= 300 )      
   mean = average(list_of_dG)
   err = std(list_of_dG)   
   print mean, err
   return err

#returns the bin in the weighted distribution that the new random belongs to
def random_from_weighted_dist(P):
   r = uniform(0,1)
   csum = 0
   for i in range(len(P)):
      csum += P[i]
      if (r < csum):
         return i
   return i
   
    	
def make_plot( fname, data1, data2, result, err, nbins, dpi ):

    figure( figsize = (8, 6) )
    mf, devf, Af = data_to_gauss( data1 )
    mb, devb, Ab = data_to_gauss( data2 )
    
    maxi = max( data1+data2 )
    mini = min( data1+data2 )
    n1, bins1, patches1 = hist(data1, range = (mini,maxi),bins=nbins, facecolor='blue', alpha=0.75, normed=True, label='0->1')
    n2, bins2, patches2 = hist(data2, range = (mini,maxi),bins=nbins, facecolor='red', alpha=0.75, normed=True, label='1->0')
  
    xlabel('W [kJ/mol]', fontsize=20)
    ylabel('Probability', fontsize=20)
    title('Work Distribution $\lambda$ 0->1 (blue) $\lambda$ 1->0 (red)')
    grid(lw = 2)
    loc, lab = yticks()
    ll = []
    for i in range(len(lab)):
        ll.append("")
    yticks( loc, ll )
    x = arange( mini, maxi, .5 )
    y1 = gauss_func( Af, mf, devf, x )
    y2 = gauss_func( Ab, mb, devb, x )
    
    plot(x, y1, 'b--', linewidth=2)
    plot(x, y2, 'r--', linewidth=2)
    
    size = max( [max(y1), max(y2)] )
    res_x = [result, result ]
    res_y = [0, size*1.2 ]
    plot( res_x, res_y, 'k--', linewidth=2, label = r'$\Delta$G = %.2f $\pm$ %.2f kJ/mol' % (result, err))
    legend(shadow=True, fancybox = True)
    ylim(0, size*1.2 )
    xl = gca()
    for val in xl.spines.values():
        val.set_lw(2)
    savefig( fname, dpi= dpi )
        

def make_W_over_time_plot( fname, data1, data2, result, err, nbins, dpi):

    figure( figsize = (8, 6) )
    x1 = range( len(data1) )
    x2 = range( len(data2) )
    if x1>x2: x = x1
    else: x = x2
    mf, devf, Af = data_to_gauss( data1 )
    mb, devb, Ab = data_to_gauss( data2 )
    
    maxi = max( data1+data2 )
    mini = min( data1+data2 )

    sm1 = smooth( array(data1) )
    sm2 = smooth( array(data2) )
    subplot(1,2,1)
    plot(x1,data1,'g-',linewidth=2 ,label="Forward (0->1)", alpha=.3)
    plot(x1,sm1,'g-',linewidth=3) 
    plot(x2,data2,'b-',linewidth=2 ,label="Backward (1->0)", alpha=.3)
    plot(x2,sm2,'b-',linewidth=3) 
    legend(shadow=True, fancybox = True, loc='upper center')
    ylabel(r'W [kJ/mol]', fontsize = 20)
    xlabel(r'# Snapshot', fontsize = 20)
    grid(lw=2)
    xlim(0,x[-1]+1)
    xl = gca()
    for val in xl.spines.values():
        val.set_lw(2)
    subplot(1,2,2)
    hist(data1,bins=nbins, orientation='horizontal', facecolor='green',alpha=.75, normed=True)
    hist(data2,bins=nbins, orientation='horizontal', facecolor='blue',alpha=.75, normed=True)

    x = arange( mini, maxi, .5 )

    y1 = gauss_func( Af, mf, devf, x )
    y2 = gauss_func( Ab, mb, devb, x )

    plot(y1, x, 'g--', linewidth=2)
    plot(y2, x, 'b--', linewidth=2)
    size = max( [max(y1), max(y2)] )
    res_x = [result, result ]
    res_y = [0, size*1.2 ]
    plot( res_y, res_x, 'k--', linewidth=2, label = r'$\Delta$G = %.2f $\pm$ %.2f kJ/mol' % (result, err))
    legend(shadow=True, fancybox = True, loc='upper center')
    xticks([])
    yticks([])
    xl = gca()
    for val in xl.spines.values():
        val.set_lw(2)
    subplots_adjust(wspace=0.0, hspace = 0.1)
    savefig(fname, dpi=dpi)
    
def smooth(x,window_len=11,window='hanning'):

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval(window+'(window_len)')
    y=convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def select_random_subset( lst, n):
    ret = []
    idx = []
    while len(ret) < n:
        rn = randint(0, len(lst)-1)
        if rn not in idx:
            idx.append( rn )
            ret.append( lst[rn])
    print idx
    return ret

def main(argv):

    version = "1.1"

    options = [
        Option( "-nbins", "int", 20, "number of histograms bins for plot"),
        Option( "-obins", "int", 50, "number of histograms bins for overlap dG calculation"),
        Option( "-T", "real", 298, "Temperature for BAR calculation"),
        Option( "-dpi", "int", 300, "plot resolution"),
        Option( "-reverseB", "bool", False, "reverse state B"),
        Option( "-firstA", "int", 0, "first trajectory to analyze (by default all values are taken)"),
        Option( "-lastA", "int", 100, "last trajectory to analyze (by default all values are taken)"),
        Option( "-firstB", "int", 0, "first trajectory to analyze (by default all values are taken)"),
        Option( "-lastB", "int", 100, "last trajectory to analyze (by default all values are taken)"),
        Option( "-rand", "int", 50, "take a random subset of trajectories"),
        Option( "-integ_only", "bool", False, "Do integration only. Skip analysis."),
        Option( "-KS", "bool", True, "Do Kolmogorov-Smirnov test"),
        ]
    
    files = [
        FileOption("-pa", "r/m",["xvg"],"dgdl.xvg", "paths to 0->1 runs"),
        FileOption("-pb", "r/m",["xvg"],"dgdl.xvg", "paths to 1->0 runs"),
        FileOption("-o", "w",["dat"],"results.dat", "results"),
        FileOption("-cgi_plot", "w",["png","eps","svg","pdf"],"cgi.png", "plot work histograms "),
        FileOption("-W_over_t", "w",["png","eps","svg","pdf"],"W_over_t.png", "plot work over time "),
        FileOption("-i0", "r/m/o",["dat"],"integ0.dat", "read integrated W (0->1)"),
        FileOption("-i1", "r/m/o",["dat"],"integ1.dat", "read integrated W (1->0)"),
        FileOption("-o0", "w",["dat"],"integ0.dat", "write integrated W (0->1)"),
        FileOption("-o1", "w",["dat"],"integ1.dat", "write integrated W (1->0)"),
        
        ]
    
    
    
    help_text = ('Calculates free energies from fast growth  ',
                 'thermodynamic integration runs.',
                 'First method: Crooks-Gaussian Intersection (CGI)',
                 'Second method: Benett Acceptance Ratio (BAR)'
                 )

    
    cmdl = Commandline( argv, options = options,
                        fileoptions = files,
                        program_desc = help_text,
                        check_for_existing_files = False, version = version)

    out = open(cmdl['-o'],'w')
    print >>out, "# analyze_crooks.py, version = %s" % version
    print >>out, "# pwd = %s" % os.getcwd()
    print >>out, "# %s (%s)" % (time.asctime(), os.environ.get('USER') )
    print >>out, "# command = %s" % ' '.join(argv)
    print >>out, "#------------------------------------------------"
    
    if not cmdl.opt['-i0'].is_set: 
        run_ab = cmdl['-pa']
        run_ba = cmdl['-pb']
        run_ab = sort_file_list( run_ab )
        run_ba = sort_file_list( run_ba )
        res_ab, ab_data = work_from_crooks( run_ab, lambda0 = 0 )
        res_ba, ba_data = work_from_crooks( run_ba, lambda0 = 1 )
        dump_integ_file( cmdl['-o0'], ab_data)
        dump_integ_file( cmdl['-o1'], ba_data)
    else:
        res_ab = []
        res_ba = []
        for fn in cmdl['-i0']:
            print '\t\tReading integrated values (0->1) from', fn
            res_ab.extend(data_from_file( fn ) )
        for fn in cmdl['-i1']:
            print '\t\tReading integrated values (1->0) from', fn
            res_ba.extend(data_from_file( fn ) )

    if cmdl['-integ_only']:
        print '\n    Integration done. Skipping analysis.'
        print '\n    ......done........\n'
        sys.exit(0) 

    firstA = 0
    lastA = len(res_ab)
    firstB = 0
    lastB = len(res_ba)
    if cmdl.opt['-firstA'].is_set:
        firstA = cmdl['-firstA']
        tee(out, '   first trajectory to read from A: %d' % firstA)
    if cmdl.opt['-lastA'].is_set:
        lastA = cmdl['-lastA']
        tee(out, '   last trajectory to read from A : %d' % lastA)
    if cmdl.opt['-firstB'].is_set:
        firstB = cmdl['-firstB']
        tee(out, '   first trajectory to read from B: %d' % firstB)
    if cmdl.opt['-lastB'].is_set:
        lastB = cmdl['-lastB']
        tee(out, '   last trajectory to read from B : %d' % lastB)

    res_ab = res_ab[firstA:lastA]
    res_ba = res_ba[firstB:lastB]


    if cmdl.opt['-rand'].is_set:
        ntraj = cmdl['-rand']
        tee(out, ' select random subset of trajectories: %d' % ntraj )
        res_ab = select_random_subset(res_ab, ntraj)
        res_ba = select_random_subset(res_ba, ntraj)
        
    mf, devf, Af = data_to_gauss( res_ab )
    mb, devb, Ab = data_to_gauss( res_ba )
    tee(out, ' --------------------------------------------------------')
    tee(out, '             ANALYSIS: NUMBER OF TRAJECTORIES:')
    tee(out, '               0->1 : %d' % len(res_ab))
    tee(out, '               1->0 : %d' % len(res_ba))
    
    tee(out, ' --------------------------------------------------------')
    tee(out, '             ANALYSIS: Crooks-Gaussian Intersection     ')
    tee(out, ' --------------------------------------------------------')
    tee(out, '  Forward  : mean = %8.3f  std = %8.3f' % ( mf, devf ))
    tee(out, '  Backward : mean = %8.3f  std = %8.3f' % ( mb, devb ))

    if cmdl['-KS']:
        tee(out, '  Running KS-test ....')
        q0, lam00, check0, bOk0 = ks(res_ab)
        q1, lam01, check1, bOk1 = ks(res_ba)
    
        tee(out, '  Forward  : gaussian quality = %3.2f' % q0)
        if bOk0:
            tee(out, '             ---> KS-Test Ok')
        else: 
            tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f, lambda0 = %4.2f' %( q0, check0 ))
        tee(out, '  Backward : gaussian quality = %3.2f' % q1)
        if bOk1:
            tee(out, '             ---> KS-Test Ok')
        else: 
            tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f, lambda1 = %4.2f' %( q1, check1 ))
       
   

    tee(out, '  Calculating Intersection...')
    cgi_result = gauss_intersection( [Af, mf, devf], [Ab, mb, devb ] )
    intersection = True
    T = cmdl['-T']
    if not cgi_result:
        tee(out, '\n  Gaussians to close for intersection calculation')
        tee(out, '   --> Taking difference of mean values')
        cgi_result = (mf+mb)*.5
        intersection = False
    tee(out, '  RESULT: dG ( CGI )  = %8.4f kJ/mol' % cgi_result)
    
    if intersection:
        cgi_err = cgi_error( 1000, mf, devf, len( res_ab), mb, devb, len(res_ba ) )
    else:
        cgi_err = cgi_error_from_mean( 1000, mf, devf, len( res_ab), mb, devb, len(res_ba ) )
    tee(out, '  RESULT: error_dG ( CGI ) = %8.4f kJ/mol' % cgi_err)
   
    # overlap result: added by Hadas
    tee(out, ' --------------------------------------------------------')
    tee(out, '             ANALYSIS: estimation from Overlap           ')
    tee(out, ' --------------------------------------------------------') 
    overlap_dG_result, oBeta, oN, unimp_dG = dG_from_overlap(res_ab, res_ba, T, "overlap_reg.png", cmdl['-dpi'], bins=cmdl['-obins'])
    overlap_err  = 0
    overlap_err2 = 0
    beta_err     = 0

    if (oN > 0):
       tee(out, '  RESULT: dG ( overlap )  = %8.4f kJ/mol' % overlap_dG_result)
       tee(out, '  RESULT: dG ( overlap unimposed beta )  = %8.4f kJ/mol' % unimp_dG)
       kb=0.00831447215
       tee(out, '  RESULT: Unimposed Beta vs. Real beta  = %3.3f %3.3f' % (oBeta, 1./(kb*T)) )
       
       overlap_err, beta_err = overlap_error_g(1000, mf, devf, len( res_ab), mb, devb, len(res_ba ),T)       
       overlap_err2 = overlap_error_hists(1000, res_ab, res_ba, T, bins=cmdl['-obins'])
 
       tee(out, '  RESULT: Number of overlapping points = %d '  % oN)
       tee(out, '  RESULT: error_dG (overlap_Gs) = %8.4f Kj/mol  error_beta: = %3.3f'  % (overlap_err, beta_err))
       tee(out, '  RESULT: error_dG (overlap) = %8.4f Kj/mol '  % overlap_err2)
    #check for huge errors resulting from too few overlap points or NaN
    if (abs(overlap_err) > 1000 or (overlap_err != overlap_err)): 
       overlap_err = 0
    if (abs(overlap_err2) > 1000 or (overlap_err2 != overlap_err2)): 
       overlap_err2 = 0 
    
    ##
    tee(out, ' --------------------------------------------------------')
    tee(out, '             ANALYSIS: Bennett Acceptance Ratio     ')
    tee(out, ' --------------------------------------------------------')
    
    tee(out, '  Solving numerical equation with Nelder-Mead Simplex algorithm.. ')
    tee(out, '  Temperature used: %8.2f K' % T)
    bar_result = BAR( res_ab, res_ba, T)
    tee(out, '  RESULT: dG (BAR ) = %8.4f kJ/mol' % bar_result)
    bar_err = BAR_err( bar_result, res_ab, res_ba, T)
    
    tee(out, '  RESULT: error_dG (BAR ) = %8.4f kJ/mol' % bar_err)
    tee(out, ' ------------------------------------------------------')
    diff = abs( cgi_result - bar_result )
    mean = (cgi_result+bar_result)*.5
    tee(out, '  Difference between BAR and CGI = %8.5f kJ/mol' % diff ) 
    tee(out, '  Mean of  BAR and CGI           = %8.5f kJ/mol' % mean )
    # CGI errCGI BAR errBAR overlap errOVERLAP err2_OVERLAP beta beta_err
    tee(out, ' HUMMER dG imposed-beta vs. unimposed-beta: %8.4f\t%8.4f' %(overlap_dG_result, unimp_dG))
    tee(out, ' REMARK: CGI\t    errCGI \tBAR\t errBAR \tHummer\t errHummer err2_Hummer \tbeta\toverlap points\tHummer unimposed-beta')
   
    tee(out, ' FINAL: %8.4f  %8.4f \t%8.4f %8.4f \t %8.4f %8.4f %8.4f \t %3.3f \t%d' %(cgi_result, cgi_err, bar_result, bar_err, overlap_dG_result, overlap_err, overlap_err2, oBeta, oN ))
    tee(out, ' ------------------------------------------------------')
    

    print '\n   Plotting histograms......'
    make_plot( cmdl['-cgi_plot'], res_ab, res_ba, cgi_result, cgi_err, cmdl['-nbins'], cmdl['-dpi'] )
    make_W_over_time_plot( cmdl['-W_over_t'], res_ab, res_ba, cgi_result, cgi_err, cmdl['-nbins'], cmdl['-dpi'])

    tee(out, '\n   ......done...........\n')
    
   
main( sys.argv )
