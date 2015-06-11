#!/home/dseelig/bin/python
############################################################
#                                                           #
#           THIS FILE IS PART OF PYMACS                     #
#   WRITTEN BY DANIEL SEELIGER (dseelig@gwdg.de)            #
#                   2006-2009                               #
#    COMPUTATIONAL BIOMOLECULAR DYNAMICS GROUP              #
#  MAX-PLANCK-INSTITUTE FOR BIOPHYSICAL CHEMISTRY           #
#                   GOETTINGEN                              #
#                                                           #
#############################################################
import sys, os
from pymacs import *
from pymacs.parser import *
from pymacs.futil import *
from numpy import *
from pylab import *
from scipy.optimize import leastsq, fsolve, fmin
from scipy.integrate import simps
from scipy.special import erf
from random import gauss, randint
import pyutil
import optparse

def error_ana(nruns, n1, mu1, sig1, n2, mu2, sig2):
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
        iq = interseq2([p1,m1,s1],[p2,m2,s2])
        iseq.append(iq)
    mean = average(iseq)
    err = std(iseq)
    return mean, err

def error_ana3(nruns, n1, mu1, sig1, n2, mu2, sig2):
    iseq = []
    m1 = mu1+sig1
    m2 = mu2+sig2
    print 'm1 = %g' % m1
    print 'm2 = %g' % m2
    p1 = 1./(sig1*sqrt(2*pi))
    p2 = 1./(sig2*sqrt(2*pi))
    iq1 = interseq([p1,m1,sig1],[p2,m2,sig2])
    peakdiff = abs(m1-m2)
    adiff=abs(abs(m1)-abs(iq1))
    bdiff=abs(abs(m2)-abs(iq1))
    if peakdiff<adiff or peakdiff<bdiff:
        print 'Switching from intersection to peak diff'
        iq1 = (m1+m2)/2.
    print 'iq1 = %g' % iq1
    m1 = mu1-sig1
    m2 = mu2-sig2
    print 'm1 = %g' % m1
    print 'm2 = %g' % m2
    iq2 = interseq([p1,m1,sig1],[p2,m2,sig2])
    peakdiff = abs(m1-m2)
    adiff=abs(abs(m1)-abs(iq2))
    bdiff=abs(abs(m2)-abs(iq2))
    if peakdiff<adiff or peakdiff<bdiff:
        print 'Switching from intersection to peak diff'
        iq2 = (m1+m2)/2.
    print 'iq2 = %g' % iq2
    print 'Estimated max error = %g' % (abs(iq1-iq2)/2.)
    return iq1, iq2

def error_ana2(nruns, n1, mu1, sig1, n2, mu2, sig2):
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
    return mean, err


def check_dgdl(fn, lda0, lda1):
    try:
        l = open(fn).readlines()
        l = kickOutComments(l,'#@')
    except:
        return 0, False, 0
    if len(l) != 0:
        entr1 = float(l[-1].split()[0])
        try:
            entr2 = float(l[-2].split()[0])
        except:
            print l
            print len(l)
            print fn
            sys.exit(1)
        timestep = entr1-entr2
        dlambda = timestep/entr1*(lda1-lda0)
#        dlambda = lda1/entr1
#        print 'dlambda = %g' % dlambda
        return l, True, dlambda
    else: return l, False, 0
    

def get_finished_runs(rundir, run='A'):
    frames = glob(rundir+'/frame*')
    dgdl = []
    for d in frames:
        if run=='A':
            if os.path.isfile(os.path.join(d,'aqoff.xvg')) and \
               os.path.isfile(os.path.join(d,'vdw.xvg')) and \
               os.path.isfile(os.path.join(d,'aqon.xvg')):
                dgdl.append([os.path.join(d,'aqoff.xvg'),os.path.join(d,'vdw.xvg'),os.path.join(d,'aqon.xvg')])
            else:
                print 'Missing files in %s' % d
        elif run=='B':
            if os.path.isfile(os.path.join(d,'bqon.xvg')) and \
               os.path.isfile(os.path.join(d,'vdw.xvg')) and \
               os.path.isfile(os.path.join(d,'bqoff.xvg')):
                dgdl.append([os.path.join(d,'bqon.xvg'),os.path.join(d,'vdw.xvg'),os.path.join(d,'bqoff.xvg')])
            else:
                print 'Missing files in %s' % d
    return dgdl


def get_dgdl(frame,dgdls):
    for p in dgdls:
        fr = p.split('/')[1]
        if frame == fr:
            return p
    print 'Error: %s not found!' % p
    sys.exit(1)

def gaussfunc(x, p):
    fac, mean, dev = p
    return fac*exp(-(((x-mean)**2)/(2.0*(dev**2))))

def residuals(p,y,x,func):
    return y-func(x,p)

def peval(x, p, func):
    return func(x,p)

def gauss_fit(x,y):
    # get initial parameters
    av = average(x)
    s = std(x)
    fac = 1
    p0 = [fac, av, s]
    lsq = leastsq(residuals,p0,args=(y,x,gaussfunc))
    newy = peval(x, lsq[0], gaussfunc)
    return lsq[0], newy
                  

def read_dgdl(lst, dl = 4e-5, lda0 = 0, asarray = True):
#    l = open(fname).readlines()
#    l = kickOutComments(l,'#@')
    l = parseList('ff',lst)
    # convert time to lambda
    x = []
    for i in range(len(l)):
        x.append( lda0+i*dl )
    y = map(lambda n: n[1], l)
    if lda0 == 1:
        x.reverse()
        y.reverse()
    if not asarray: return x, y
    else: return array(x), array(y)


def get_list(d, fname = 'dhdl.part0001.xvg', endtime=200):
   if not os.path.isdir(d):
       return []
   dl = listDirs(d)
   fl = []
   for x in dl:
      if (float(pyutil.run_command("tail -1 %s" % os.path.join(x,fname)).split()[0]) == endtime):
         fl.append(os.path.join(x,fname))
   return fl


def integrate_dgdl(lst, fname, fp, dl = 4e-5, lda0=0, bReverse = False):
    x, y = read_dgdl(lst, dl = dl, lda0=lda0)
    #res = trapz(y,x)
    res = simps(y,x)
    if bReverse:
        res*=-1
    if fp:
        print >>fp, "%g %s" % (res,fname)
    return res

def interseq(g1, g2):
    def func(x):
        r = A1*exp(-(x-m1)**2/(2*s1**2)) - A2*exp(-(x-m2)**2/(2*s2**2))
        return r
    A1, m1, s1 = g1
    A2, m2, s2 = g2
    x0 = (m1+m2)/2.
    result=fsolve(func,x0)
    return result

def interseq2(g1, g2):
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
    
def histo_plot(data1, data2, nbins, fp,show_plot=True, plname = 'plot.png'):
    maxi = max(data1+data2)
    mini = min(data1+data2)
    size = maxi-mini
#    nbins = int(size/2.)
    print data1
    print len(data1)
    print data2
    print len(data2)
    print mini, maxi, nbins
   # n,b,p=hist(data1, range=(mini, maxi), bins=30, facecolor='green', alpha=0.75, normed=True, label='0->1')
    n1, bins1, patches1 = hist(data1, range = (mini,maxi),bins=int(nbins), facecolor='green', alpha=0.75, normed=True, label='0->1')
    n2, bins2, patches2 = hist(data2, range = (mini,maxi),bins=int(nbins), facecolor='red', alpha=0.75, normed=True, label='1->0')
    xlabel('W [kJ/mol]')
    ylabel('Probability')
    title(r'Work Distribution $\lambda$ 0->1 (green) $\lambda$ 1->0 (red)')
    grid(True)
    b1 = array(map(lambda xx: xx, bins1)[:-1])
    mx = max(map(lambda xx: xx, bins1)[:-1]+map(lambda xx: xx, bins2)[:-1])
    mi = min(map(lambda xx: xx, bins1)[:-1]+map(lambda xx: xx, bins2)[:-1])
    n1 = array(map(lambda xx: xx, n1))
    b2 = array(map(lambda xx: xx, bins2)[:-1])
    n2 = array(map(lambda xx: xx, n2))


    avA = average(data1)
    avB = average(data2)
    stdA = std(data1)
    stdB = std(data2)
    par1 = [1./(stdA*sqrt(2*pi)),avA,stdA]
    par2 = [1./(stdB*sqrt(2*pi)),avB,stdB]

##     par1, newy1 = gauss_fit(b1, n1)
##     par2, newy2 = gauss_fit(b2, n2)
##     fp0 = open('histo0.dat','w')
##     print >>fp0, "# %g %g" % (avA,stdA)
##     for i in range(len(n1)):
##         print >>fp0, "%g %g" % (b1[i],n1[i])

##     fp0 = open('histo1.dat','w')
##     print >>fp0, "# %g %g" % (avB,stdB)
##     for i in range(len(n2)):
##         print >>fp0, "%g %g" % (b2[i],n2[i])

    xr = arange(mi, mx, 0.5)
    y1 = gaussfunc(xr,par1)
    y2 = gaussfunc(xr,par2)
    # check gaussion behaviour
    q0, lam00, check0, bOk0 = ks(data1)
    q1, lam01, check1, bOk1 = ks(data2)
    if bOk0:
        print 'Data 0->1 survive ks-test. Quality = %3.2f' % q0
        print >>fp, '(0->1) gaussian quality: = %3.2f (KS-test Ok)' % q0
    else:
        print 'Data 0->1 fail ks-test!!!. Quality = %3.2f; sqrt(N)*Dmax = %4.2f; lambda0 = %4.2f' % (q0,check0,lam00) 
        print >>fp, '(0-1) gaussian quality: = %3.2f (KS-test Not Ok) sqrt(N)*Dmax = %4.2f; lambda0 = %4.2f' % (q0,check0,lam00) 
    if bOk1:
        print 'Data 1->0 survive ks-test. Quality = %3.2f' % q1
        print >>fp, '(1->0) gaussian quality: = %3.2f (KS-test Ok)' % q1
    else:
        print 'Data 1->0 fail ks-test!!!. Quality = %3.2f; sqrt(N)*Dmax = %4.2f; lambda0 = %4.2f' % (q1,check1,lam01) 
        print >>fp, '(1->0) gaussian quality: = %3.2f (KS-test Not Ok) sqrt(N)*Dmax = %4.2f; lambda0 = %4.2f' % (q1,check1,lam01) 

#    av_q = average([q0,q1])
#    print 'Average gaussian quality: %3.2f' % av_q
                   

#    plot(b1, newy1, 'g--', linewidth=2, label='fit 0->1')
#    plot(b2, newy2, 'r--', linewidth=2, label='fit 1->0')
#    plot(xr, y1, 'g--', linewidth=2, label='fit 0->1')
#    plot(xr, y2, 'r--', linewidth=2, label='fit 1->0')
    plot(xr, y1, 'g--', linewidth=2)
    plot(xr, y2, 'r--', linewidth=2)
#    r1=interseq(par1,par2)
    res = interseq2(par1,par2)
    if res == None:
        print 'Gaussians are too close to calculate intersection -> Taking peak difference'
        res = (par1[1]+par2[1])/2.
        
#    print 'interseq = %g' % r1
    print 'interseq2 = %g ' % res
    peakdiff = abs(par1[1]-par2[1])
    adiff=abs(abs(par1[1])-abs(res))
    bdiff=abs(abs(par2[1])-abs(res))
    if peakdiff<adiff or peakdiff<bdiff:
        print 'Gaussians are too close to calculate intersection -> Taking peak difference'
        res = (par1[1]+par2[1])/2.
        mean, error = error_ana2(1000,len(data1), par1[1], par1[2], len(data2), par2[1], par2[2])
    else:
        mean, error = error_ana(1000,len(data1), par1[1], par1[2], len(data2), par2[1], par2[2])

#    print 'Error without ks-test correction: %4.2f' % error
#    error/=av_q
#    print 'Error after ks-test correction: %4.2f' % error
    
#    err1, err2 = error_ana3(1000,len(data1), par1[1], par1[2], len(data2), par2[1], par2[2])
#    err_diff = abs(err1-err2)
    m = max(y1+y2)
    ylim(0,m*1.2)
#    par3 = [3*m*1./(error*sqrt(2*pi)),mean,error]
    par3 = [m,res,error]
    y3 = gaussfunc(xr, par3)
    yy = [0, m]
#    yy2 = [0,m]
#    yy3 = [0,m]
#    xx2 = [err1,err1]
#    xx3 = [err2,err2]
    xx = [res, res]
    plot(xx, yy, 'k--', linewidth=2, label = r'$\Delta$G = %8.3f kj/mol' % res)
#    plot(xr, y3, 'k--', linewidth=2, label='Error = %8.3f kj/mol' % error)
#    plot(xx2, yy2, 'b--', linewidth=2)#, label='Error2 = %8.3f kj/mol' % error)
#    plot(xx3, yy3, 'b--', linewidth=2)
    legend()
    savefig(plname)
    if show_plot:
        show()
    return res,error


def read_integ_file(fn):
    l = open(fn).readlines()
    dg = map(lambda x: float(x.split()[0]), l)
#    dg.sort()
    return dg



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
#    print ' The data is gaussian distributed within a level of significance (alpha=%4.3f) of %g percent!'%(q,(1-q))

    
def edf(data):
    '''Calculates the empirical distribution function '''

    edf=[]
    ndata=[]
    data.sort()
    N=float(len(data))
    cnt=0
    for item in data:
        cnt+=1
        edf.append(cnt/N)
        ndata.append(item)
    ndata=array(ndata)
    edf=array(edf)

    return ndata,edf

def cdf(data):
    mean = average(data)
    sig = std(data)
    cdf=0.5*(1+erf((data-mean)/float(sig*sqrt(2))))
    return cdf


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
#        NA = float(len(res_ab))
#        NB = float(len(res_ba))
#        r = 1/NA * (1./(1+exp(beta*(res_ab - x)))) - 1/NB * (1./(1+exp(-beta*(res_ba - x))))
        return r**2

    avA = average(res_ab)
    avB = average(res_ba)
    x0 = (avA+avB)/2.
#    print x0
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
        err+=  1./(2+2*cosh(beta*(M+v-dG)))
    for v in res_ba:
        err+=  1./(2+2*cosh(beta*(M+v-dG)))
    N = nf+nr
    err/=float(N)
    
    tot = 1/(beta**2*N)*(1./err-(N/nf+N/nr))
    return sqrt(tot)

    

        


def main(argv):

    usage = "Usage: %prog [options]"
	
    parse = optparse.OptionParser(usage=usage)
    
    parse.add_option("-a",dest="Adir",default="AB",help="Directory with runs A->B")
    parse.add_option("-b",dest="Bdir",default="BA",help="Directory with runs B->A")
    parse.add_option("--plot",dest="showPlot",action="store_true",help="show plot")
    parse.add_option("--nbins",dest="nbins", default="30", help="# Histogram bins")
    parse.add_option("--lambda_start",dest="lmb_start", default=0 ,help="lambda starts at this value")
    parse.add_option("--lambda_end",dest="lmb_end",default=1,help="lambda ends at this value")
    parse.add_option("--first",dest="first_n", default=1, help="first trajectory")
    parse.add_option("--last",dest="last_n",default=200, help="last trajectory")
    parse.add_option("-r",dest="rand", help="take # random trajectories")
    parse.add_option("-T",dest="temp",default=300,help="Temperature [K]")
    parse.add_option("--rev",dest="rev", action="store_true",help="reverse state B")
    parse.add_option("--integ_only",dest="integ_only",action="store_true",help="Only integrate. No analysis")

    #files
    parse.add_option('-o',dest="res_out", default='results.dat', help="write results output file")
    parse.add_option('--o0',dest="integ0_o",default='integ0.dat', help="Output integration file 0->1")
    parse.add_option('--o1',dest="integ1_o",default='integ1.dat', help="Output integration file 1->0")
    parse.add_option('--i0',dest="integ0_r",help="Input integration file 0->1")
    parse.add_option('--i1',dest="integ1_r",help="Input integration file 0->1")
    parse.add_option('--name0',dest="xvg0",default='dhdl.part0001.xvg',help="Read xvg0 file")
    parse.add_option('--name1',dest="xvg1",default='dhdl.part0001.xvg',help="Read xvg1 file")
    parse.add_option('--pic',dest="pic",default='plot.png',help="Output file for figure")

    (opt, args) = parse.parse_args()
    print opt, args
	
    lda0 = opt.lmb_start
    lda1 = opt.lmb_end
	
    bPlot = opt.showPlot
    nbins = opt.nbins
    bRandom = opt.rand #(bset??)
    bSubset = opt.first_n or opt.last_n or bRandom #(bset?)
    plname = opt.pic
    nrand = opt.rand
    first = opt.first_n
    last = opt.last_n
    name0 = opt.xvg0
    name1 = opt.xvg1
    bReverse = opt.rev
    bInteg = opt.integ_only

    T = opt.temp
	
    if opt.integ0_r :
        print 'Reading intgrated dgdl (%4.2f->%4.2f) from %s' % (lda0, lda1, opt.integ0_r)
        print 'Reading intgrated dgdl (%4.2f->%4.2f) from %s' % (lda1, lda0, opt.integ1_r)
        res_ab = read_integ_file(opt.integ0_r)
        res_ba = read_integ_file(opt.integ1_r)

        if bSubset:
            print 'Will analyze subset of trajectories'
            if not bRandom:
                print 'Selecting trajectories from %d to %d' % (first,last)
                res_ab = res_ab[first:last]
                res_ba = res_ba[first:last]
            else:
                print 'Selecting %d random trajectories' % nrand
                tmp1 = []
                lst1 = []
                ntraj = 0
                while ntraj < nrand:
                    idx = randint(0,len(res_ab)-1)
                    if idx not in tmp1:
                        lst1.append(res_ab[idx])
                        tmp1.append(idx)
                        ntraj+=1
                ntraj = 0
                lst2 = []
                tmp2 = []
                while ntraj < nrand:
                    idx = randint(0,len(res_ba)-1)
                    if idx not in tmp2:
                        lst2.append(res_ba[idx])
                        tmp2.append(idx)
                        ntraj+=1
                res_ab = lst1
                res_ba = lst2
                
    else:
        res_ab  = []
        res_ba = []
        abdir = opt.Adir
        badir = opt.Bdir

        # we read single dgdl files
        ab = get_list(abdir, name0)
        if ab:
            fp0 = open(opt.integ0_o,'w')
        ba = get_list(badir, name1)
        if ba:
            fp1 = open(opt.integ1_o,'w')

        for dgdl in ab:
            l, flag, dlambda = check_dgdl(dgdl, lda0, lda1)
            if flag:
                fname = dgdl.split('/')[-2]
                sys.stdout.write('\r# Integrating %s (0->1)' % dgdl.split('/')[-2])
                res_ab.append( integrate_dgdl(l, fname, fp0, dl = dlambda, lda0=0) )
                sys.stdout.flush()
            else:
                print '\nSkipping %s (not present or invalid)' % dgdl
        print
        for dgdl in ba:
            l, flag, dlambda = check_dgdl(dgdl, lda0, lda1)
            if flag:
                sys.stdout.write('\r# Integrating %s (1->0)' % dgdl.split('/')[-2])
                fname = dgdl.split('/')[-2]
                r = integrate_dgdl(l, fname, fp1, dl = -dlambda, lda0=1, bReverse = bReverse)
                res_ba.append( r )
                sys.stdout.flush()
            else:
                print '\nSkipping %s (not present or invalid)' % dgdl
        print


    # do analysis
    print
    if bInteg:
        sys.exit(0)
    fp = open(opt.res_out,'w')
    meana = average(res_ab)
    meanb = average(res_ba)
    stda = std(res_ab)
    stdb = std(res_ba)
    mean = (meana+meanb)/2.
    print >>fp, '#==============================='
    print >>fp, '#        mean        std'
    print >>fp, '(0->1):   %g          %g' % (meana, stda)
    print >>fp, '(1->0):   %g          %g' % (meanb, stdb)
    print >>fp, '#==============================='
    print len(res_ab)
    print len(res_ba)
    res, error = histo_plot(res_ab, res_ba, nbins, fp, show_plot=bPlot, plname = plname)
    print >>fp, 'Result_CGI: %g # gaussian intersection' % res
    print >>fp, 'Error_CGI: %g' % error
    print >>fp, 'Result_mean: %g # simple mean' % mean
    print 'Running Nelder-Mead simplex algorithm for BAR'
    bar = BAR(res_ab, res_ba,T)
#    print bar[0]
    print >>fp, 'Result_BAR: %g # Bennet Acceptance Ratio' % bar
    err = BAR_err(bar, res_ab, res_ba,T)
    print >>fp, 'Error_BAR: %g ' % err
    fp.close()


if __name__=='__main__':
    main(sys.argv)
    
