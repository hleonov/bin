import sys, os
from glob import glob
from numpy import *
from pylab import *
import random as rnd
def simple_xy_plot(x,y,xlab, ylab, name, result, err, ylimit = False):
    figure(figsize=(8,8))
    plot(x,y,'rd-', lw=2)
    xlabel(xlab, fontsize=20)
    ylabel(ylab, fontsize=20)
    title(r"$\int_0^1\langle\frac{dG}{d\lambda}\rangle_{\lambda}d\lambda$ = %.2f +/- %.2f kJ/mol" % (result, err))
    grid(lw=2)
    xx = gca()
    if ylimit:
        ylim( ylimit)
    for x in xx.spines.values():
        x.set_lw(2)
    savefig(name, dpi=600)


def plot_random_blocks( run_dic, block_file,avx, avy, result, err, ylimit = False):
    figure(figsize=(8,8))
    outf = 'random_blocks.png'
    xlab = r'$\lambda$'
    ylab = r"dG/d$\lambda$"
    block_res = []
    for lda, p in run_dic.items():
        blck = read_blocks(p, block_file)
        block_res.append( (lda, blck ))
    block_res.sort(lambda a, b: cmp(a[0],b[0]) )
    xx = []
    for i in range(1000):
        dG, x, y = random_value_from_block( block_res, with_data = True )
        plot(x,y,'k-', lw=.5, alpha = .01)        
    plot(avx,avy,'rd-', lw=2)
    if ylimit:
        ylim( ylimit)
    xlabel(xlab, fontsize=20)
    ylabel(ylab, fontsize=20)
    title(r"$\int_0^1\langle\frac{dG}{d\lambda}\rangle_{\lambda}d\lambda$ = %.2f +/- %.2f kJ/mol" % (result, err))

#    title("Curves from random blocks")
    grid(lw=2)
    xx = gca()
    for x in xx.spines.values():
        x.set_lw(2)        
    savefig(outf, dpi=600)
    

def check( lst ):
    
    files = ['fe_convergence.txt',
             'fe_avg.txt',
             'fe_block_1000.txt']#,
             #'fe_block_2000.txt']
    for d in lst:
        for f in files:
            p = os.path.join(d, f)
            if not os.path.isfile( p ):
                print 'file %s missing in directory %s -> Exiting' % (d,f)
                sys.exit(1)

#reads fe_avg, only the avg itself, not number of points        
def read_result(p):
    f = os.path.join(p,'fe_avg.txt')
    arr = open(f).read().split()
    r,s = (float(arr[0]), float(arr[1]))
    return r,s

#read fe_block file, stores only second column in a list
def read_blocks(p, block_file):
    f = os.path.join( p, block_file)
    l = open(f).readlines()
    r = []
    for line in l:
        r.append( float(line.split()[1] ))
    return r

def do_dgdl( results ):
    x = map(lambda a: a[0], results)	#list of x's - dimension to integrate along
    y = map(lambda a: a[1], results)	#list of y's - values to integrate
    dG =  trapz(y,x)
    return dG, x, y

#generates a random list of (lambda, avg_val) from block_avg, and integrates over it.
def random_value_from_block( block_res, with_data = False ):
    #print "in random value"
    data = []
    n_blocks = len(block_res)				#number of lambdas
	#accumulate one  average value per lambda at random [(0.0, rand_fe0), (0.1, rand_fe1)..]
    for i in range(n_blocks):				#range 1..number of lambdas
        size = len(block_res[i][1]) 	#block_res[i][1]=list of block avg     
        ri = rnd.randint(0,size-1)			#choose one of these points
        data.append( (block_res[i][0], block_res[i][1][ri]) )	#lmb_i, + tuple i -> list[ri]
    dg, x, y = do_dgdl( data )
    if with_data:
        return dg, x, y
    else:
        return dg

#generate 1000 values of dG from random block averages (corresponding to their lambda)
#Their std =error. The error will vary more if the block averages differ much over time
def error_from_block_aver(run_dic, block_file):
    block_res = []
    for lda, p in run_dic.items():
        blck = read_blocks(p, block_file)	#list of averages
        block_res.append( (lda, blck ))	#list of tuples where 2nd element=list
    block_res.sort(lambda a, b: cmp(a[0],b[0]) )    
    xx = []
    for i in range(100000):
        dG = random_value_from_block( block_res )
        #print  "random dG="+str(dG)
        xx.append( dG )  
    print str(average(xx)), str(std(xx))  
    return std(xx)
    

#reads fe_avg for all lambdas, stores and sorts according to lambda value (first)
#the lambda keyword in the sort defines an anonymous sort function
def read_dGdl( run_dic) :
    results_all = []
    std_all = []
    for lda, p in run_dic.items():
        r,s = read_result(p)
        results_all.append( (lda, r ))
        std_all.append((lda, s))
    results_all.sort(lambda a, b: cmp(a[0],b[0]) )
    std_all.sort(lambda a, b: cmp(a[0],b[0]) )
    return results_all, std_all

def read_convergence_file( f ):
    l = open(f).readlines()
    r = []
    for line in l:
        entr = line.split()
        r.append( (float(entr[0]), float(entr[1]) ) )
    return r
                
def dG_over_time( run_dic, infile, outfile ):
    fp = open(outfile,'w')
    dG_vs_time = []
    for lda, p in run_dic.items():
        r = read_convergence_file( os.path.join(p, infile) )
        dG_vs_time.append( (lda, r) )
    dG_vs_time.sort(lambda a, b: cmp(a[0],b[0]) )
	 
	 #dG_vs_time = [(0.0, [(1000.0, 12), (2000.0, 13), (3000.0, 13)]), 
	 #					 (0.05 [(1000.0, 94), (2000.0, 94), (3000.0, 93)])...
    # the values here are the average over the cumulative dgdl up until the specified time-point.
    # for each time point, we do an integration over lambda and get the dG convergence over time.
    # The average in each time point is not independent from the previous, and so is the dG. (averaged backwards -done before in DTI_analysis.py)
    min_size = 100
    for lda, lst in dG_vs_time:
        if len(lst) < min_size:
            min_size = len(lst)
    for i in range(min_size):
        time = dG_vs_time[0][1][i][0]		
        lda_vals = map(lambda a: a[0], dG_vs_time)
        dgdl_vals = map(lambda a: a[1][i][1], dG_vs_time)
        dG = trapz( dgdl_vals, lda_vals)
        print >>fp, time, dG 

def get_max_number_of_blocks(run_dic, block_file):
    block_res = []
    for lda, p in run_dic.items():
        blck = read_blocks(p, block_file)
        block_res.append( (lda, blck ))
    block_res.sort(lambda a, b: cmp(a[0],b[0]) )
    min_size = min( map(lambda a: len(a[1]), block_res ) )
    return min_size

def dG_from_last_blocks( run_dic, block_file, nblocks ):
    block_res = []
	 #block_res = [(0.0, [1290, 1280, 1276, ..]), 
    #             (0.05,[1240, 1282, 1278, ..]), 
    # The values here are the 1000 (or 2000) block averages computed for separate blocks from the data
    # For each lambda, we average the dgdl of the last X blocks (depends on the nblocks input). 
    # Averaged forward.-- done here.
    for lda, p in run_dic.items():
        blck = read_blocks(p, block_file)
        block_res.append( (lda, blck ))
    block_res.sort(lambda a, b: cmp(a[0],b[0]) )
    min_size = min( map(lambda a: len(a[1]), block_res ) )
    first_block = min_size - nblocks
#    first_block = 15
    res = []
    for b in block_res:
        lda = b[0]
        dGdl = average( b[1][first_block:] )
        res.append( (lda, dGdl) )
    dG = trapz( map(lambda a:a[1], res), map(lambda a:a[0], res) )
    return dG, nblocks
                        
def calc_big_std(std_all):
   lda_vals = map(lambda a: a[0], std_all)
   err_vals = map(lambda a: a[1], std_all)
   #print err_vals
   err_sqr = (array(err_vals))**2
   sig = sqrt(trapz(err_sqr, lda_vals))
   return sig

#Sig_lmb = sqrt[(1/T-1)*SUM_t[Bavg_t-Avg_lmb]^2]/sqrt(T)
def calc_std_from_block_avg(run_dic, block_file, all_avg):
   var_vals = []
   lda_vals = map(lambda a: a[0], all_avg)
   avg_vals =  map(lambda a: a[1], all_avg)
   #for lda, p in run_dic.items():
   for i in range(len(lda_vals)):
      blck = read_blocks(run_dic[lda_vals[i]], block_file)
      s = 0
      for B_avg in blck: 
         s=s+(B_avg-avg_vals[i])**2    
      sig=sqrt(s*(1/(float)(len(blck))))/(float)(sqrt(len(blck)+1))
      #print lda_vals[i], sig
      var_vals.append(sig**2)
   #print var_vals
   sig = sqrt(trapz(var_vals, lda_vals))
   return sig
        
#============================== MAIN ================================        
base_name = sys.argv[1]
#lst = glob('lambda_*.*')
lst = glob("%s*.*" % base_name)
check(lst)

run_dic = {}

for d in lst:
    lda = float(d.split('_')[1])
    run_dic[lda] = d

#read fe_avg from all lambdas
results_all, std_all = read_dGdl( run_dic )

#numerical integration over averages from the entire trajectory    
dgdl_all, x, y = do_dgdl(results_all)

#Error estimation
sig1 = calc_big_std(std_all)
sig2 = calc_std_from_block_avg(run_dic, 'fe_block_1000.txt', results_all)
err = error_from_block_aver( run_dic, 'fe_block_1000.txt')
err2 = error_from_block_aver( run_dic, 'fe_block_2000.txt')

print "delta G = "+str(dgdl_all)
#print "Sigma from entire trj: "+str(sig1)
print "Sigma from block avg: "+str(sig2)
outf = 'lda_vs_dGdl.txt'
fp = open(outf, 'w')
for i in range(len(x)):
    print >>fp, x[i], y[i]

outf = 'lda_vs_dGdl.png'
try:
    yl = ( float(sys.argv[2]), float(sys.argv[3]) )
except:
    yl = False
           
simple_xy_plot(x,y,r'$\lambda$', r"dG/d$\lambda$",outf, dgdl_all, err, yl)
if err != 0:
    plot_random_blocks(run_dic, 'fe_block_1000.txt', x, y, dgdl_all, err, yl)

fp = open('fe_result.txt','w')
print >>fp, 'Result dG = ', round(dgdl_all,2), 'kJ/mol', 'Err(1000) = ', round(err,2), 'kJ/mol', 'Err(2000) = ', round(err2,2), 'kJ/mol'
print >>fp, "Sigma from block avg: "+str(sig2)
dG_over_time( run_dic, 'fe_convergence.txt', 'dG_vs_time_cum.txt' )
dG_over_time( run_dic, 'fe_block_1000.txt', 'dG_vs_time.txt' )

fp2 = open("dG_from_last_blocks.txt",'w')
max_block = get_max_number_of_blocks( run_dic, 'fe_block_1000.txt')
for i in range(1, max_block+1):
    dG, ndata = dG_from_last_blocks( run_dic, 'fe_block_1000.txt', i)
    print >>fp, 'dG (block) = ', round(dG,2), 'kJ/mol', 'from last ', ndata, 'blocks'
    print >>fp2, (max_block-ndata+1)*1000 , round(dG,2)
