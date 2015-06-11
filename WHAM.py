# Calculate PMF using the WHAM equations.
#
# conventions follow Cumar et al. In particular the choice of R for n_sims is particularly miserable...
# Note that the gas constant is also denoted by R (global)
#
# INPUT:
# ======
# E[i][j] = total energy of step j in simulation i.
# zeta[i][j] = reaction coordinate at step j, simulation i.
# V[l] = restraining potential *function* l; V[l] = V[l](zeta)
# lam[i][l] = weight of potential l at simulation i.
# T[i] = temperature at simulation i.
# bins = binning of the reaction coordinate; bins [0,1,2] correspond to zeta values 0,1.
#
# CONVENTIONS:
# R = number of simulations 
# n[i] = number of snapshots in simulation i.
# L = number of *restraining* potentials 
#
# CALCULATED VALUES:
# ==================
# E[i][j][l] = total energy of step j in simulation i under Hamiltonian l. Only E[i][j] is an array.
# f[i] = free energy at simulation i. #exp(-A[i]) where A[i] is the free energy of system at simulation i.
# N[i][b] = value of the histogram at simulation i and bin b

from numpy import *

R = .0083144 # kJ/mol*K

class WHAM:
    def __init__(self, zeta, V, lam, T = "Room", f = "zeros"):
        self.zeta = zeta
	self.V = V
	self.lam = lam
	self.R = len(zeta)
	if T =="Room":
	    T = [300.] * self.R
	self.beta = array([1./(R * t) for t in T])
	self.n = map(len, zeta)
	self.L = len(lam[0])
	if f=="zeros":
	    self.f = [0.] * self.R
	else:
	    self.f = f

    def set_energies(self, E):
        self.E = E
	
    def get_energies(self, E):
        self.E = []
	for i in range(self.R):
            self.E.append([])
	    for j in range(self.n[i]):
	        V = array([self.V[l](self.zeta[i][j]) for l in range(self.L)])
		self.E[i].append(array([sum((self.lam[l] - self.lam[i])*V) + E[i][j] for l in range(self.R)]))
 	return
    
    # Perform one iteration of free energy calculation.
    # Follows notation of Cumar et al.
    def getFiter(self):
        newf = [0.] * self.R
	for i in range(self.R):
	    for k in range(self.R):
	        for t in range(self.n[k]):
		    lognum = -self.beta[i]*self.E[k][t][i]
		    den = sum([self.n[m] * exp(self.f[m]-self.E[k][t][m]-lognum) for m in range(self.R)])
		    newf[i] += 1./den
	self.f = array([exp(-x) for x in newf])  
	
    # bin Get WHAM results for zeta binned according to bins, with lam values test_lam	
    def bin(self, bins):
        self.bins = bins
	self.N = []
	for i in range(self.R):
	    self.N.append(array([len([z for z in self.zeta[i] if bins[b] <= z < bins[b+1]]) for b in range(len(bins) - 1)]))
	
    # Get Probability density function for lam values test_lam (a natural default is [0.]*L) and temperature T	    
    def WHAM(self, test_lam, T = 300):
        beta = 1./ (R * T)
	self.P = {} 
	for b in range(len(self.bins)-1):
	    z = self.bins[b]
	    num = sum([self.N[k][b] for k in range(self.R)]) * \
	          exp(-beta*sum([test_lam[j] * self.V[j](z) for j in range(self.L)]))
	    den = sum([self.n[m] * exp(self.f[m]-beta*sum([self.lam[m][j]*self.V[j](z) for j in range(self.L)])) for m in range(self.R)])
            self.P[z] = num / den
	
	return

    def checkConvergence(self):
        self.conv = []
        for i in range(self.R):
	    self.WHAM(self.lam[i], 1./(R*self.beta[i]))
	    self.conv.append(exp(-self.f[i]) - sum(self.P.values()))

    def WHAMiter(self):
        newF = array([0.]*self.R)
	for i in range(self.R):
	    self.WHAM(self.lam[i], 1./(R*self.beta[i]))
	    newF[i] = -log(sum(self.P.values()))
	self.f = newF
	
	
        
