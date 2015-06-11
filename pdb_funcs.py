from Numeric import *
from math import *
from string import split
import os

execfile("/home/daniel/proj/hook/src/funcs.py")

class Atom:
	def __init__(self, name, coor):
		self.name = name[:]
		while self.name[0]==' ':
			self.name = self.name[1:]
		while self.name[-1]==' ':
			self.name = self.name[:-1]
		self.c = coor
		self.x, self.y, self.z = self.c[0], self.c[1], self.c[2]
	
	def update(self):
		self.x, self.y, self.z = self.c[0], self.c[1], self.c[2]
		return
		
	def __str__(self):
		return fstr(self.x, 7) + ' ' + fstr(self.y, 7) + ' ' + fstr(self.z, 8) + '\n'

	def __repr__(self):
		return fstr(self.x, 7) + ' ' + fstr(self.y, 7) + ' ' + fstr(self.z, 8) + '\n'		 

	def __getitem__(self, i):
		return self.c[i]
	
	def remake(self):
		self.c = array(self.c)
		self.x, self.y, self.z = self.c[0], self.c[1], self.c[2]

	def locate_origin(self, c):
		if type(c) != type(array(range(10))):
			c = array(c)
		self.c -= c
		self.remake()

	def rotx(self, theta):
		self.c = rotx(self.c, theta)
		self.remake()

	def roty(self, theta):
		self.c = roty(self.c, theta)
		self.remake()

	def rotz(self, theta):
		self.c = rotz(self.c, theta)
		self.remake()

			
class AtomInHelix:
	def __init__(self, n, atom, aa, res):
		self.n = n
		self.atom = atom
		self.aa = aa
		self.res = res
	def __str__(self):
		return 'ATOM'+fstr(self.n+1, 7)+fstr(self.atom.name, 5)+' '+fstr(self.aa,5,'r')+fstr(self.res, 4)+'    '+str(self.atom)
	def __repr__(self):
		return 'ATOM'+fstr(self.n+1, 7)+fstr(self.atom.name, 5)+' '+fstr(self.aa,5,'r')+fstr(self.res, 4)+'    '+str(self.atom)

class Helix:
	def __init__(self,hlx=[]):
		self.copy(hlx)
		
	def remake(self):
		self.backbone = [x for x in self.l if x.atom.name in ['CA','C','N']]
		self.c_alpha =  [x for x in self.l if x.atom.name == 'CA']

	def process(self):
            self.seq = self.sequence()
            self.n_res = self.seq
            self.res = [self.get_res(x) for x in self.res_nums()]
		
	def load(self, file):
		lines = open(file).readlines()
		layout = array([0,6,12,16,17,21,22,26,30,38,46,54])
		pdb = [map(lambda x,y,w=w: w[x:y], layout[:-1], layout[1:]) for w in lines]
		pdb = [x for x in pdb if x and (x[0]=='ATOM  ' or x[0]=='HETATM')]
		self.l = [AtomInHelix(int(x[1])-1,Atom(x[2], map(float, x[8:11])),x[4],int(x[6])) for x in pdb]
		self.remake()
		return
	
	def read(self, string):
		lines = split(string, '\n')
		layout = array([0,6,12,16,17,21,22,26,30,38,46,54])
		pdb = [map(lambda x,y,w=w: w[x:y], layout[:-1], layout[1:]) for w in lines]
		pdb = [x for x in pdb if x and (x[0]=='ATOM  ' or x[0]=='HETATM')]
		self.l = [AtomInHelix(int(x[1])-1,Atom(x[2], map(float, x[8:11])),x[4],int(x[6])) for x in pdb]
		self.remake()

	def copy(self, hlx):
		self.l = [AtomInHelix(x.n, Atom(x.atom.name, array(x.atom.c).copy()), x.aa, x.res) for x in hlx]
		self.remake()
	
	def append(self, hlx):
		self.l += [AtomInHelix(x.n, Atom(x.atom.name, array(x.atom.c).copy()), x.aa, x.res) for x in hlx]
		self.remake
		
        def __add__(self, hlx):
            h = Helix(self)
            h.append(hlx)
            return h
        
	def __getitem__(self, i):
		return self.l[i]

        def __getslice__(self,i,j):
            return Helix(self.l[i:j])
        
	def __len__(self):
		return len(self.l)
		
	def remove_hidrogens(self):
		self.l = [x for x in self if x.atom.name[0]!='H']
		return
		
	def renum(self):
		begres = self[0].res
		for i in range(len(self.l)):
			self[i].n = i
			self[i].res -= begres - 1
		self.remake()
	
	def locate_origin(self, c):
		for i in self:
			i.atom.locate_origin(c)
		self.remake()
		
	def rotx(self, theta):
		for i in self:
			i.atom.rotx(theta)
		self.remake()
			
	def roty(self, theta):
		for i in self:
			i.atom.roty(theta)
		self.remake()

	def rotz(self, theta):
		for i in self:
			i.atom.rotz(theta)
		self.remake()
	
	def align_z(self):
		print "Better not to trust align_z!"
		self.rotz(y_ang(proj_xy(self.director())))
		self.rotx(z_ang(self.director()))

	def rasmol(self):
		self.save("tmp_helix.pdb")
		os.system("/usr/local/rasmol/bin/rasmol-32 tmp_helix.pdb")
		os.system("rm tmp_helix.pdb")
	def vmd(self):
		self.save("tmp_helix.pdb")
		os.system("~/Desktop/VMD\ 1.8/VMD\ 1.8.app/Contents/Resources/vmd/vmd_MACOSX tmp_helix.pdb")
		os.system("rm tmp_helix.pdb")
		
	def __str__(self):
		return reduce(lambda x,y: x+y, map(str, self), '')

	def __repr__(self):
		return reduce(lambda x,y: x+y, map(str, self), '')

	def save(self, file):
		open(file, "w").write(str(self))

	def beg(self):
		return sum([x.atom.c for x in self.c_alpha[:4]])/4

	def end(self):
		return sum([x.atom.c for x in self.c_alpha[-4:]])/4

	def director(self):
		return self.end() - self.beg()
		
	def sequence(self):
                if not self.l:
                    return []
		res = [self[0].aa]
		for i in range(1,len(self)):
			if self[i].res != self[i-1].res:
				res.append(self[i].aa)
		return res

	def res_nums(self):
            if not self.l:
                return []
            res = [self[0].res]
            for i in range(1,len(self)):
                if self[i].res != self[i-1].res:
                    res.append(self[i].res)
            return res
        
	def get_res(self, n):
            return Helix([x for x in self if x.res==n])
        
	def phi(self):
		myat = [[y.atom.c for y in self.backbone[x:x+4]] for x in range(len(self.backbone)-4) if self.backbone[x].atom.name=='C']
		return map(dihedral, myat)
		
	def psi(self):
		myat = [[y.atom.c for y in self.backbone[x:x+4]] for x in range(len(self.backbone)-4) if self.backbone[x].atom.name=='N']
		return map(dihedral, myat)
		
