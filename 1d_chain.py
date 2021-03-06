import numpy as np
import matplotlib.pyplot as plt
# from __future__ import print_function
from ipywidgets import interact, interactive, fixed, interact_manual
def odc(l=50, t0=-2., t1=-2.5, e1=-2.):
    k = int(0.5*l) - 1
    a = np.zeros((l,l))
    a[[tuple(list(range(l))+list(range(-1,l-1))),tuple(list(range(-1,l-1))+list(range(l)))]] = t0
    a[[(k+1,k-1,k,k),(k,k,k+1,k-1)]] = t1 
    a[k,k]=e1
    f = np.array([ [ ir*ic  for ir in [ np.any(ja!=0) for ja in a.T ]  ] for ic in [ np.any(ia!=0) for ia in a] ])
    d = [i for i in range(l) if f[i,i]]
    x = a[f].reshape((int(np.sqrt(np.shape(a[f])[0])),int(np.sqrt(np.shape(a[f])[0]))))
    e,q = np.linalg.eig(x)
    uu  = (q*q).T
    kk  = int(np.ceil(np.shape(x)[0]/2)) - 1 
    if kk != 0:
        if len(d)%2 == 0:
            rho =[2.0*sum(itu) for itu in uu[:kk+1,:].T]
        else :
            rho =[2.0*sum(itu[:-1]) + itu[-1] for itu in uu[:kk+1,:].T]
    else:
        rho = [1]
    den = np.ones(l)
    for ind, dd in enumerate(d):
        den[dd] = rho[ind] 
    plt.plot(den)
    plt.show()

if __name__=='__main__':
    odc()
    
    
    
 '''
 import numpy as np
import matplotlib.pyplot as plt 
#Simulate One Dimansional impurity Atom Model 
class Line_Model:
    def __init__(self,natm=11,nimp=1,orgeng=0.,offset=-1.0,conj_ii=0.5,conj_ij=1.0):
        self.num_atm = natm
        self.imp_atm = nimp
        self.off_set = offset
        self.org_eng = orgeng
        self.cnj_i2j = conj_ij
        self.cnj_i2i = conj_ii
        self.eff_hmt = np.matrix(np.zeros((natm,natm)))
        self.eff_hmt[range(natm),range(-1,natm-1)] = conj_ii
        self.eff_hmt[range(-1,natm-1),range(natm)] = conj_ii
        #Hybrid impurity atom in one dimansional atom chain.
        self.eff_hmt[range(nimp),range(nimp)] = offset
        self.eff_hmt[range(nimp+1),range(-1,nimp)] = conj_ij
        self.eff_hmt[range(-1,nimp),range(nimp+1)] = conj_ij
        print(self.eff_hmt)
        self.eig_val, self.eig_vec = np.linalg.eig(self.eff_hmt)
        print(self.eig_val)
        self.occ = np.array([ for l]).T
        self.occ_val = np.array([np.sum(iline) for iline in self.occ])
        print(self.occ_val)
        
test = Line_Model()
'''
