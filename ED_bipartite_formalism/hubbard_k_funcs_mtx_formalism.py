import numpy as np
import matplotlib.pyplot as plt

Lx = 4
Ly = 1
L = Lx*Ly
#N = 8
#S = 6
dt = np.dtype(np.complex128)

t = -1.
U = 4.


class BoundaryCondition:
      OBC, PBC = range(2)

xBC = BoundaryCondition.OBC
#yBC = BoundaryCondition.PBC


hflip = np.zeros(4)  # off-diagonal terms

hflip[0] = 0.
hflip[1] = t
hflip[2] = t
hflip[3] = 0.

def find_rep(state,basis):
    idx = -1
    for i in range(len(basis)):
        if(state == basis[i]):
            idx= i
        else:
            continue
    return(idx)

def spinless_r_space(basis):
    dim = len(basis)
    H = np.zeros(shape=(dim, dim))
    for i in range(dim):
        state = basis[i]
        print(bin(state))
        for site_i in range(L-1):
            site_j = (site_i + 1 + L) % L
       
            # Off-diagonal term
            mask = (1 << site_i) | (1 << site_j)
            two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
            #print(bin(state),site_i,site_j,two_sites)
            fermion_num = 0
            for num_i in range(site_i+1,L):
                if(IBITS(state,num_i)==1):fermion_num+=1
            for num_i in range(site_j+1,L):
                if(IBITS(state,num_i)==1):fermion_num+=1
            if(site_i>site_j and IBITS(state,site_i)==1 and IBITS(state,site_j)==0):fermion_num+=1#i<j CdjCi
            if(site_i<site_j and IBITS(state,site_i)==0 and IBITS(state,site_j)==1):fermion_num+=1 #i<j CdjCi
            value = hflip[two_sites]
            if (value != 0.):
                new_state = (state ^ mask)
                j = find_rep(new_state, basis)
                H[i, j] += value * ((-1)**fermion_num)       
    print("Aup: ")
    print(H)
    #d, v = np.linalg.eigh(H)  # Diagonalize the matrix           
    return H#d,v





class BipartiteSystem(object): 

    def __init__(self, _nsites,basisUp,dt,U):
        # the bipartite here is based on SPIN, not position!
        self.dt = dt
        self.nsites = _nsites
        self.dim_basis = len(basisUp) #assume Nup=Ndn
        self.u0 = np.zeros(shape=(self.dim_basis,self.dim_basis),dtype=self.dt) # 
        self.Aup = spinless_r_space(basisUp)
        #print("Spinless results: ",np.linalg.eigh(self.Aup))
        self.Adn = np.copy(self.Aup)
        for i in range(self.dim_basis):
            for j in range(self.dim_basis):
                for n in range(self.nsites):
                    self.u0[i,j] += U * IBITS(basisUp[i],n) * IBITS(basisUp[j],n)
        #print("D: ", self.u0)
        self.psi = np.zeros(shape=(self.dim_basis,self.dim_basis),dtype=self.dt) # g.s. wave function
    def product(self,psi): #dot product bwtween H and psi vector
        npsi = np.zeros(shape=(self.dim_basis,self.dim_basis))
        for i in range(len(psi)):
            for j in range(len(psi[0])):
                npsi[i][j] = self.u0[i][j] * psi[i][j]
                for k in range(len(psi)):
                    npsi[i][j] += self.Aup[i][k] * psi[k][j]
                    npsi[i][j] += self.Adn[j][k] * psi[i][k]
        return npsi


def psi_dot_psi(psi1, psi2): #<psi1|psi2>
    x = 0.
    for i in range(psi1.shape[0]):
        for j in range(psi2.shape[1]):
            x += np.conjugate(psi1[i,j])*psi2[i,j]
    print(x)
    return x

def lanczos_customize(m, seed, maxiter, tol, use_seed = False, force_maxiter = False):
    x1 = seed
    x2 = seed
    gs = seed
    a = np.zeros(100,dtype=m.dt)
    b = np.zeros(100,dtype=m.dt)
    z = np.zeros((100,100),dtype=m.dt)
    lvectors = []
    control_max = maxiter
    e0 = 9999

    if(maxiter == -1):
        force_maxiter = False

    if(control_max == 0):
        gs = 1
        maxiter = 1
        return(e0,gs)
    
    x1[:,:] = 0
    x2[:,:] = 0
    gs[:,:] = 0
    a[:] = 0.0
    b[:] = 0.0
    if(use_seed):
        x1 = seed
    else:
        for i in range(x1.shape[0]):
            for j in range(x1.shape[1]):
                x1[i,j] = (2*np.random.random()-1.)

#    x1[:,:] = 1
    b[0] = psi_dot_psi(x1,x1)
    b[0] = np.sqrt(b[0])
    x1 = x1 / b[0]
    x2[:] = 0
    b[0] = 1.

    e0 = 9999
    nmax = min(99, maxiter)
    print("nmax: ", nmax)
    for iter in range(1,nmax+1):
        eini = e0
        if(b[iter - 1] != 0.):
            aux = x1
            x1 = -b[iter-1] * x2
            x2 = aux / b[iter-1]
        #print("state: ",x2)
        aux = m.product(x2)

        x1 = x1 + aux
        a[iter] = psi_dot_psi(x1,x2)
        x1 = x1 - x2*a[iter]

        b[iter] = psi_dot_psi(x1,x1)
        b[iter] = np.sqrt(b[iter])
        lvectors.append(x2)                                                  
#        print ("Iter =",iter,a[iter],b[iter])
        z.resize((iter,iter))
        z[:,:] = 0
        for i in range(0,iter-1):
            z[i,i+1] = b[i+1]
            z[i+1,i] = b[i+1]
            z[i,i] = a[i+1]
        z[iter-1,iter-1]=a[iter]
        d, v = np.linalg.eig(z)

        col = 0
        n = 0
        e0 = 9999
        for e in d:
            if(e < e0):
                e0 = e
                col = n
            n+=1
        e0 = d[col]
        
       
        print ("Iter = ",iter," Ener = ",e0) #check the convergence of the gs
        if((force_maxiter and iter >= control_max) or (iter >= gs.shape[0]*gs.shape[1] or iter == 99 or abs(b[iter]) < tol) or ((not force_maxiter) and abs(eini-e0) <= tol)):
            # converged
            gs[:,:] = 0.
            for n in range(0,iter):
                gs += v[n,col]*lvectors[n]

            print ("E0 = ", e0)
            maxiter = iter
            return(e0,gs) # We return with ground states energy

    return(e0,gs)




def IBITS(n,i): #count from right site from 0.
    return ((n >> i) & 1)


def build_all_basis(param_L, Nup): # Using quantum number m to build bases
    L = param_L
    maxdim = 2 ** L
    basis = []
    dim = 0
    for state in range(int(maxdim)):
        n_ones = 0
        for bit in range(L):
            if ((state & (1 << bit)) != 0):  # the bit(th) bit of state is 1
                n_ones += 1
        if (n_ones != Nup):
            continue
        else:
            basis.append(state)
            dim += 1
    print("dim: ",dim)
    return dim,basis




