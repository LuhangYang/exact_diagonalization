import numpy as np
from hubbard_k_funcs_mtx_formalism import *



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

def spinless_r_space_V4(basis,V4):
    dim = len(basis)
    H = np.zeros(shape=(dim, dim))
    for i in range(dim):
        state = basis[i]
        print(bin(state))
        for site_i in range(L-1):
            site_j = (site_i + 1 + L) % L
            # Diagonal term
            site_m = (site_i + 2 + L) % L 
            site_n = (site_i + 3 + L) % L 
            if( IBITS(state, site_i)==1 and IBITS(state, site_j)==1 and IBITS(state, site_m)==1  and IBITS(state, site_n)==1):
                H[i, i] += V4
              
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
        self.Aup = spinless_r_space_V4(basisUp,0)
        print("Spinless results: ",np.linalg.eigh(self.Aup))
        self.Adn = np.copy(self.Aup)
        for i in range(self.dim_basis):
            for j in range(self.dim_basis):
                for n in range(self.nsites):
                    self.u0[i,j] += U * IBITS(basisUp[i],n) * IBITS(basisUp[j],n)
        print("D: ", self.u0)
        self.psi = np.zeros(shape=(self.dim_basis,self.dim_basis),dtype=self.dt) # g.s. wave function
    def product(self,psi): #dot product bwtween H and psi vector
        npsi = np.zeros(shape=(self.dim_basis,self.dim_basis))
        for i in range(len(psi)):
            for j in range(len(psi[0])):
                npsi[i][j] = self.u0[i][j] * psi[i][j]
                for k in range(len(psi)):
                    npsi[i][j] += self.Aup[i][k] * psi[k][j]
                    npsi[i][j] += self.Adn[j][k] * psi[i][k]
        #npsi = np.dot(self.Aup, psi)
        #npsi += np.dot(psi, self.Adn.transpose())
        print("check psi: ", psi)
        print("check npsi: ",npsi)
        return npsi


######################## build basis for each k sector

dim, basisup = build_all_basis(Lx*Ly,int(Lx*Ly/2)) 
basisdn = np.copy(basisup)
print(dim)
dt = float

###############3# calculate eigen energy, eigen vectors, and H matrix
S = BipartiteSystem(L,basisup,dt,U)

seed = np.eye(dim)
print(seed)
eigenvaluehami = lanczos_customize(S, seed, 20, 1e-8, use_seed = False, force_maxiter = False) 
print("energys :", eigenvaluehami)



################################ Compute correction vector in r-space ##############################

#gs = eigenmatrix_hami.transpose()[0]*1.0


#-1.9531014996351634 #t=-1, U=4
#-4.472135954999584  #t=-1, U=0



