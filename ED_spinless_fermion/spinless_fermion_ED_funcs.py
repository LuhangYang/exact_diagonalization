import numpy as np
from math import cos
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigs

######################### DEFINE GLOBAL VARIABLES HERE ##################
dt = np.dtype(np.complex128)
L=16
Nocu = 4#int(L/2)


hflip = np.zeros(4)  # off-diagonal terms

hflip[0] = 0.
hflip[1] = -1.0
hflip[2] = -1.0
hflip[3] = 0.

#################commonly used functions #################################
def IBITS(n, i): #from right to left
    return ((n >> i) & 1)


def find_rep(state,basis):
    idx = -1
    for i in range(len(basis)):
        if(state == basis[i]):
            idx= i
        else:
            continue
    return(idx)


################# prepare basis #######################################3

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
    return basis




def spinless_r_space_V4(basis,V4):
    dim = len(basis)
    H = np.zeros(shape=(dim, dim))
    for i in range(dim):
        state = basis[i]
        for site_i in range(L):
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
    print("H: ")
    print(H)
    d, v = np.linalg.eigh(H)  # Diagonalize the matrix           
    return d,v


def spinless_r_space_V4_rise_k(basis,V4,mu):
    dim = len(basis)
    H = np.zeros(shape=(dim, dim),dtype=dt)
    for i in range(dim):
        state = basis[i]
        for site_i in range(L):
            site_j = (site_i + 1 + L) % L
            # Diagonal term
            site_m = (site_i + 2 + L) % L 
            site_n = (site_i + 3 + L) % L 
            if( IBITS(state, site_i)==1 and IBITS(state, site_j)==1 and IBITS(state, site_m)==1  and IBITS(state, site_n)==1):
                H[i, i] += V4
            # Diagonal terms from mu |k><k|
            for kout in range(round(L/2+1),L):
                ko = -np.pi+ kout*2.*np.pi/L
                if(IBITS(state, site_i)==1): 
                    H[i, i] += mu
                    #print(site_i,site_i,i,i, mu ) 
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
                #print("nearest: ",site_i,site_j,i,j,value*((-1)**fermion_num) ) 

            # Off-diagonal terms from mu |k><k|   
            for js in range(L):
                if(js==site_i): continue
                for kout in range(round(L/2+1),L):
                    ko = -np.pi+ kout*2.*np.pi/L 

                    fermion_num = 0
                    for num_i in range(site_i+1,L):
                        if(IBITS(state,num_i)==1):fermion_num+=1
                    for num_i in range(js+1,L):
                        if(IBITS(state,num_i)==1):fermion_num+=1
                    
                    if(IBITS(state,site_i)==1 and IBITS(state,js)==0):#i<j CdjCi
                        if(site_i>js): fermion_num+=1
                        mask = (1 << site_i) | (1 << js)
                        two_sites = IBITS(state, site_i) | (IBITS(state, js) << 1)
                        value = hflip[two_sites]
                        if (value != 0.):
                            new_state = (state ^ mask)
                            j = find_rep(new_state, basis)
                            H[i, j] += value*-1.*np.exp(-1j*ko*(js-site_i))* mu * ((-1)**fermion_num) 
                            
                    
    print("H: ")
    print([H[ii][ii] for ii in range(len(H))])
    d, v = np.linalg.eigh(H)  # Diagonalize the matrix           
    return d,v


if __name__ == "__main__":
    build_basis(63, 5)
