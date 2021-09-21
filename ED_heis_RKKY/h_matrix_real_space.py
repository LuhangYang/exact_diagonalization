from basis import *
import numpy as np
from math import cos
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigs

######################### DEFINE GLOBAL VARIABLES HERE ###############3333333
dt = np.dtype(np.complex128)

hdiag = np.zeros(4)  # Diagonal energies
hflip = np.zeros(4)  # off-diagonal terms

hdiag[0] = +0.25
hdiag[1] = -0.25
hdiag[2] = -0.25
hdiag[3] = +0.25

hflip[0] = 0.
hflip[1] = 0.5
hflip[2] = 0.5
hflip[3] = 0.

#################  commonly used functions ############################################
def IBITS(n, i):
    return ((n >> i) & 1)

def find_rep(state,basis):
    idx = -1
    for i in range(len(basis)):
        if(state == basis[i]):
            idx= i
        else:
            continue
    return(idx)


################### r-space H matrix ####################################################
def heis_rkky_r_space(basis,L,lamda,alpha):
    dim = len(basis)
    H = np.zeros(shape=(dim, dim))
    for i in range(dim):
        state = basis[i]
        for delta in range(1, int(L / 2)):
            # Diagonal term
            for site_i in range(L):
                site_j = (site_i - delta + L) % L
                
                if (delta == 1):
                    two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
                    value = hdiag[two_sites]
                    H[i, i] += value
                else:
                    two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
                    value = hdiag[two_sites]
                    
                    H[i, i] += value * (-lamda * ((-1) ** delta) / (delta ** alpha))
        for site_i in range(int(L / 2)):
            site_j = (site_i - int(L / 2) + L) % L
            two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
            value = hdiag[two_sites]
            H[i, i] += value * (-lamda * ((-1) ** (L / 2)) / ((L / 2) ** alpha))

    for i in range(dim):#[1,2,4,5]:
        state = basis[i]
        for delta in range(1, int(L / 2)):
            # Diagonal term
            for site_i in range(L):
                site_j = (site_i - delta + L) % L
                # Off-diagonal term
                if (delta == 1):
                    mask = (1 << site_i) | (1 << site_j)
                    two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
                    value = hflip[two_sites]
                    if (value != 0.):
                        new_state = (state ^ mask)
                        j = find_rep(new_state, basis)
                        H[i, j] += value
                else:
                    mask = (1 << site_i) | (1 << site_j)
                    two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
                    value = hflip[two_sites]
                    if (value != 0.):
                        new_state = (state ^ mask)
                        j = find_rep(new_state, basis)
                        # delta = ((site_i - site_j + L) % L)
                        # delta = min(((site_i-site_j+L)%L),((site_j-site_i+L)%L))
                        H[i, j] += value * (-lamda * ((-1) ** delta) / (delta ** alpha))
        for site_i in range(int(L / 2)):
            site_j = (site_i - int(L / 2) + L) % L
            mask = (1 << site_i) | (1 << site_j)
            two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
            value = hflip[two_sites]
            if (value != 0.):
                new_state = (state ^ mask)
                j = find_rep(new_state, basis)
                H[i, j] += value * (-lamda * ((-1) ** (L / 2)) / ((L / 2) ** alpha))
    
    #print(H)
    d, v = np.linalg.eigh(H)  # Diagonalize the matrix           
    return d,v


def heis_r_space(basis):
    dim = len(basis)
    H = np.zeros(shape=(dim, dim))
    for i in range(dim):
        state = basis[i]
        for delta in range(1, 2):
            # Diagonal term
            for site_i in range(L):
                site_j = (site_i - delta + L) % L
                
                if (delta == 1):
                    two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
                    value = hdiag[two_sites]
                    H[i, i] += value
                else:
                    two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
                    value = hdiag[two_sites]
                    H[i, i] += value * (-lamda * ((-1) ** delta) / (delta ** alpha))

    for i in range(dim):#[1,2,4,5]:
        state = basis[i]
        for delta in range(1, 2):
            # Diagonal term
            for site_i in range(L):
                site_j = (site_i - delta + L) % L
                # Off-diagonal term
                if (delta == 1):
                    mask = (1 << site_i) | (1 << site_j)
                    two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
                    value = hflip[two_sites]
                    if (value != 0.):
                        new_state = (state ^ mask)
                        j = find_rep(new_state, basis)
                        H[i, j] += value
                else:
                    mask = (1 << site_i) | (1 << site_j)
                    two_sites = IBITS(state, site_i) | (IBITS(state, site_j) << 1)
                    value = hflip[two_sites]
                    if (value != 0.):
                        new_state = (state ^ mask)
                        j = find_rep(new_state, basis)
                        H[i, j] += value * (-lamda * ((-1) ** delta) / (delta ** alpha))
      
    #print(H)
    d, v = np.linalg.eigh(H)  # Diagonalize the matrix           
    return d



if __name__ == "__main__":
    L=8
    Nup=4
    basisq = build_all_basis(L, Nup)
    valu,vect = heis_rkky_r_space(basisq,L,1.,2.)
    print(valu)
    







