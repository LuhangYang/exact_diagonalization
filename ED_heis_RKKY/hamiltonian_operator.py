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

#################commonly used
def IthBITinN(n,i):
    return ((n >> i) & 1)

def spintrans(state,L,bit): #translate the left #n bit to the right
    b = state >> L-bit
    lastbit=0
    for i in range(bit-1):
        lastbit += 2**i
    state = state << bit
    state = state & (2**L-1)
    state = state | b
    return state


def ch(vector,k):#kth subblock ranges from 1 to L, corresponding to [0,2pi)
    thetak = 2 * np.pi * (k - 1) / L
    c = np.exp(-1j * (vector) * thetak)
    return(c)

def find_rep(state,basis):
    idx = -1
    for i in range(len(basis)):
        if(state == basis[i]):
            idx= i
        else:
            continue
    return(idx)

############################ H operator ###########################################

def hami_heis(state_in, term,L): # H Heis operator act on state_in,
    x = 0.
    if(term == 0):
        # diagonal term
        state_out = state_in
        for site_i in range(L):
            site_j = (site_i - 1 + L) % L  # (sitei-1+L)mod L,site_j is the right side bit on site_i
            ijnotsame = ((IthBITinN(state_in, site_i) ^ IthBITinN(state_in, site_j)))
            x += ijnotsame * (-0.25) + (1 - ijnotsame) * (0.25)
    else:
        # off-diagonal term
        # do the spin flip at position x
        state_out = []
        x = []
        for site_i in range(L):
            site_j = (site_i - 1 + L) % L  # (sitei-1+L)mod L,site_j is the right side bit on site_i
            mask = (1 << site_i) | (1 << site_j)  # 1 on site i and 1 on site j
            state_outi = (state_in ^ mask)  # flip the two adjacent
            xi = (IthBITinN(state_in, site_i) ^ IthBITinN(state_in, site_j)) * 0.5
            state_out.append(state_outi)
            x.append(xi)
    return (state_out,x) 


def hami_heis_rkky(state_in, term,lamda,alpha,L): # H operator act on state_in,term=0:diag,term=1:spin flip
    if(term == 0):
        # diagonal term
        x = 0.
        state_out = state_in
        for delta in range(1,int(L/2)):
            for site_i in range(L):
                site_j = (site_i - delta + L) % L
                ijnotsame = ((IthBITinN(state_in, site_i) ^ IthBITinN(state_in, site_j)))
                if(delta==1):  # (sitei-1+L)mod L,site_j is the right side bit on site_i
                    x += ijnotsame * (-0.25) + (1 - ijnotsame) * (0.25)
                else:
                    x += (ijnotsame * (-0.25) + (1 - ijnotsame) * (0.25)) * (-lamda*((-1)**delta)/(delta**alpha))
        for site_i in range(int(L/2)):
            site_j = (site_i - int(L/2) + L) % L
            ijnotsame = ((IthBITinN(state_in, site_i) ^ IthBITinN(state_in, site_j)))
            x += (ijnotsame * (-0.25) + (1 - ijnotsame) * (0.25)) * (-lamda * ((-1) ** (L/2)) / ((L/2) ** alpha))
    elif(term==1):
        # off-diagonal term
        # do the spin flip at position x
        state_out = []
        x = []
        for delta in range(1,int(L/2)):
            for site_i in range(L):
                site_j = (site_i - delta + L) % L
                mask = (1 << site_i) | (1 << site_j)  # 1 on site i and 1 on site j
                state_outi = (state_in ^ mask)  # flip the two adjacent
                if(delta==1):
                    xi = (IthBITinN(state_in, site_i) ^ IthBITinN(state_in, site_j)) * 0.5
                    state_out.append(state_outi)
                    x.append(xi)
                else:  # (sitei-1+L)mod L,site_j is the right side bit on site_i
                    xi = (IthBITinN(state_in, site_i) ^ IthBITinN(state_in, site_j)) * 0.5 * (-lamda*((-1)**delta)/(delta**alpha))
                    state_out.append(state_outi)
                    x.append(xi)
        for site_i in range(int(L/2)):
            site_j = (site_i - int(L/2) + L) % L
            mask = (1 << site_i) | (1 << site_j)  # 1 on site i and 1 on site j
            state_outi = (state_in ^ mask)  # flip the two adjacent
            xi = (IthBITinN(state_in, site_i) ^ IthBITinN(state_in, site_j)) * 0.5 * (
                        -lamda * ((-1) ** (L/2)) / ((L/2) ** alpha))
            state_out.append(state_outi)
            x.append(xi)
    return (state_out,x) # after the operator's action, becomes stateout with coefficient x



if __name__ == "__main__":
    L=8
    Nup=4
    







