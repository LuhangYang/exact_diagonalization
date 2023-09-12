import numpy as np
from math import cos
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigs


################# some commonly used functions #####################################
def IBITS(n, i):
    return ((n >> i) & 1)

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

def build_all_trans(state,L):
    all_state = []
    for l in range(L):
        all_state.append(spintrans(state,L,l))
    return all_state


def transin(state, basis, L):
    #This function has two purpose: 
    #    1)check whether state is already in basis, if yes-->s=1, no-->s=0
    #    2)count the period of the translation symmetry
    ntrans = 1
    for bit in range(L):
        s = 0
        statenew = spintrans(state, L, bit+1)
        if (statenew in basis):
            s = 1
            break
        ntrans += 1
    return (s, ntrans)


def find_rep(state,basis):
    idx = -1
    for i in range(len(basis)):
        if(state == basis[i]):
            idx= i
        else:
            continue
    return(idx)


################# create basis states##########################################

def build_all_basis(param_L, Nup): # Used in Real space hamiltonian
    # Using quantum number m to build bases
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


def build_basis(param_L, Nup): # Used in k space hamiltonian
    # Using quantum number Nup to build bases, 
    # and then only save the basis by applying translational symmetry. 
    # For example: 1010 is in the same group of 0101, so we only save one of them.
    L = param_L
    maxdim = 2 ** L
    basis = []
    dim = 0
    for state in range(int(maxdim / 2) + 1):

        n_ones = 0
        for bit in range(L):
            if ((state & (1 << bit)) != 0):  # the bit(th) bit of state is 1
                n_ones += 1
        if (n_ones != Nup):
            continue
        s, ntrans = transin(state, basis, L)
        if (s):
            basis = basis
        else:
            basis.append(state)
            dim += 1
    Nf = []  # normalization factor for each basis
    period = [] #the period of translation, for example 1010 will return the same state be translating the left 2 bits to the right, so the period is 2.
    for i in range(dim):
        s, ntrans = transin(basis[i], basis, L)
        if (ntrans < L):
            normal = 1 / np.sqrt((L / (ntrans)) ** 2 * (ntrans))
            period.append(ntrans)
        else:
            normal = 1 / np.sqrt(L)
            period.append(ntrans)
        Nf.append(normal)
    #print('Basis:')
    #print([bin(basis[i]) for i in range(len(basis))])
    #print(period)
    #print(len(basis))
    print('Dimention of basis:')
    print(dim)
    return (basis,Nf,dim,period)


def build_basis_Nup2(param_L): 
    # Create the basis with Nup=2, Ndn=L-2, Sz=-(L//2-2); using translational symmetry
    L = param_L
    Nup1=2
    maxdim = 2 ** L
    basis = []
    dim = 0
    basis_large = []
    for bit1 in range(L-1):
        for bit2 in range(bit1,L):
            basis_large.append(2**bit1+2**bit2)
    for state in basis_large:
        n_ones = 0
        for bit in range(L):
            if ((state & (1 << bit)) != 0):  # the bit(th) bit of state is 1
                n_ones += 1
        if (n_ones != Nup1):
            continue
        s, ntrans = transin(state, basis, L)
        if (s):
            basis = basis
        else:
            basis.append(state)
            dim += 1
    Nf = []  # normalization factor for each basis
    period = []
    for i in range(dim):
        s, ntrans = transin(basis[i], basis, L)
        if (ntrans < L):
            normal = 1 / np.sqrt((L / (ntrans)) ** 2 * (ntrans))
            period.append(ntrans)
        else:
            normal = 1 / np.sqrt(L)
            period.append(ntrans)
        Nf.append(normal)
    #print('Basis:')
    #print([bin(basis[i]) for i in range(len(basis))])
    print('Dimention:')
    print(dim)
    return (basis,Nf,dim,period)


def total_weight(k,periodR,L):
    if(abs(( 2 * np.pi * (k - 1) / L)*periodR/np.pi/2. - int(( 2 * np.pi * (k - 1) / L)*periodR/np.pi/2.))<1e-6):
        return L/periodR
    else:
        #print(abs(k*periodR/np.pi/2. - int(k*periodR/np.pi/2.)))
        return 0.


def delete_redundant_basis(basis_in,Nf_in,dimbasis,period_in,k,L):
    basis = []
    Nf = []
    period = []
    for nni in range(dimbasis):
        if(period_in[nni]<L and abs(total_weight(k,period_in[nni],L))<1e-6):
            continue
        else:
            basis.append(basis_in[nni])
            Nf.append(Nf_in[nni])
            period.append(period_in[nni])
    dimbasis = len(basis)
    return basis,Nf,dimbasis,period 



if __name__ == "__main__":
    L=8
    Nup=4
    basisq = build_all_basis(L, Nup)
    basis,Nf,dim,period =build_basis(L, Nup)
    


    







