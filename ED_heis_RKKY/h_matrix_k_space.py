from basis import *
from hamiltonian_operator import *
import numpy as np
from math import cos
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigs

######################### DEFINE GLOBAL VARIABLES HERE ###################
dt = np.dtype(np.complex128)


################# commonly used functions #################################
def IBITS(n, i):
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

def ch(vector,k): #kth subblock ranges from 1 to L, corresponding to [0,2pi)
    thetak = 2 * np.pi * (k - 1) / L
    c = np.exp(-1j * (vector) * thetak)
    return c

def find_rep(state,basis):
    idx = -1
    for i in range(len(basis)):
        if(state == basis[i]):
            idx= i
        else:
            continue
    return idx



################# Solve H matrix ##############################################
def total_weight(k,periodR,L):
    if(abs(( 2 * np.pi * (k - 1) / L)*periodR/np.pi/2. - int(( 2 * np.pi * (k - 1) / L)*periodR/np.pi/2.))<1e-6):
        return L/periodR
    else:
        #print(abs(k*periodR/np.pi/2. - int(k*periodR/np.pi/2.)))
        return 0.
def sum_szsz_pbc(s_in,L):
    su = 0.
    for isite in range(L):
        i1site = (isite+1)%L
        if(IBITS(s_in,isite)==IBITS(s_in,i1site)):
            su += 0.25
        else:
            su += -0.25
    return su


def hami_mtx_heis(basis_in,Nf_in,dimbasis,period_in,k,L):
    # This is conventional Heisenberg model
    nvector=[1,L]
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
    hmatrix = np.zeros(shape = (dimbasis,dimbasis),dtype=dt)
    for term in range(2): # 0: SzSz // 1: Spin flip
        for i in range(dimbasis):
            out_state, x = hami_heis(basis[i], term,L) ##apply hami on |a0>
            if(term==0):
                hmatrix[i, i] = hmatrix[i, i] + sum_szsz_pbc(basis[i],L)
            else:
                for bi in range(len(out_state)):  
                    find_bi = 0
                    for vector in range(nvector[term]):
                        new_x = spintrans(out_state[bi],L,vector) #trasx(bi, vector) from 0 to L-1
                        if((new_x in basis) and find_bi==0):
                            idx = find_rep(new_x, basis)
                            if (abs(total_weight(k,period[i],L))>1e-6 and abs(total_weight(k,period[idx],L))>1e-6):
                                coef = np.exp(-1j*vector*((2*np.pi*(k-1)/L)))#ch(vector,k)
                                hmatrix[i,idx]=hmatrix[i,idx]+x[bi]*coef*np.sqrt(period[i]/period[idx])
                                find_bi=1
    values,vectors = np.linalg.eigh(hmatrix)
    return values,vectors


def hami_mtx_heis_rkky(basis_in,Nf_in,dimbasis,period_in,k,L,lamda,alpha):
    # This is Heisenberg model with RKKY type interactions
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
    print("dim after : ", len(basis))
    nvector = [1, L]
    hmatrix = np.zeros(shape=(len(basis), len(basis)), dtype=dt)
    
    for term in range(2): # 0: SzSz // 1: Spin flip
        for i in range(len(basis)):
            coef = 0.
            for vector in range(nvector[term]):
                new_x = spintrans(basis[i],L,vector) 
                out_state, x = hami_heis_rkky(new_x, term,lamda,alpha,L)
                if (term == 0):
                    hmatrix[i, i] = hmatrix[i, i] + x
                if (term == 1):
                    for xi in range(len(x)):
                        if (np.abs(x[xi]) > 1.e-10):
                            idx = find_rep(out_state[xi], basis)
                            if (idx >= i):
                                coef = ch(vector,k) #e^-i*vector*k, kth subblock
                                c2=0.
                                if(period[idx] == L):
                                    c2=1
                                else:
                                    for p in range(int(L/period[idx])):
                                        thetak = 2 * np.pi * (k - 1) / L
                                        c2 += np.exp(1j * p*(period[idx]) * thetak)
                                if(np.abs(c2)<1.e-7 or np.abs(coef)<1.e-7):
                                    hmatrix = hmatrix
                                else:
                                    hmatrix[i,idx]=hmatrix[i,idx]+(x[xi])*coef/c2*(Nf[i]/Nf[idx]) #*L
                                hmatrix[idx,i]=np.conjugate(hmatrix[i,idx])

    values,vectors = np.linalg.eigh(hmatrix)
    return values,vectors



if __name__ == "__main__":
    L=8
    Nup=4
    basis,Nf,dim,period =build_basis(L, Nup)
    vals_rkky = []
    vals_regular = []
    for k in range(L):
        va,vb = hami_mtx_heis_rkky(basis,Nf,dim,period,k+1,L,0.,2.)
        vals_rkky.append(va)
        #print("k: ",k+1,[va[i] for i in range(len(va))])
    for k in range(L):
        values,vectors = hami_mtx_heis(basis,Nf,dim,period,k+1,L)
        vals_regular.append(values)
        #print("k: ",k+1,[values[i] for i in range(len(values))])
    #for i in range(len(vals_rkky)):
    #    for j in range(len(vals_rkky[i])):
    #        if(abs(vals_rkky[i][j]-vals_regular[i][j])>1e-5): print("Wrong!!")
    







