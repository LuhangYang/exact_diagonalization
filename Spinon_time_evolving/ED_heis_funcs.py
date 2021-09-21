import numpy as np
import matplotlib.pyplot as plt
from math import cos
from scipy.special import gamma
from matplotlib import cm
import decimal
import Permutation2


L=48
num_dis_max = int((L-4.)/2. + 1.)
alpha_all=[1.8,2.0,2.2,2.4]
dt = np.dtype(np.complex128)

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

def move_spinons_basis_obc(L):
    move_basis = []
    spinons = 1 + 2 + 4
    for xbit in range(int((L - 4) / 2)):
        spinons += 2 ** (3 + xbit * 2 + 1)
    print("first spinon: ",spinons,bin(spinons))
    pos_spinon = [0, 1, 2]
    move_basis.append(spinons)
    spinons = 1 + 2 + 8+16
    for xbit in range(int((L - 6) / 2)):
        spinons += 2 ** (5 + xbit * 2 + 1)
    print("second spinon: ", spinons, bin(spinons))
    pos_spinon = [0, 1, 2]
    move_basis.append(spinons)
    for ith in range(int((L - 4) / 2.)):
        new_spinons = 1+2
        for ii in range(ith+1):
            new_spinons += 2**(2+ii*2+1)
        new_spinons += 2**(5+ith*2)+2**(6+ith*2)
        for jj in range(int((L-6-ith*2-1)/2)):
            new_spinons += 2**(6+ith*2+1+jj*2+1)
        if(new_spinons not in move_basis):
            all_trans = build_all_trans(new_spinons, L)
            notin = 0
            for y in all_trans:
                if(y in move_basis):
                    notin = 1
            notin =0 #pbc: comment out
            if(notin==0):
                move_basis.append(new_spinons)
                print(new_spinons,bin(new_spinons))
                
    return move_basis

def IthBITinN(n,i):
    return ((n >> i) & 1)


def one_anti(state):
    num_anti = 0
    for m in range(L-1):
        if(IthBITinN(state,m)==0 and IthBITinN(state,m+1)==0):
            num_anti +=1
    return num_anti


def one_spinon(state):
    num_spinon = 0
    for m in range(L-1):
        if(IthBITinN(state,m)==1 and IthBITinN(state,m+1)==1):
            num_spinon +=1
    return num_spinon



def find_rep(state,basis):
    idx = -1
    for i in range(len(basis)):
        if(state == basis[i]):
            idx= i
        else:
            continue
    return idx



def hami(state_in, term,lamda,alpha): # H operator act on state_in,term=0:diag,term=1:spin flip
    if(term == 0):
        # diagonal term
        x = 0.
        state_out = state_in
        for site_i in range(L-1):
            for site_j in range(site_i+1,L):
                delta = site_j - site_i
                if(delta>4): continue
                ijnotsame = ((IthBITinN(state_in, site_i) ^ IthBITinN(state_in, site_j)))
                if(delta==1):  # (sitei-1+L)mod L,site_j is the right side bit on site_i
                    x += ijnotsame * (-0.25) + (1 - ijnotsame) * (0.25)
                else:
                    x += (ijnotsame * (-0.25) + (1 - ijnotsame) * (0.25)) * (-lamda*((-1)**delta)/(delta**alpha))
    elif(term==1):
        # off-diagonal term
        # do the spin flip at position x
        state_out = []
        x = []
        for site_i in range(L-1):
            for site_j in range(site_i+1,L):
                delta = site_j - site_i
                if(delta>4): continue
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
    return (state_out,x) 

def hami_mtx(basis,L,lamda,alpha):
    nvector = [1, L]
    hmatrix = np.zeros(shape=(len(basis), len(basis)), dtype=dt)
    for term in range(2): # 0: SzSz // 1: Spin flip
        for i in range(len(basis)):
            out_state, coef_x = hami(basis[i], term,lamda,alpha)
            if (term == 0):
                hmatrix[i, i] = hmatrix[i, i] + coef_x
            if (term == 1):
                for xi in range(len(coef_x)):
                    if (np.abs(coef_x[xi]) > 1.e-10):
                        idx = find_rep(out_state[xi], basis)
                        if (idx >= i):
                            hmatrix[i,idx]=hmatrix[i,idx]+coef_x[xi]
                            hmatrix[idx,i]=np.conjugate(hmatrix[i,idx])

    values,vectors = np.linalg.eigh(hmatrix)
    eigenv = values
    eigenm = vectors
    return values,vectors


