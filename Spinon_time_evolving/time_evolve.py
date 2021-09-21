import numpy as np
import matplotlib.pyplot as plt
from math import cos
from scipy.special import gamma
from matplotlib import cm
import decimal
import Permutation2
from ED_heis_funcs import *


basis = []

get_basis = move_spinons_basis_obc(L)
len_basis = len(get_basis)
print("length of basis: ", len_basis)


basis.append(spintrans(get_basis[0],L,int(L/2.-1)))
print(basis[0],bin(basis[0]))


res = Permutation2.permutation(L, 2, 0)
for i in range(len(res)):
    if(int(res[i],2) % 2 ==0):
        basis.append(int(res[i],2))
print(len(basis))

def one_magnon(state):
    num_magnon = 0
    for m in range(L-2):
        if(IthBITinN(state,m)==1 and IthBITinN(state,m+1)==1 and IthBITinN(state,m+2)==1):
            num_magnon +=1
    return num_magnon

res = Permutation2.permutation(L, 3, 1)
for i in range(len(res)):
    if(int(res[i],2) % 2 ==0):
        if(one_magnon(int(res[i],2)) > 0):
            basis.append(int(res[i],2))
print(len(basis))


for yy in range(1,len(basis)-1):
    if(basis[yy] == basis[0]):
        del basis[yy]
print("The length of the basis: ", len(basis))
#print([bin(basis[i]) for i in range(len(basis))])





na=0
for alpha in alpha_all:
    eigenvaluehami = hami_mtx(basis,L,1.0,alpha)[0]
    eigenmatrix_hami = hami_mtx(basis,L,1.0,alpha)[1]
    print("energys :", eigenvaluehami)
    print(np.linalg.inv(eigenmatrix_hami))
    b0 = np.zeros(shape=(len(basis), 1), dtype=dt)
    b0[0][0] = 1.0
    #print("b0 ",b0)

    t_interval = 0.1
    
    exph = np.zeros(shape=(len(basis),len(basis)),dtype=dt)
    for i in range(len(basis)):
        exph[i,i] = np.exp(eigenvaluehami[i]*-1J*t_interval)

    a0 = np.matmul(exph,np.linalg.inv(eigenmatrix_hami))
    exph = np.matmul(eigenmatrix_hami,a0)
    
    #print("a0 ",a0)
    #print("check matrix solve: ",np.matmul(eigenmatrix_hami,a0))
    
    nn_all = []
    nnn_all = []
    szsz_all = []
    max_nn = []
    szt0 = np.zeros(shape=L, dtype=dt)
    
    f1= open("nn_L=%d_a=%.1f.txt" % (L, alpha),"w+")
    f2= open("nnn_L=%d_a=%.1f.txt" % (L, alpha),"w+")
    #f3= open("szsz_L=%d_a=%.1f.txt" % (L, alpha),"w+")
    for ts in range(400):
        tf = ts*t_interval
        b0 = np.matmul(exph,b0)
        nnt = np.zeros(shape=(L-1), dtype=dt)
        szszt = np.zeros(shape=(L-1), dtype=dt)
        nnnt = np.zeros(shape=(L-2), dtype=dt)
        
        for l in range(L-1):
            for i in range(len(basis)):
                 if(IthBITinN(basis[i],l) == 1 and IthBITinN(basis[i],l+1) == 1):
                     nnt[l] +=  b0[i]*np.conjugate(b0[i])
            f1.write("n(i)n(i+1) t=%.1f L=%d %.8f\r\n" % (tf, l, nnt[l].real))
            nnt[l] = nnt[l].real
        for l in range(L-2):
            for i in range(len(basis)):
                 if(IthBITinN(basis[i],l) == 1 and IthBITinN(basis[i],l+1) == 1 and IthBITinN(basis[i],l+2) == 1):
                     nnnt[l] +=  b0[i]*np.conjugate(b0[i])
                     #print("l: ",l,"i: ",i,basis[i],nnnt[l])
            f2.write("n(i)n(i+1)n(i+2) t=%.1f L=%d %.8f\r\n" % (tf, l, nnnt[l].real))
            nnnt[l] = nnnt[l].real

        max_nn.append(max(nnt))
        nn_all.append(nnt)
        nnn_all.append(nnnt)
    f1.close()
    f2.close()
    na+=1




