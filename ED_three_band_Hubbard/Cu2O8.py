import numpy as np
import matplotlib.pyplot as plt
from math import cos
from scipy.special import gamma
from matplotlib import cm
import decimal
from hubbard_funcs import *


############## Cu2O8 ## cdup: 1,3,5,7....; cup: -1,-3,-5,... cddn: 2,4,6,8....; cdn: -2,-4,-6...

def create_terms(tpd, tpp, tdd, tppP, Udd, Upp, delta_pd):
    cc_terms = []
    cccc_terms = []
    all_terms = []
    ######### delta_pd ######### 
    listall= [[0,0],[1,1],[3,3],[4,4],[5,5],[6,6],[8,8],[9,9]]
    for li in range(len(listall)):
        list1 = listall[li]
        for i in range(1):
            ##### list1
            # cdup cup
            j1 = 2*(list1[i])
            j2 = 2*(list1[i+1])
            cdup_cddn = j1+1
            cup_cdn = -(j2+1)
            all_terms.append(c_terms([cdup_cddn, cup_cdn],delta_pd))
            
            # cddn cdn
            j1 = 2*(list1[i])+1
            j2 = 2*(list1[i+1])+1
            cdup_cddn = j1+1
            cup_cdn = -(j2+1) 
            all_terms.append(c_terms([cdup_cddn, cup_cdn],delta_pd))
            

    ######### tpd ######### a(b) means "a hop to b"
    listall= [[0,2],[1,2],[2,3],[2,5],[4,7],[7,8],[6,7],[7,9]]
    for li in range(len(listall)):
        list1 = listall[li]
        for i in range(1):
            ##### list1
            # cdup cup
            j1 = 2*(list1[i])
            j2 = 2*(list1[i+1])
            cdup_cddn = j1+1
            cup_cdn = -(j2+1)
            all_terms.append(c_terms([cdup_cddn, cup_cdn],tpd))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],tpd)) #h.c.
            # cddn cdn
            j1 = 2*(list1[i])+1
            j2 = 2*(list1[i+1])+1
            cdup_cddn = j1+1
            cup_cdn = -(j2+1) 
            all_terms.append(c_terms([cdup_cddn, cup_cdn],tpd))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],tpd)) #h.c.
    ######### tpp #########
    listall= [[0,1],[1,5],[3,5],[0,3],[1,4],[5,8],[4,8],[4,6],[6,9],[8,9]]
    for li in range(len(listall)):
        list1 = listall[li]
        for i in range(1):
            ##### list1
            # cdup cup
            j1 = 2*(list1[i])
            j2 = 2*(list1[i+1])
            cdup_cddn = j1+1
            cup_cdn = -(j2+1)
            all_terms.append(c_terms([cdup_cddn, cup_cdn],tpp))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],tpp)) #h.c.
            # cddn cdn
            j1 = 2*(list1[i])+1
            j2 = 2*(list1[i+1])+1
            cdup_cddn = j1+1
            cup_cdn = -(j2+1) 
            all_terms.append(c_terms([cdup_cddn, cup_cdn],tpp))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],tpp)) #h.c.
    ######### tdd #########
    
    ######### tpp' ######### 
    listall= [[0,5],[1,8],[4,9],[1,3],[4,5],[6,8]]
    for li in range(len(listall)):
        list1 = listall[li]
        for i in range(1):
            ##### list1
            # cdup cup
            j1 = 2*(list1[i])
            j2 = 2*(list1[i+1])
            cdup_cddn = j1+1
            cup_cdn = -(j2+1)
            all_terms.append(c_terms([cdup_cddn, cup_cdn],tppP))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],tppP)) #h.c.
            # cddn cdn
            j1 = 2*(list1[i])+1
            j2 = 2*(list1[i+1])+1
            cdup_cddn = j1+1
            cup_cdn = -(j2+1) 
            all_terms.append(c_terms([cdup_cddn, cup_cdn],tppP))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],tppP)) #h.c.
    
    
    for r1 in [2,7]:
        cccc_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],Udd))
        all_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],Udd))
    for r1 in [0,1,3,4,5,6,8,9]:
        cccc_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],Upp))
        all_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],Upp))

   
    return all_terms    





######################## build basis for each k sector
dim, basis = build_basis_useQN(L,1,1)
print("check basis: ",len(basis))
all_term = create_terms(tpd, tpp, tdd, tppP, Udd, Upp, delta_pd)# this is for real space
h_hubbard_mtx = hami_mtx(basis,L,all_term)
eigenvaluehami,eigenmatrix_hami = np.linalg.eigh(h_hubbard_mtx)
singlet0 = eigenvaluehami[0]
#print("singlet energys :", eigenvaluehami)


dim, basis = build_basis_useQN(L,2,0)
print("check basis: ",len(basis))
all_term = create_terms(tpd, tpp, tdd, tppP, Udd, Upp, delta_pd)# this is for real space
h_hubbard_mtx = hami_mtx(basis,L,all_term)
eigenvaluehami,eigenmatrix_hami = np.linalg.eigh(h_hubbard_mtx)
triplet0 = eigenvaluehami[0]
print("2 holes energys :", eigenvaluehami[0])


print("estimated Jnnn: ",triplet0 - singlet0)


dim, basis = build_basis_useQN(L,2,1)
print("check basis: ",len(basis))
all_term = create_terms(tpd, tpp, tdd, tppP, Udd, Upp, delta_pd)# this is for real space
h_hubbard_mtx = hami_mtx(basis,L,all_term)
eigenvaluehami,eigenmatrix_hami = np.linalg.eigh(h_hubbard_mtx)
singlet0 = eigenvaluehami[0]
print("3 holes energys :", eigenvaluehami[0])
print("estimated t': ", (eigenvaluehami[1] - eigenvaluehami[0])/2.0)




'''
eivec = np.transpose(eigenmatrix_hami)
gs = eivec[0]
excited = eivec[1]
for i in range(len(gs)):
    if(abs(gs[i])>1e-4):
        print("gs: ",bin(basis[i]), gs[i])
    if(abs(excited[i])>1e-4):
        print("excited: ",bin(basis[i]), excited[i])
'''
