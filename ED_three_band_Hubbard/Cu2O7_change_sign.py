import numpy as np
import matplotlib.pyplot as plt
from math import cos
from scipy.special import gamma
from matplotlib import cm
import decimal
from hubbard_funcs import *


############## Cu2O7 ## cdup: 1,3,5,7....; cup: -1,-3,-5,... cddn: 2,4,6,8....; cdn: -2,-4,-6...

def create_terms(tpd, tpp, tdd, tppP, Udd, Upp, delta_pd):
    cc_terms = []
    cccc_terms = []
    all_terms = []
    ######### delta_pd ######### 
    listall= [[0,0],[1,1],[3,3],[4,4],[5,5],[7,7],[8,8]]
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
    listall= [[0,2],[1,2],[5,6],[4,6]]
    for li in range(len(listall)):
        list1 = listall[li]
        for i in range(1):
            ##### list1
            # cdup cup
            j1 = 2*(list1[i])
            j2 = 2*(list1[i+1])
            cdup_cddn = j1+1
            cup_cdn = -(j2+1)
            all_terms.append(c_terms([cdup_cddn, cup_cdn],-tpd))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],-tpd)) #h.c.
            # cddn cdn
            j1 = 2*(list1[i])+1
            j2 = 2*(list1[i+1])+1
            cdup_cddn = j1+1
            cup_cdn = -(j2+1) 
            all_terms.append(c_terms([cdup_cddn, cup_cdn],-tpd))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],-tpd)) #h.c.
    listall= [[2,3],[2,4],[6,8],[6,7]]
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
    listall= [[0,1],[3,4],[4,5],[7,8]]
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
    listall= [[1,4],[0,3],[4,7],[5,8]]
    for li in range(len(listall)):
        list1 = listall[li]
        for i in range(1):
            ##### list1
            # cdup cup
            j1 = 2*(list1[i])
            j2 = 2*(list1[i+1])
            cdup_cddn = j1+1
            cup_cdn = -(j2+1)
            all_terms.append(c_terms([cdup_cddn, cup_cdn],-tpp))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],-tpp)) #h.c.
            # cddn cdn
            j1 = 2*(list1[i])+1
            j2 = 2*(list1[i+1])+1
            cdup_cddn = j1+1
            cup_cdn = -(j2+1) 
            all_terms.append(c_terms([cdup_cddn, cup_cdn],-tpp))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],-tpp)) #h.c.
    ######### tdd #########
    # 2(6)
    # cdup cup
    j1 = 2*(2)
    j2 = 2*(6)
    cdup_cddn = j1+1
    cup_cdn = -(j2+1)
    all_terms.append(c_terms([cdup_cddn, cup_cdn],tdd))
    all_terms.append(c_terms([-cup_cdn,-cdup_cddn],tdd)) #h.c.
    # cddn cdn
    j1 = 2*(2)+1
    j2 = 2*(6)+1
    cdup_cddn = j1+1
    cup_cdn = -(j2+1) 
    all_terms.append(c_terms([cdup_cddn, cup_cdn],tdd))
    all_terms.append(c_terms([-cup_cdn,-cdup_cddn],tdd)) #h.c.
    ######### tpp' ######### 
    listall= [[1,3],[5,7],[1,5],[0,4],[4,8],[3,7]]
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
    
    
    for r1 in [2,6]:
        cccc_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],Udd))
        all_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],Udd))
    for r1 in [0,1,3,4,5,7,8]:
        cccc_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],Upp))
        all_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],Upp))

   
    return all_terms    

def build_basis_3holes(L):
    dim1, basis1 = build_basis_useQN(5,1,1)
    print("check basis1: ",len(basis1))
    dim2, basis2 = build_basis_useQN(4,1,0)
    print("check basis2: ",len(basis2))
    left_basis = []
    right_basis = []
    for i1 in basis1:
        for i2 in basis2:
            if((IBITS(i1,2*2)==1 and IBITS(i1,2*2+1)==0) or (IBITS(i1,2*2)==0 and IBITS(i1,2*2+1)==1)):
                left_basis.append(i1*(2**8)+i2)
                right_basis.append(i2*(2**10)+i1)
            #print(bin(i1), bin(i2), bin(i1*(2**8)+i2), bin(i2*(2**10)+i1))
    print(len(left_basis))
    norm_factor = 1.0/np.sqrt(len(left_basis)*2)
    print(norm_factor)
    coef_left_E = [norm_factor for ie in range(len(left_basis))]
    coef_right_E = [norm_factor for ie in range(len(right_basis))]
    coef_left_O = [-norm_factor for ie in range(len(left_basis))]
    coef_right_O = [norm_factor for ie in range(len(right_basis))]    
    all_basis = [left_basis[ia] for ia in range(len(left_basis))]
    coef_E = [coef_left_E[ia]  for ia in range(len(left_basis))]
    coef_O = [coef_left_O[ia]  for ia in range(len(left_basis))]
    for ib in range(len(right_basis)):
        all_basis.append(right_basis[ib])
        coef_E.append(coef_right_E[ib])
        coef_O.append(coef_right_O[ib])
    return all_basis, coef_E,coef_O


######################## build basis for each k sector
dim, basis = build_basis_useQN(L,1,1)
print("check basis: ",len(basis))
all_term = create_terms(tpd, tpp, tdd, tppP, Udd, Upp, delta_pd)# this is for real space
h_hubbard_mtx = hami_mtx(basis,L,all_term)
eigenvaluehami,eigenmatrix_hami = np.linalg.eigh(h_hubbard_mtx)
singlet0 = eigenvaluehami[0]
print("singlet energys :", eigenvaluehami[0])


dim, basis = build_basis_useQN(L,2,0)
print("check basis: ",len(basis))
all_term = create_terms(tpd, tpp, tdd, tppP, Udd, Upp, delta_pd)# this is for real space
h_hubbard_mtx = hami_mtx(basis,L,all_term)
eigenvaluehami,eigenmatrix_hami = np.linalg.eigh(h_hubbard_mtx)
triplet0 = eigenvaluehami[0]
print("triplet energys :", eigenvaluehami[0])


print("estimated Jnn: ",triplet0 - singlet0)


dim, basis = build_basis_useQN(L,2,1)#
all_term = create_terms(tpd, tpp, tdd, tppP, Udd, Upp, delta_pd)# this is for real space
h_hubbard_mtx = hami_mtx(basis,L,all_term)
eigenvaluehami,eigenmatrix_hami = np.linalg.eigh(h_hubbard_mtx)
singlet0 = eigenvaluehami[0]
#print("3 holes energys :", eigenvaluehami)
print("estimated t: ", -(eigenvaluehami[1] - eigenvaluehami[0])/2.0, eigenvaluehami[1], eigenvaluehami[0])
'''
print("##########################################################")
##### Directly build the states
basis,coef_E, coef_O = build_basis_3holes(L)
print("check basis: ",len(basis), np.linalg.norm(coef_E))
EE = hami_energy(basis,coef_E,L,all_term)#
EO = hami_energy(basis,coef_O,L,all_term)#
print("estimated t: ", -(EE - EO)/2.0, EE,EO)
'''





