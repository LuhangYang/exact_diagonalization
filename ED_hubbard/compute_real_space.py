import numpy as np
from hubbard_k_funcs import *

xBC = BoundaryCondition.OBC
#yBC = BoundaryCondition.PBC

############## Hubbard OBC

def create_terms_lri_2D():
    cc_terms = []
    cccc_terms = []
    all_terms = []
  
    for iy in range(Ly-1):
        for ix in range(Lx-1):
            # cdup cup
            i = 2*(iy*Lx+ix)    
            #list_n[i] += t
            cdup_cddn = i+1
            #if(i+3<=L*2):
            cup_cdn = -(i+1+2) # hoping to the right
            all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
            cup_cdn = -(i+1+2*Lx) #hoping to the top
            all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
            # cddn cdn
            i = 2*(iy*Lx+ix)+1   
            #list_n[i] += t
            cdup_cddn = i+1
            #if(i+3<=L*2):
            cup_cdn = -(i+1+2) # hoping to the right
            all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
            cup_cdn = -(i+1+2*Lx) #hoping to the top
            all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
    #bottom row:
    for ix in range(Lx-1):
        # cdup cup
        i = 2*(Lx*(Ly-1)+ix)    
        #list_n[i] += t
        cdup_cddn = i+1
        cup_cdn = -(i+1+2) # hoping to the right
        all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
        all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
        # cddn cdn
        i = 2*(Lx*(Ly-1)+ix)+1   
        #list_n[i] += t
        cdup_cddn = i+1
        cup_cdn = -(i+1+2) # hoping to the right
        all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
        all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
   
    #right end column:
    for iy in range(Ly-1):
        # cdup cup
        i = 2*(iy*Lx+Lx-1)    
        #list_n[i] += t
        cdup_cddn = i+1
        cup_cdn = -(i+1+2*Lx) #hoping to the bottom
        all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
        all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
        
        # cddn cdn
        i = 2*(iy*Lx+Lx-1)+1   
        #list_n[i] += t
        cdup_cddn = i+1
        cup_cdn = -(i+1+2*Lx) #hoping to the bottom
        all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
        all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
        
    
    
    if(xBC == BoundaryCondition.PBC):
        for iy in range(Ly):
            i = 2*(iy*Lx+Lx-1)
            #list_n[i] += t
            cdup_cddn = i+1
            cup_cdn = -(2*iy*Lx+1) # hoping to the right
            all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
           
            i = 2*(iy*Lx+Lx-1)+1
            #list_n[i] += t
            cdup_cddn = i+1
            cup_cdn = -(2*iy*Lx+2) # hoping to the right
            all_terms.append(c_terms([cdup_cddn, cup_cdn],t))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],t)) #h.c.
            
    

    for r1 in range(L):
        cccc_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],U))
        all_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],U))

   
    return all_terms    





######################## build basis for each k sector
dim, basis = build_basis_useQN(Lx*Ly,int(Lx*Ly/2),int(Lx*Ly/2)) #basis = build_HD_basis(L,True)
print("check flag 2: ",len(basis))

#print("energy ",eigenvaluehami)
## op == 1: Sz
## op == 2: j

################# creat all the terms ####### when using the H in k-space, we automatically assume pbc
all_term = create_terms_lri_2D()# this is for real space
print(all_term)
###############3# calculate eigen energy, eigen vectors, and H matrix
h_hubbard_mtx = hami_mtx(basis,L,all_term)
#for i in range(len(h_hubbard_mtx)):
#    for j in range(len(h_hubbard_mtx)):
#        if(h_hubbard_mtx[i][j]!=0):
#            print(i,j,bin(basis[i]),bin(basis[j]),h_hubbard_mtx[i][j])
print(h_hubbard_mtx)
eigenvaluehami,eigenmatrix_hami = np.linalg.eigh(h_hubbard_mtx)
print("energys :", eigenvaluehami)
#print(np.linalg.inv(eigenmatrix_hami))


################################ Compute correction vector in r-space ##############################


gs = eigenmatrix_hami.transpose()[0]*1.0







