import numpy as np
import matplotlib.pyplot as plt
from hubbard_k_funcs import *

############## Hubbard PBC
'''
cc_terms = []
list_n = np.zeros(int(L*2))
cccc_terms = []
all_terms = []
for i in range(int(L*2)):
    list_n[i] += t
    cdup_cddn = i+1
    if(i+3<=L*2):
        cup_cdn = -(i+3)
    else:
        cup_cdn = -((i+3+L*10)%(L*2))
    cc_terms.append(c_terms([cdup_cddn, cup_cdn],list_n[i]))
    cc_terms.append(c_terms([-cup_cdn,-cdup_cddn],list_n[i]))
    print(" cc terms ",[cdup_cddn, cup_cdn],list_n[i])
    print(" cc terms ",[-cup_cdn,-cdup_cddn],list_n[i])
    all_terms.append(c_terms([cdup_cddn, cup_cdn],list_n[i]))
    all_terms.append(c_terms([-cup_cdn,-cdup_cddn],list_n[i]))
'''
############## Hubbard OBC

def create_terms_lri(lamda):
    cc_terms = []
    list_n = np.zeros(int(L*2))
    cccc_terms = []
    all_terms = []
    for i in range(int(L*2)):
        list_n[i] += t
        cdup_cddn = i+1
        if(i+3<=L*2):
            cup_cdn = -(i+3)   
            cc_terms.append(c_terms([cdup_cddn, cup_cdn],list_n[i]))
            cc_terms.append(c_terms([-cup_cdn,-cdup_cddn],list_n[i]))
            print(" cc terms ",[cdup_cddn, cup_cdn],list_n[i])
            print(" cc terms ",[-cup_cdn,-cdup_cddn],list_n[i])
            all_terms.append(c_terms([cdup_cddn, cup_cdn],list_n[i]))
            all_terms.append(c_terms([-cup_cdn,-cdup_cddn],list_n[i]))
    
    
    for r1 in range(L):
        cccc_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],U))
        all_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r1 + 2, -(2 * r1 + 2)],U))

    for r1 in range(L-1):
        r2=r1+1
        all_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r2 + 1, -(2 * r2 + 1)],V))
        all_terms.append(c_terms([2 * r1+1, -(2 * r1+1), 2 * r2 + 2, -(2 * r2 + 2)],V))
        all_terms.append(c_terms([2 * r1+2, -(2 * r1+2), 2 * r2 + 1, -(2 * r2 + 1)],V))
        all_terms.append(c_terms([2 * r1+2, -(2 * r1+2), 2 * r2 + 2, -(2 * r2 + 2)],V))


    for i in range(L-1):
        for j in range(i+1,L):
            if(i==j): continue
            coef = ((-1)**abs(j-i+1))*lamda/(j-i)**2
            print(i,j,coef)
            all_terms.append(c_terms([2 * i+1, -(2 * i+2), 2 * j + 2, -(2 * j + 1)],coef*0.5))#S+iS-j
            all_terms.append(c_terms([2 * i+2, -(2 * i+1), 2 * j + 1, -(2 * j + 2)],coef*0.5))#S-iS+j
            all_terms.append(c_terms([2 * i+1, -(2 * i+1), 2 * j + 1, -(2 * j + 1)],coef*0.25))#NupiNupj
            all_terms.append(c_terms([2 * i+1, -(2 * i+1), 2 * j + 2, -(2 * j + 2)],coef*(-0.25)))#NupiNdnj
            all_terms.append(c_terms([2 * i+2, -(2 * i+2), 2 * j + 1, -(2 * j + 1)],coef*(-0.25)))#NdniNupj
            all_terms.append(c_terms([2 * i+2, -(2 * i+2), 2 * j + 2, -(2 * j + 2)],coef*0.25))#NdniNdnj
    return all_terms    


        


basis_hubbard = []
dimr,basisr = build_basis_useQN(L,int(L/2),int(L/2))#a spin up #int(N/2),int(N/2))
for i in basisr:
    basis_hubbard.append(i)
print([bin(basis_hubbard[i]) for i in range(len(basis_hubbard))])


f1= open("overlap_L=%d_a=%.1f_U=%.1f.txt" % (L, 2.0,U),"w+")

for l in range(50):
    lamda=l*0.2-10.
    all_term = create_terms_lri(lamda)
    h_hubbard_mtx = hami_mtx(basis_hubbard,L,all_term)
    #print(h_hubbard_mtx)
    values,vectors = np.linalg.eigh(h_hubbard_mtx)
    print(values)#,np.transpose(vectors)[1])
    vec = np.transpose(vectors)[23]
    for m in range(len(basis_hubbard)):
        if(abs(vec[m])>0.00000001): print(bin(basis_hubbard[m]), vec[m],eta_or_not(basis_hubbard[m],L))
    vectors = np.transpose(vectors)

    overlap = np.zeros(len(values),dtype=dt)
    for im in range(len(values)):
        for iba in range(len(basis_hubbard)):
            if(eta_or_not(basis_hubbard[iba],L)==1):
                overlap[im] += vectors[im][iba] * np.conjugate(vectors[im][iba])
    print(overlap)
    for i in range(len(basis_hubbard)):
        plt.scatter(lamda,values[i],s=50.*overlap[i],color='blue')
        f1.write("lamda= %.1f %d th eigenstate= %.8f overlap= %.8f\r\n" % (lamda, i, values[i],overlap[i]))
f1.close()
plt.xlabel('$\lambda$',fontsize=18)
plt.ylabel('E',fontsize=18)
plt.title('L=%d, alpha=2, U=%.1f' % (L,U))
plt.show()
            



