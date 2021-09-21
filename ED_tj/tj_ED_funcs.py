import numpy as np
import matplotlib.pyplot as plt
from math import cos
#from math import comb #only exist in python 3.8+
from scipy.special import gamma
from matplotlib import cm
import decimal
import Permutation2
from itertools import combinations


# fermion sites, always count from right to left
#L>=4
#L=7
#num_dis_max = int((L-4.)/2. + 1.)
#alpha_all=[1.8]
dt = np.dtype(np.complex128)

###########################################################
###########  Build  basis  ###############################
##########################################################
def IBITS(n,i): #count from right site from 0.
    return ((n >> i) & 1)

def rSubset(arr, r):
    # return list of all subsets of length r
    # to deal with duplicate subsets use 
    # set(list(combinations(arr, r)))
    return list(combinations(arr, r))

def Nup_num_fermi(sttate,Li):
    numr = 0
    for si in range(Li):
        if(IBITS(sttate,2*si)):
            numr +=1
    return numr

def Nup_num_heis(sttate,Li):
    numr = 0
    for si in range(Li):
        if(IBITS(sttate,si)):
            numr +=1
    return numr

def spinon_num_fermi(state,Li):#spinon means (dn dn) in this case
    numr = 0
    for si in range(Li-1):
        if(IBITS(state,2*si+1) and IBITS(state,2*(si+1)+1)):
            numr += 1
    return numr

def anti_spinon_num_fermi(state,Li):#anti spinon means (up up) in this case
    numr = 0
    for si in range(Li-1):
        if(IBITS(state,2*si) and IBITS(state,2*(si+1))):
            numr += 1
    return numr
 

def build_basis_useQN(L1,Nup1,Ndn1): #count site from right to left; NO DOUBLE OCCUPANCY for tJ!!!!!!!
    # Choose #Nup sites for spin up, multiply by #Ndn sites for spin dn, the maximum is L!/(L/2)!
    basis1 = []
    site_list = [i for i in range(L1)]
    Nup_site_list = rSubset(site_list,Nup1)
    #print(Nup_site_list)
    for nupl in Nup_site_list:
        #print(nupl)
        site_list_rest = [ele for ele in site_list if ele not in nupl]
        #print(site_list_rest)
        Ndn_site_list = rSubset(site_list_rest,Ndn1)
        for ndnl in Ndn_site_list:
            #print(nupl,ndnl)
            state=0.
            for siteup in range(Nup1):
                state += 2**(2*nupl[siteup])
            for sitedn in range(Ndn1):
                state += 2**(2*ndnl[sitedn]+1)
            basis1.append(int(state))  
            #print(int(state),bin(int(state)))
    dim1 = len(basis1)
    #print ("Dim=",dim1)
    #print ("Basis:")
    #print ([bin(basis1[i]) for i in range(dim1)])
    return dim1,basis1

def basis_spinon_antispinon_hole(Nspinon,Nanti,Li):
    dim0,basis0 = build_basis_useQN(Li,Li//2-1,Li//2)
    basis_re = []
    for basi in basis0:
        if((spinon_num_fermi(basi,Li) == Nspinon) and (anti_spinon_num_fermi(basi,Li) == Nanti)):
            basis_re.append(basi)
    return basis_re


def spintrans(state,L,bit): #translate the left #n bit to the right
    b = state >> L-bit
    lastbit=0
    for i in range(bit-1):
        lastbit += 2**i
    state = state << bit
    state = state & (2**L-1)
    state = state | b
    return state

#print(spintrans(123,8,2))
def build_all_trans(state,L):
    all_state = []
    for l in range(L):
        all_state.append(spintrans(state,L,l))
    return all_state

def first_spin_hole_basis_obc(L):## 2*n means the nth site with spin up, 2*n+1 means the nth site with spin dn, sites count from right. At beginning it's ...up dn (up dn up dn). Totally 2L bits.
    move_basis = []
    
    spinons = 2**(2*0+1) + 2**(2*2+1) +2**(2*3) 
    for xbit in range(4,L):
        if(xbit%2==0):#spin dn
            spinons += 2 ** (2*xbit + 1)
        else:#spinup
            spinons += 2 ** (2*xbit)
    print("first state: ",spinons,bin(spinons))
    move_basis.append(spinons)
    for trans in range(1,L-2):
            new_state = spintrans(spinons,2*L,2*trans)
            if(new_state not in move_basis): 
                move_basis.append(new_state)
                print(new_state,trans,bin(new_state))
    
    return move_basis

########################################################################################
###########  Hamiltonian operator and Hamiltonian matrix  ##############################
#######################################################################################

def IthBITinN(n,i):
    return ((n >> i) & 1)

def indexinlist(List,num):
    for x in range(len(List)):
        if(List[x]==num or List[x]==-num):
            return x

def swapPositions(list, r, l): #move pos1(r) element to l(pos2)
    new_list = []
    for i in range (0,l):
        new_list.append(list[i])
    new_list.append(list[r])
    for i in range(l,r):
        new_list.append(list[i])
    for i in range(r+1,len(list)):
        new_list.append(list[i])
    return new_list

def fermion_sign(lst,state):
    lstnew = [lst[0],-lst[1],lst[2],-lst[3]]
    lstnew.sort(reverse=True)  #descenting order
    #print("new list: ",lstnew)
    indexnew = [indexinlist(lst,lstnew[0]),indexinlist(lst,lstnew[1]),indexinlist(lst,lstnew[2]),indexinlist(lst,lstnew[3])]
    #print(indexnew)
    sign=1.0
    for i in range(3):
        sign *= (-1)**(indexinlist(indexnew,i)-i)
        indexnew = swapPositions(indexnew,indexinlist(indexnew,i),i)
        #print("i: ",i, sign, indexnew)
    lstn = [lstnew[xi]-1 for xi in range(4)]
    for ii in range(2*L):
        if(IBITS(state,2*L-ii)): lstn.append(2*L-ii)
        #print(lstn)
    for i in ([3,2,1,0]):
        for j in range(4,len(lstn)):
            if(lstn[i]< lstn[j]): 
                sign *= -1
                
    return sign
    
#print("check fermion sign: ", fermion_sign([1, -3, 2, -12],228))#fermion_sign([3, -7, 2, -6]))

def hami_obc(state_in, term,t,J,lamda,alpha,L):#obc  ##term=1,2-->diagonal, tem=0,3-->off-diagonal
    # H operator act on state_in, 
    #  term=0: -t*(Cd_iC_i+1+hc); #off-diagonal
    #  term=1: -J/4*(N_i N_i+1);  #diagonal
    #  term=2: J*(Sz_i Sz_i+1) #diagonal
    #  term=3: J*(S+_i S-+i+1 + hc) #off-fiagonal
    #  term=4: lamda*(Sz_i Sz_i+delta)*(-1)^(delta+1)/abs(delta)^alpha, sum over delta=[2,L-1] for obc
    #  term=5: lamda*(S+_i S-_i+delta + hc)*(-1)^(delta+1)/abs(delta)^alpha
    if(term==0):
        # off-diagonal term, Cdup_i*Cup_i+1 - Cup_i*Cdup_i+1,Cddn_i Cdn_i+1,and hc., for all "i"
        #i+1=j
        state_out = []
        x = []
        for site_i in range(1,L):#obc
        #for site_i in range(L):#pbc
            site_j = (site_i - 1 + L) % L  # (sitei-1+L)mod L,site_j is the right side bit on site_i
            mask_up = (1 << 2*site_i) | (1 << 2*site_j)  # 1 on site i and 1 on site j
            mask_dn = (1 << (2*site_i+1)) | (1 << (2*site_j+1))
            state_out_up = (state_in ^ mask_up) # flip the two adjacent up sites
            state_out_dn = (state_in ^ mask_dn)
            x_up = t #* (1-IBITS(state_in, 2*site_i)) * (IBITS(state_in, 2*site_j))
            ##the folowing 2 lines are for fermion sign, actually not necessary in this case
            #for m in range(min([2*site_i,2*site_j])+1,max([2*site_i,2*site_j])):
            #    if(IBITS(state_in, m)): x_up *= -1
            state_out.append(state_out_up)
            x.append(x_up)
            x_dn = t #* (1-IBITS(state_in, 2*site_i+1)) * (IBITS(state_in, 2*site_j+1))
            #for m in range(min([2*site_i+1,2*site_j+1])+1,max([2*site_i+1,2*site_j+1])):
            #    if(IBITS(state_in, m)): x_dn *= -1
            state_out.append(state_out_dn)
            x.append(x_dn)
    if(term==1):
    #  term=1: -J/4*(N_i N_i+1);  #diagonal
        x = 0.
        state_out = int(state_in+0.)
        for site_i in range(1,L):#obc
        #for site_i in range(L):#pbc
            site_j = (site_i - 1 + L) % L 
            x += (IBITS(state_in, 2*site_i) + IBITS(state_in, 2*site_i+1)) * (IBITS(state_in, 2*site_j) + IBITS(state_in, 2*site_j+1)) * (-J/4.0)

    if(term==2):
    #  term=2: J*(Sz_i Sz_i+1) #diagonal Sz_i Sz_i+1 = 1/2(Nupi -Ndni)1/2(Nupi+1 -Ndni+1)=1/4(NupiNupi+1 - NupiNdni+1 -NdniNupi+1 + NdniNdni+1)
    #  term=4: lamda*(Sz_i Sz_i+delta)*(-1)^(delta+1)/abs(delta)^alpha, sum over delta=[2,L-1] for obc
        x = 0.
        state_out = int(state_in+0.)
        
        for site_i in range(L-1):#obc site_i < site_j, site_i at right side of site_j
        #for site_i in range(L):#pbc
            site_j = site_i + 1
            x += J * (0.25*(IBITS(state_in, 2*site_i)*IBITS(state_in, 2*site_j))-0.25*(IBITS(state_in, 2*site_i+1)*IBITS(state_in, 2*site_j))-0.25*(IBITS(state_in, 2*site_i)*IBITS(state_in, 2*site_j+1))+0.25*(IBITS(state_in, 2*site_i+1)*IBITS(state_in, 2*site_j+1)))
            
            for delta in range(2,L):#obc
                site_j = site_i + delta
                if(site_j<=L-1): 
                    coef_rkky = lamda *((-1)**(delta+1))/((abs(delta))**alpha)
                    x +=  ( 0.25*(IBITS(state_in, 2*site_i)*IBITS(state_in, 2*site_j))-0.25*(IBITS(state_in, 2*site_i+1)*IBITS(state_in, 2*site_j))-0.25*(IBITS(state_in, 2*site_i)*IBITS(state_in, 2*site_j+1))+0.25*(IBITS(state_in, 2*site_i+1)*IBITS(state_in, 2*site_j+1)) ) * coef_rkky
            
        
    if(term==3):      
    #  term=3: J*(S+_i S-+i+1 + hc) #off-fiagonal
    #  term=5: lamda*(S+_i S-_i+delta + hc)*(-1)^(delta+1)/abs(delta)^alpha
    #  S+_i = Cdup_i Cdn_i; S-_i+1 =Cddn_i+1 Cup_i+1
        x = []
        state_out = []
        for site_i in range(L-1):#obc site_i < site_j, site_i at right side of site_j
        #for site_i in range(L):#pbc
            site_j = site_i + 1
            if(IBITS(state_in, 2*site_j+1)*(1-IBITS(state_in, 2*site_j))*(1-IBITS(state_in, 2*site_i+1))*IBITS(state_in, 2*site_i)):
                mask_1001 = (1 << 2*site_j+1) | (1 << 2*site_j) | (1 << 2*site_i+1) | (1 << 2*site_i)  # 1 on site j with dn and up and 1 on site i with dn and up
                state_outi = (state_in ^ mask_1001)
                x.append(0.5*J)
                state_out.append(state_outi)            
            elif((1-IBITS(state_in, 2*site_j+1))*IBITS(state_in, 2*site_j)*IBITS(state_in, 2*site_i+1)*(1-IBITS(state_in, 2*site_i))):
                mask_0110 = (1 << 2*site_j+1) | (1 << 2*site_j) | (1 << 2*site_i+1) | (1 << 2*site_i)  # 1 on site j with dn and up and 1 on site i with dn and up
                state_outi = (state_in ^ mask_0110)
                x.append(0.5*J)
                state_out.append(state_outi)       
            
            for delta in range(2,L):#obc
                site_j = site_i + delta
                if(site_j<=L-1):
                    coef_rkky = lamda *((-1)**(delta+1))/((abs(delta))**alpha)
                    if(IBITS(state_in, 2*site_j+1)*(1-IBITS(state_in, 2*site_j))*(1-IBITS(state_in, 2*site_i+1))*IBITS(state_in, 2*site_i)):
                        mask_1001 = (1 << 2*site_j+1) | (1 << 2*site_j) | (1 << 2*site_i+1) | (1 << 2*site_i)  # 1 on site j with dn and up and 1 on site i with dn and up
                        state_outi = (state_in ^ mask_1001)
                        x.append(0.5*coef_rkky)
                        state_out.append(state_outi)            
                    elif((1-IBITS(state_in, 2*site_j+1))*IBITS(state_in, 2*site_j)*IBITS(state_in, 2*site_i+1)*(1-IBITS(state_in, 2*site_i))):
                        mask_0110 = (1 << 2*site_j+1) | (1 << 2*site_j) | (1 << 2*site_i+1) | (1 << 2*site_i)  # 1 on site j with dn and up and 1 on site i with dn and up
                        state_outi = (state_in ^ mask_0110)
                        x.append(0.5*coef_rkky)
                        state_out.append(state_outi)
            
    
    return state_out,x

def hami_pbc(state_in, term,t,J,lamda,alpha,L):#obc  ##term=1,2-->diagonal, tem=0,3-->off-diagonal
    # H operator act on state_in, 
    #  term=0: -t*(Cd_iC_i+1+hc); #off-diagonal
    #  term=1: -J/4*(N_i N_i+1);  #diagonal
    #  term=2: J*(Sz_i Sz_i+1) #diagonal
    #  term=3: J*(S+_i S-+i+1 + hc) #off-fiagonal
    #  term=4: lamda*(Sz_i Sz_i+delta)*(-1)^(delta+1)/abs(delta)^alpha, sum over delta=[2,L-1] for obc
    #  term=5: lamda*(S+_i S-_i+delta + hc)*(-1)^(delta+1)/abs(delta)^alpha
    if(term==0):
        # off-diagonal term, Cdup_i*Cup_i+1 - Cup_i*Cdup_i+1,Cddn_i Cdn_i+1,and hc., for all "i"
        #i+1=j
        state_out = []
        x = []
        #for site_i in range(1,L):#obc
        for site_i in range(L):#pbc
            site_j = (site_i - 1 + L) % L  # (sitei-1+L)mod L,site_j is the right side bit on site_i
            mask_up = (1 << 2*site_i) | (1 << 2*site_j)  # 1 on site i and 1 on site j
            mask_dn = (1 << (2*site_i+1)) | (1 << (2*site_j+1))
            state_out_up = (state_in ^ mask_up) # flip the two adjacent up sites
            state_out_dn = (state_in ^ mask_dn)
            x_up = t #* (1-IBITS(state_in, 2*site_i)) * (IBITS(state_in, 2*site_j))
            ##the folowing 2 lines are for fermion sign, actually not necessary in this case
            #for m in range(min([2*site_i,2*site_j])+1,max([2*site_i,2*site_j])):
            #    if(IBITS(state_in, m)): x_up *= -1
            state_out.append(state_out_up)
            x.append(x_up)
            x_dn = t #* (1-IBITS(state_in, 2*site_i+1)) * (IBITS(state_in, 2*site_j+1))
            #for m in range(min([2*site_i+1,2*site_j+1])+1,max([2*site_i+1,2*site_j+1])):
            #    if(IBITS(state_in, m)): x_dn *= -1
            state_out.append(state_out_dn)
            x.append(x_dn)
    if(term==1):
    #  term=1: -J/4*(N_i N_i+1);  #diagonal
        x = 0.
        state_out = int(state_in+0.)
        #for site_i in range(1,L):#obc
        for site_i in range(L):#pbc
            site_j = (site_i - 1 + L) % L 
            x += (IBITS(state_in, 2*site_i) + IBITS(state_in, 2*site_i+1)) * (IBITS(state_in, 2*site_j) + IBITS(state_in, 2*site_j+1)) * (-J/4.0)

    if(term==2):
    #  term=2: J*(Sz_i Sz_i+1) #diagonal Sz_i Sz_i+1 = 1/2(Nupi -Ndni)1/2(Nupi+1 -Ndni+1)=1/4(NupiNupi+1 - NupiNdni+1 -NdniNupi+1 + NdniNdni+1)
    #  term=4: lamda*(Sz_i Sz_i+delta)*(-1)^(delta+1)/abs(delta)^alpha, sum over delta=[2,L-1] for obc
        x = 0.
        state_out = int(state_in+0.)
        
        #for site_i in range(L-1):#obc site_i < site_j, site_i at right side of site_j
        for site_i in range(L):#pbc
            site_j = (site_i + 1) % L
            x += J * (0.25*(IBITS(state_in, 2*site_i)*IBITS(state_in, 2*site_j))-0.25*(IBITS(state_in, 2*site_i+1)*IBITS(state_in, 2*site_j))-0.25*(IBITS(state_in, 2*site_i)*IBITS(state_in, 2*site_j+1))+0.25*(IBITS(state_in, 2*site_i+1)*IBITS(state_in, 2*site_j+1)))
            
            #for delta in range(2,L):#obc
            #    site_j = site_i + delta
            for site_j in range(site_i+2,L):
                delta = abs(site_j - site_i) 
                if(delta>L//2-1+L%2): delta = L-delta
                if(delta == 1): 
                    #print("wrong!!!!!!!!!!!!1",site_i, site_j)
                    continue
                if(site_j<=L-1): 
                    coef_rkky = lamda *((-1)**(delta+1))/((abs(delta))**alpha)
                    x +=  ( 0.25*(IBITS(state_in, 2*site_i)*IBITS(state_in, 2*site_j))-0.25*(IBITS(state_in, 2*site_i+1)*IBITS(state_in, 2*site_j))-0.25*(IBITS(state_in, 2*site_i)*IBITS(state_in, 2*site_j+1))+0.25*(IBITS(state_in, 2*site_i+1)*IBITS(state_in, 2*site_j+1)) ) * coef_rkky
            
        
    if(term==3):      
    #  term=3: J*(S+_i S-+i+1 + hc) #off-fiagonal
    #  term=5: lamda*(S+_i S-_i+delta + hc)*(-1)^(delta+1)/abs(delta)^alpha
    #  S+_i = Cdup_i Cdn_i; S-_i+1 =Cddn_i+1 Cup_i+1
        x = []
        state_out = []
        #for site_i in range(L-1):#obc site_i < site_j, site_i at right side of site_j
        for site_i in range(L):#pbc
            site_j = (site_i + 1) % L
            if(IBITS(state_in, 2*site_j+1)*(1-IBITS(state_in, 2*site_j))*(1-IBITS(state_in, 2*site_i+1))*IBITS(state_in, 2*site_i)):
                mask_1001 = (1 << 2*site_j+1) | (1 << 2*site_j) | (1 << 2*site_i+1) | (1 << 2*site_i)  # 1 on site j with dn and up and 1 on site i with dn and up
                state_outi = (state_in ^ mask_1001)
                x.append(0.5*J)
                state_out.append(state_outi)            
            elif((1-IBITS(state_in, 2*site_j+1))*IBITS(state_in, 2*site_j)*IBITS(state_in, 2*site_i+1)*(1-IBITS(state_in, 2*site_i))):
                mask_0110 = (1 << 2*site_j+1) | (1 << 2*site_j) | (1 << 2*site_i+1) | (1 << 2*site_i)  # 1 on site j with dn and up and 1 on site i with dn and up
                state_outi = (state_in ^ mask_0110)
                x.append(0.5*J)
                state_out.append(state_outi)       
            
            #for delta in range(2,L):#obc
            #    site_j = site_i + delta
            for site_j in range(site_i+2,L):
                delta = abs(site_j - site_i)
                if(delta>L//2-1+L%2): delta = L-delta
                if(delta == 1): continue
                if(site_j<=L-1):
                    coef_rkky = lamda *((-1)**(delta+1))/((abs(delta))**alpha)
                    if(IBITS(state_in, 2*site_j+1)*(1-IBITS(state_in, 2*site_j))*(1-IBITS(state_in, 2*site_i+1))*IBITS(state_in, 2*site_i)):
                        mask_1001 = (1 << 2*site_j+1) | (1 << 2*site_j) | (1 << 2*site_i+1) | (1 << 2*site_i)  # 1 on site j with dn and up and 1 on site i with dn and up
                        state_outi = (state_in ^ mask_1001)
                        x.append(0.5*coef_rkky)
                        state_out.append(state_outi)            
                    elif((1-IBITS(state_in, 2*site_j+1))*IBITS(state_in, 2*site_j)*IBITS(state_in, 2*site_i+1)*(1-IBITS(state_in, 2*site_i))):
                        mask_0110 = (1 << 2*site_j+1) | (1 << 2*site_j) | (1 << 2*site_i+1) | (1 << 2*site_i)  # 1 on site j with dn and up and 1 on site i with dn and up
                        state_outi = (state_in ^ mask_0110)
                        x.append(0.5*coef_rkky)
                        state_out.append(state_outi)
            
    
    return state_out,x

def hami_mtx_block_pbc(k_total,L,basis_set,t,J,lamda,alpha): #k_total=0,1,2,...L-1
    basislen = len(basis_set[k_total])
    hmatrix = np.zeros(shape = (basislen,basislen),dtype=dt)
    #print("lenth of basis: ", basislen,[bin(basis_set[k_total][i]) for i in range(basislen)])
    for i in range(basislen):
        for term in [1,2]:# diagonal elements
            state_out,x0 = hami_pbc(basis_set[k_total][i], term,t,J,lamda,alpha,L)
            if(state_out!=basis_set[k_total][i]): continue
            else:
                hmatrix[i,i] += x0
        for term in [0,3]:# off-diagonal elements
            state_out,x0 = hami_pbc(basis_set[k_total][i], term,t,J,lamda,alpha,L)
            for si in range(len(x0)):
                if(state_out[si] in basis_set[k_total]):
                    out_idx = indexinlist(basis_set[k_total],state_out[si])
                    hmatrix[i,out_idx] += x0[si]
        
    values,vectors = np.linalg.eigh(hmatrix)
    #print(hmatrix)
    return values
def hami_whole_mtx_pbc(L,basis_all,t,J,lamda,alpha): #k_total=0,1,2,...L-1
    #short range terms are with coef "t" and "J", long range terms are only with coef "lamda"
    basislen = len(basis_all)
    hmatrix = np.zeros(shape = (basislen,basislen),dtype=dt)
    #print("lenth of basis: ", basislen,[bin(basis_set[k_total][i]) for i in range(basislen)])
    for i in range(basislen):
        for term in [1,2]:# diagonal elements
            state_out,x0 = hami_pbc(basis_all[i], term,t,J,lamda,alpha,L)
            if(state_out!=basis_all[i]): continue
            else:
                hmatrix[i,i] += x0
                #print(bin(basis_all[i]),x0)
        for term in [0,3]:# off-diagonal elements
            state_out,x0 = hami_pbc(basis_all[i], term,t,J,lamda,alpha,L)
            for si in range(len(x0)):
                if(state_out[si] in basis_all):
                    out_idx = indexinlist(basis_all,state_out[si])
                    hmatrix[i,out_idx] += x0[si]
        
    values,vectors = np.linalg.eigh(hmatrix)
    #print(hmatrix)
    return values,vectors


def hami_whole_mtx_obc(L,basis_all,t,J,lamda,alpha): #k_total=0,1,2,...L-1
    #short range terms are with coef "t" and "J", long range terms are only with coef "lamda"
    basislen = len(basis_all)
    hmatrix = np.zeros(shape = (basislen,basislen),dtype=dt)
    #print("lenth of basis: ", basislen,[bin(basis_set[k_total][i]) for i in range(basislen)])
    for i in range(basislen):
        for term in [1,2]:# diagonal elements
            state_out,x0 = hami_obc(basis_all[i], term,t,J,lamda,alpha,L)
            if(state_out!=basis_all[i]): continue
            else:
                hmatrix[i,i] += x0
                #print(bin(basis_all[i]),x0)
        for term in [0,3]:# off-diagonal elements
            state_out,x0 = hami_obc(basis_all[i], term,t,J,lamda,alpha,L)
            for si in range(len(x0)):
                if(state_out[si] in basis_all):
                    out_idx = indexinlist(basis_all,state_out[si])
                    hmatrix[i,out_idx] += x0[si]
    values,vectors = np.linalg.eigh(hmatrix)
    #print(hmatrix)
    return values,vectors

############### Calculate each momentum sector #########################33
def spintrans(state,L,bit): #translate the left #n bit to the right
    b = state >> L-bit
    lastbit=0
    for i in range(bit-1):
        lastbit += 2**i
    state = state << bit
    state = state & (2**L-1)
    state = state | b
    return state

def transin_heis(state, basis, L):
    ntrans = 1
    for bit in range(L):
        s = 0
        statenew = spintrans(state, L, (bit+1))
        if (statenew in basis):
            s = 1
            break
        ntrans += 1
    return (s, ntrans)

def transin(state, basis, L):
    ntrans = 1
    for bit in range(L):
        s = 0
        statenew = spintrans(state, 2*L, 2*(bit+1))
        if (statenew in basis):
            s = 1
            break
        ntrans += 1
    return (s, ntrans)

def ch(vector,k,L):#kth subblock ranges from 1 to L, corresponding to [0,2pi)
    thetak = 2 * np.pi * (k - 1) / L
    c = np.exp(1j * (vector) * thetak)
    return(c)

def find_rep(state,basis):
    idx = -1
    for i in range(len(basis)):
        if(state == basis[i]):
            idx= i
        else:
            continue
    return(idx)

def build_basis_tj(param_L, Nup,Ndn): # Using quantum number m to build bases
    L = param_L
    #maxdim = 2 ** (2*L)
    dimall,basisall = build_basis_useQN(L,Nup,Ndn)
    #print("Dimension of combenitorics: ",dimall)
    basis = []
    dim = 0
    for state in basisall:
        n_ups = 0
        #for bit in range(L):
        #    if ((state & (1 << (2*bit))) != 0):  # the bit(th) bit of state is 1
        #        n_ones += 1
        #if (n_ups != Nup):
        #    continue
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
    #print(period)
    #print(len(basis))
    print('Basis Dimention:')
    print(dim)
    return (basis,Nf,dim,period)

def total_weight(k,periodR,L):
    if(abs(( 2 * np.pi * (k - 1) / L)*periodR/np.pi/2. - int(( 2 * np.pi * (k - 1) / L)*periodR/np.pi/2.))<1e-6):
        return L/periodR
    else:
        #print(abs(k*periodR/np.pi/2. - int(k*periodR/np.pi/2.)))
        return 0.

def hami_mtx_tj_pbc(basis_in,Nf_in,dimbasis,period_in,k,t,J,lamda,alpha,L):
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
    nvector=[L,1,1,L]
    hmatrix = np.zeros(shape = (dimbasis,dimbasis),dtype=dt)
    for term in range(4): # 0:CdC // 1: NN //2: SzSz // 3: Spin flip
        for i in range(dimbasis):
            x0=Nf[i]
            coef = 0.
            for vector in range(nvector[term]):
                new_x = spintrans(basis[i],2*L,2*vector) #trasx(basis(i), vector) from 0 to L-1
                out_state, x = hami_pbc(new_x, term,t,J,lamda,alpha,L)
                #print("instate: ",new_x," outstate: ",out_state,"x: ",x)
                if (term ==1 or term ==2):
                    hmatrix[i, i] = hmatrix[i, i] + x
                if (term == 0 or term == 3):
                    for xi in range(len(x)):
                        if (np.abs(x[xi]) > 1.e-10):
                            idx = find_rep(out_state[xi], basis)
                            if (idx >= i):
                                coef = ch(vector,k,L)
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
                                    hmatrix[i,idx]=hmatrix[i,idx]+x[xi]*coef*x0/Nf[idx]/c2 #*L
                                hmatrix[idx,i]=np.conjugate(hmatrix[i,idx])

    values,vectors = np.linalg.eigh(hmatrix)
    #seed = np.zeros(hmatrix.shape[0])
    #values,vectors = lanczos(hmatrix,seed,6,1.e-5,False,False)
    return values,vectors




if __name__ == "__main__":
    L=9
    lamd = [2.0]
    J=1.0
    alpha=2.0
    S=0
    N=L+0
    Nup = 4#int((N+S)/2.) 
    Ndn = 4#int((N-S)/2.)
    t=-1.

    '''
    ##first_spin_hole_basis_obc(L)
    ##first state should be move_basis[L/2-2] for even L; move_basis[int((L-3)/2)] for odd L
    for lamda in lamd:
        #f1= open("e0_L=%d_Nup=%d_U=%.1f.txt" % (L, Nup[ns],U),"w+")

        #basis_set = [ [] for _ in range(L)]
        basis_set = []
        for l in range(L):
            basis_s = []
            for ii in range(dimq):
                #print(basisq[ii],bin(basisq[ii]))
                k=0.
                for si in range(L):
                    if(IBITS(basisq[ii], 2*si)): k+=si
                    if(IBITS(basisq[ii], 2*si+1)): k+=si
                #print("k: ",k, "      ",k%L)
                if(int(k%L)==l): basis_s.append(basisq[ii])
            basis_set.append(basis_s)
            print("length of the ",l,"th basis set: ",len(basis_s))
        #print("eigenvalues: ", hami_mtx(build_basis_useQN(L,Nup,Ndn)[1],L,all_term_list))


        gs = []
        gs0 = 0.0
        for mx in range(L):
            ek = hami_mtx_block(mx,L,basis_set,t,J,lamda,alpha)
            print("lambda=%d, blokrized k=%d: " %(lamda,mx),ek[0])#ek[0] only print minimum energy
            if(min(ek)<gs0): gs0 = min(ek)
            gs.append(ek)
            #for my in range(len(ek)):
                #f1.write("blokrized k=%d ek[%d]= %.8f\r\n" %(mx,my,ek[my]))
        gs.append(hami_mtx_block(0,L,basis_set,t,J,lamda,alpha))
        #np.savetxt("e0_L=%d_Nup=%d.txt" % (L, Nup[ns]),gs, fmt="%s")
        ks = np.linspace(0, 2 * (L) / (L), L+1)
        gs_plot = []
        #gs=np.transpose(gs)
        print(len(gs),len(gs[0]), len(ks))
        #gs_all0.append(gs0)
        #gs_all.append(gs)
        #f1.close()
    '''
    #dimq,basisq = build_basis_useQN(L,Nup,Ndn)
    #whole_eigen,whole_vec = hami_whole_mtx_pbc(L,basisq,t,J,0.0,2.)
    basis,Nf,dim,period =build_basis_tj(L, Nup,Ndn)#build_basis_heis(L, 4)#
    #print("The gs states: ",len(whole_eigen), [whole_eigen[mm] for mm in range(len(whole_eigen))])
    for k in range(L):
        values,vectors = hami_mtx_tj_pbc(basis,Nf,dim,period,k+1,t,J,0.,2.,L)
        print([values[i] for i in range(len(values))])
    print(bin(101), bin(spintrans(101,2*8,2*5)))
    
   





























        
  
