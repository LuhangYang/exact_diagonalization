import numpy as np
import matplotlib.pyplot as plt


##################### Global variables ########################################
L = 4
#N = 1
#S = 6
dt = np.dtype(np.complex128)
#lamda = 1.0

t = -1.
U = -1.0
V=0.0

##################### Commonly used  functions ###############################

def IBITS(n,i): #count from right site from 0.
    return ((n >> i) & 1)
#print(bin(19),IBITS(19,0),IBITS(19,1),IBITS(19,2),IBITS(19,3),IBITS(19,4),IBITS(19,5),IBITS(19,6))

def IthBITinN(n,i):
    return ((n >> i) & 1)

#translate the spin chain
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

def transin(state, basis, L):
    #    print("state: ",state)
    ntrans = 1
    for bit in range(L):
        s = 0
        statenew = spintrans(state, 2*L, 2*(bit+1))
        #        print("statenew: ", statenew,"L: ",L,"bit:",bit)
        if (statenew in basis):
            s = 1
            break
        ntrans += 1
        #    print("state: ",state,"basis: ",basis,"return value: ",s)
    return (s, ntrans)


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
#list = [0,1,2,3,4,5]
#print("check swap: ",swapPositions(list,3,0))

def ch(vector,k):#kth subblock ranges from 1 to L, corresponding to [0,2pi)
    thetak = 2 * np.pi * (k - 1) / L
    #print("thetak: ",thetak)
    c = np.exp(1j * (vector) * thetak)
    #print("c: ",c)
    return(c)

##################### Bulid Basis ########################################
def build_basis_useQN(L1,Nup1,Ndn1): #count site from right to left
    maxdim = 4 ** L1
    basis1 = []
    dim1=0
    for state in range(maxdim):
        nup_ones = 0
        ndn_ones = 0
        for bit in range(L1):
            nup_ones += IBITS(state,2*bit) #(dn up) (dn up) ...
            ndn_ones += IBITS(state,2*bit+1)
        #print (state,nup_ones,ndn_ones)
        if(nup_ones == Nup1 and ndn_ones == Ndn1):
            basis1.append(state)
            dim1 += 1
    print ("Dim=",dim1)
    print ("Basis:")
    print (basis1)
    print ([bin(basis1[i]) for i in range(dim1)])
    return dim1,basis1

def build_eta_basis(basis_k,L1):
    basis_et = []
    for state_basis in basis_k:
        is_eta = 1
        for bit in range(L1):
            if(IBITS(state_basis,2*bit) != IBITS(state_basis,2*bit+1)): 
                is_eta = 0
                continue
        if(is_eta==1): basis_et.append(state_basis)
    return basis_et
        
            
    



def build_basis_k(L1, Nup1,Ndn1): # Using quantum number m to build bases
    maxdim = 4 ** L1
    basis1 = []
    dim1=0
    for state in range(int(maxdim / 2) + 1):
        nup_ones = 0
        ndn_ones = 0
        for bit in range(L1):
            nup_ones += IBITS(state,2*bit)
            ndn_ones += IBITS(state,2*bit+1)
        if(nup_ones != Nup1 or ndn_ones != Ndn1):
            continue
        s, ntrans = transin(state, basis1, L1)
        if (s):
            basis1 = basis1
        else:
            basis1.append(state)
            dim1 += 1
    Nf1 = []  # normalization factor for each basis
    period1 = []
    for i in range(dim1):
        s, ntrans = transin(basis1[i], basis1, L1)
        if (ntrans < L1):
            normal = 1 / np.sqrt((L1 / (ntrans)) ** 2 * (ntrans))
            period1.append(ntrans)
        else:
            normal = 1 / np.sqrt(L1)
            period1.append(ntrans)
        Nf1.append(normal)
    
    return (basis1,Nf1,dim1,period1)



##################### CC CCCC terms Class ########################################


class c_terms(object):##odd is spin up on site (num-1)/2, even is spin down on site (num-2)/2
    def __init__(self, lst=[], coef=1.0):
        self.lst = lst
        self.coef = coef

    def __mul__(self, other): #self*other
        if(type(other)==c_terms):
            new_lst = []
            for i in range(len(self.lst)):
                new_lst.append(self.lst[i])
            for i in range(len(other.lst)):
                new_lst.append(other.lst[i])
            return c_terms(new_lst,self.coef*other.coef)
        elif(type(other)==float):
            return c_terms(self.lst, self.coef * other)
    def __rmul__(self, other): #other*self
        if(type(other)==c_terms):
            new_lst = []
            for i in range(len(other.lst)):
                new_lst.append(other.lst[i])
            for i in range(len(self.lst)):
                new_lst.append(self.lst[i])
            return c_terms(new_lst,self.coef*other.coef)
        elif (type(other) == float):
            return c_terms(self.lst, self.coef * other)




##################### Fermion sign ########################################


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
    
def fermion_sign_cdc(lst,state):
    lstnew = [lst[0],-lst[1]]
    lstnew.sort(reverse=True)  #descenting order
    #print("new list: ",lstnew)
    indexnew = [indexinlist(lst,lstnew[0]),indexinlist(lst,lstnew[1])]
    #print(indexnew)
    sign=1.0
    for i in range(2):
        sign *= (-1)**(indexinlist(indexnew,i)-i)
        indexnew = swapPositions(indexnew,indexinlist(indexnew,i),i)
        #print("i: ",i, sign, indexnew)
    lstn = [lstnew[xi]-1 for xi in range(2)]
    for ii in range(2*L):
        if(IBITS(state,2*L-ii)): lstn.append(2*L-ii)
        #print(lstn)
    for i in ([1,0]):
        for j in range(2,len(lstn)):
            if(lstn[i]< lstn[j]): 
                sign *= -1               
    return sign

##################### Bulid Terms ########################################
def creat_v_terms(ud1,ud2,ud3,ud4,V):#dagger notdagger dagger notdagger, ud1=1--down, 0--up
    cccc_termss=[] 
    for k1 in range(L):
        for k2 in range(L):
            for k3 in range(L):
                k4=(k1-k2+k3+100*L)%L
                #print(2 * k1, 2 * k2, 2 * k3 + 1, 2 * k4 + 1)   
                cccc_termss.append(c_terms([2 * k1+1, -(2 * k2+1), 2 * k3 + 2, -(2 * k4 + 2)],np.exp(2.*np.pi*(k4-k3)*1.j/L)*V/L))
    return cccc_termss
                


#print("check fermion sign: ", fermion_sign([1, -3, 2, -12],228))#fermion_sign([3, -7, 2, -6]))
##################### H operator ########################################
def hami(state_in, term): # H operator act on state_in, term =c_term([],coef)
    x=0.
    if(len(term.lst) == 2):
        # diagonal term #ck+ck term or nk term
        [s1,s2] = term.lst
        if(s1==-s2):
            state_out = state_in
            site_k = term.lst[0]-1
            x = term.coef*IBITS(state_in, site_k)
        #off-diag cdicj term
        else:
            mask1 = (1 << (-s2-1)) | (1 << (s1-1))  # 1 on site i and 1 on site j
            state_out = (state_in ^ mask1)  # flip the two 
            x = (1-IBITS(state_in, s1-1))*(IBITS(state_in, -s2-1))*term.coef*fermion_sign_cdc(term.lst,state_in)
        
    elif(len(term.lst) == 4):
        # off-diagonal term
        [s1,s2,s3,s4] = term.lst
        if(s1==-s2 and s3==-s4):
            state_out = state_in
            x = term.coef*IBITS(state_in, s1-1)*IBITS(state_in, s3-1)
        else:
            mask1 = (1 << (-s2-1)) | (1 << (s1-1))  # 1 on site i and 1 on site j
            mask2 = (1 << (-s4-1)) | (1 << (s3-1))
            state_outi = (state_in ^ mask1)  # flip the two 
            state_out = state_outi ^ mask2
            x = (1-IBITS(state_in, s1-1))*(IBITS(state_in, -s2-1))*(1-IBITS(state_in, s3-1))*(IBITS(state_in, -s4-1))*term.coef*fermion_sign(term.lst,state_in)
            #print(state_in, [s1,s2,s3,s4],(1-IBITS(state_in, s1-1)),(IBITS(state_in, -s2-1)),(1-IBITS(state_in, s3-1)),(IBITS(state_in, -s4-1)) , fermion_sign(term.lst,state_in))
    return state_out,x
 
#print(bin(228))
'''
print(hami(228, cc_terms[0]))
print(hami(228, cc_terms[1]))
print(hami(228, cc_terms[2]))
print(hami(228, cc_terms[3]))
print(hami(228, cc_terms[4]))
print(hami(228, cc_terms[5]))
print(hami(228, cc_terms[6]))
print(hami(228, cc_terms[7]))
'''
#
##################### Calculate Hamiltonian ################################
def hami_mtx(basis,L,all_term_list):
    hmatrix = np.zeros(shape = (len(basis),len(basis)),dtype=dt)
    basislen = len(basis)
    print("lenth of basis: ", basislen)
    termlen = len(all_term_list)
    for term in all_term_list: 
        for i in range(basislen):
            state_out,x0 = hami(basis[i], term)
            if(state_out in basis):
                out_idx = indexinlist(basis,state_out)
                if(out_idx==i):
                    hmatrix[i,out_idx] += x0
                elif(out_idx>i):
                    hmatrix[i,out_idx] += x0
                    hmatrix[out_idx,i] = np.conjugate(hmatrix[i,out_idx])
    
    return hmatrix#values,np.transpose(vectors)[0]

def hami_mtx_lri(basis,L,all_term_list):
    hmatrix = np.zeros(shape = (len(basis),len(basis)),dtype=dt)
    hszsz = np.zeros(shape = (len(basis),len(basis)),dtype=dt)
    basislen = len(basis)
    print("lenth of basis: ", basislen)
    termlen = len(all_term_list)
    for term in all_term_list: 
        for i in range(basislen):
            state_out,x0 = hami(basis[i], term)
            if(state_out in basis):
                #print(state_out,indexinlist(basis,state_out))
                out_idx = indexinlist(basis,state_out)
                if(out_idx==i):
                    hmatrix[i,out_idx] += x0
                elif(out_idx>i):
                    hmatrix[i,out_idx] += x0
                    hmatrix[out_idx,i] = np.conjugate(hmatrix[i,out_idx])
    for iba in range(basislen):
        for i in range(L-1):
            for j in range(i+1,L):
                coef = ((-1)**abs(j-i+1))*lamda/(j-i)**2
                hmatrix[iba,iba] += coef*0.25*(IBITS(basis[iba], 2*i)-IBITS(basis[iba], 2*i+1))*(IBITS(basis[iba], 2*j)-IBITS(basis[iba], 2*j+1))
                hszsz[iba,iba] += coef*0.25*(IBITS(basis[iba], 2*i)-IBITS(basis[iba], 2*i+1))*(IBITS(basis[iba], 2*j)-IBITS(basis[iba], 2*j+1))
    print("hszsz:   ")
    for i in range(basislen):
        for j in range(basislen):
            if(hszsz[i,j]*np.conjugate(hszsz[i,j]) > 0.00001):
                print(i,j,hszsz[i,j])
    #values,vectors = np.linalg.eigh(hmatrix)
    #print(hmatrix)
    return hmatrix#values,np.transpose(vectors)[0]

def hami_mtx_block(k_total,L,all_term_list,basis_set): #k_total=0,1,2,...L-1
    basislen = len(basis_set[k_total])
    print("length of the basis_ste:", len(basis_set))
    hmatrix = np.zeros(shape = (basislen,basislen),dtype=dt)
    #print("lenth of basis: ", basislen,[bin(basis_set[k_total][i]) for i in range(basislen)])
    termlen = len(all_term_list)
    for term in all_term_list: 
        for i in range(basislen):
            state_out,x0 = hami(basis_set[k_total][i], term)
            if(state_out in basis_set[k_total]):
                #if(abs(x0)>1e-8): print(basis_set[k_total][i], term.lst,state_out,indexinlist(basis_set[k_total],state_out),x0)
                out_idx = indexinlist(basis_set[k_total],state_out)
                if(out_idx>=i):
                    hmatrix[i,out_idx] += x0
                    hmatrix[out_idx,i] = np.conjugate(hmatrix[i,out_idx])
    values,vectors = np.linalg.eigh(hmatrix)
    #print(hmatrix)
    return values


def hami_k(state_in, term,t,U): # H operator act on state_in,
    x = 0.
    if(term == 0):
        # diagonal term
        state_out = state_in
        for site_i in range(L): 
            up_site = 2*site_i  #spin up
            dn_site = up_site+1  # spin dn
            x += IBITS(state_in, up_site) * IBITS(state_in, dn_site) *U
    else:
        # off-diagonal term
        state_out = []
        x = []
        for site_i in range(L):
            site_j = (site_i - 1 + L) % L  # (sitei-1+L)mod L,site_j is the right side bit on site_i

            mask_up = (1 << 2*site_i) | (1 << 2*site_j)  # 1 on site i and 1 on site j
            mask_dn = (1 << (2*site_i+1)) | (1 << (2*site_j+1))
            state_out_up = (state_in ^ mask_up) # flip the two adjacent
            state_out_dn = (state_in ^ mask_dn)
            x_up = (1-IBITS(state_in, 2*site_i)) * (IBITS(state_in, 2*site_j)) * t
            for m in range(min([2*site_i,2*site_j])+1,max([2*site_i,2*site_j])):
                if(IBITS(state_in, m)): x_up *= -1
            
            state_out.append(state_out_up)
            x.append(x_up)
            x_dn = (1-IBITS(state_in, 2*site_i+1)) * (IBITS(state_in, 2*site_j+1)) * t
            state_out.append(state_out_dn)
            for m in range(min([2*site_i+1,2*site_j+1])+1,max([2*site_i+1,2*site_j+1])):
                if(IBITS(state_in, m)): x_dn *= -1
            
            x.append(x_dn)

    return (state_out,x) 

def hami_mtx_k(k,L,Nup,Ndn,t,U):
    nvector=[1,L]
    basis,Nf,dimbasis,period = build_basis_k(L,Nup,Ndn)
    print([bin(basis[i]) for i in range(dimbasis)])
    print([period[i] for i in range(dimbasis)])
    print([Nf[i] for i in range(dimbasis)])
    hmatrix = np.zeros(shape = (len(basis),len(basis)),dtype=dt)
    for term in range(2): # 0 is U*nupndn , 1 is -t*(CdupiCupi+1 CddniCdni+1)
        for i in range(dimbasis):
            x0=Nf[i]
            coef = 0.
            for vector in range(nvector[term]):
                new_x = spintrans(basis[i],2*L,2*vector) #trasx(basis(i), vector) from 0 to L-1
                out_state, x = hami_k(new_x, term,t,U)
                if(term==1):
                    print("instate: ",new_x,bin(new_x)," outstate: ",[bin(out_state[xxx]) for xxx in range(len(out_state)) if (out_state[xxx] in basis)],"x: ",[x[xxx] for xxx in range(len(out_state)) if (out_state[xxx] in basis)])
                if (term == 0):#out_state=statein
                    hmatrix[i, i] = hmatrix[i, i] + x
                if (term == 1):
                    for xi in range(len(x)):
                        if (np.abs(x[xi]) > 1.e-10 and (out_state[xi] in [spintrans(basis[ii],2*L,2*vector) for ii in range(dimbasis)])):
                            idx = indexinlist([spintrans(basis[ii],2*L,2*vector) for ii in range(dimbasis)],out_state[xi]) 
                            #print("idx: ",idx)
                            if (idx >= i):
                                #print("vector:",vector)
                                #print(Nf[idx]-1/np.sqrt(L))
                                coef = ch(vector,k)
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
                            #print("xi: ",x[i],"x0: ",x0,"Nf[idx]: ",Nf[idx])
    #print(hmatrix)
    values,vectors = np.linalg.eigh(hmatrix)
    return values


##################### Create mixed basis ########################################


def mix_basis(hubbard_basis, heis_basis,L):#the left 2L bits are hubb, the right L buts are heis
    basis_al = []
    for i in hubbard_basis:
        for j in heis_basis:
            basis_al.append(i*(2**(L)) + j)
    return basis_al


#################### Check whether eta pairing states #####################3

def eta_or_not(statei,L):
    is_or_not=1 #0 not, 1 is
    for ibi in range(L):
        if(IBITS(statei,2*ibi) != IBITS(statei,2*ibi+1)):
            is_or_not=0
    return is_or_not











