import numpy as np
import matplotlib.pyplot as plt

L = 10
N = 10
#S = 6
dt = np.dtype(np.complex128)

t = -1.
U = 4.


class c_terms(object):
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


def IBITS(n,i): #count from right site from 0.
    return ((n >> i) & 1)


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
    return dim1,basis1


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
    ntrans = 1
    for bit in range(L):
        s = 0
        statenew = spintrans(state, 2*L, 2*(bit+1))
        if (statenew in basis):
            s = 1
            break
        ntrans += 1
    return (s, ntrans)

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
    

def hami(state_in, term): # H operator act on state_in, term =c_term([],coef)
    x=0.
    if(len(term.lst) == 2):#ck+ck term
        # diagonal term
        state_out = state_in
        site_k = term.lst[0]-1
        x = term.coef*IBITS(state_in, site_k)
        
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
            
    return state_out,x
 

def hami_mtx(basis,L,all_term_list):
    hmatrix = np.zeros(shape = (len(basis),len(basis)),dtype=dt)
    basislen = len(basis)
    print("lenth of basis: ", basislen)
    termlen = len(all_term_list)
    for term in all_term_list: 
        for i in range(basislen):
            state_out,x0 = hami(basis[i], term)
            if(state_out in basis):
                #print(state_out,indexinlist(basis,state_out))
                out_idx = indexinlist(basis,state_out)
                hmatrix[i,out_idx] += x0
    values,vectors = np.linalg.eigh(hmatrix)
    #print(hmatrix)
    return values


def hami_mtx_block(k_total,L,all_term_list,basis_set): #k_total=0,1,2,...L-1
    basislen = len(basis_set[k_total])
    hmatrix = np.zeros(shape = (basislen,basislen),dtype=dt)
    #print("lenth of basis: ", basislen,[bin(basis_set[k_total][i]) for i in range(basislen)])
    termlen = len(all_term_list)
    for term in all_term_list: 
        for i in range(basislen):
            state_out,x0 = hami(basis_set[k_total][i], term)
            if(state_out in basis_set[k_total]):
                #if(abs(x0)>1e-8): print(basis_set[k_total][i], term.lst,state_out,indexinlist(basis_set[k_total],state_out),x0)
                out_idx = indexinlist(basis_set[k_total],state_out)
                hmatrix[i,out_idx] += x0
    values,vectors = np.linalg.eigh(hmatrix)
    #print(hmatrix)
    return values


def hami_real_space(state_in, term,t,U): # H operator act on state_in,
    x = 0.
    if(term == 0):
        # diagonal term, U*niup*nidn
        state_out = state_in
        for site_i in range(L): 
            up_site = 2*site_i  #spin up
            dn_site = up_site+1  # spin dn
            x += IBITS(state_in, up_site) * IBITS(state_in, dn_site) *U
    else:
        # off-diagonal term, Cdup_i Cup_j
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
            #if(IBITS(state_in, (2*site_i-1+2*L)%(2*L))):
            #    x_up *= -1
            state_out.append(state_out_up)
            x.append(x_up)
            x_dn = (1-IBITS(state_in, 2*site_i+1)) * (IBITS(state_in, 2*site_j+1)) * t
            state_out.append(state_out_dn)
            for m in range(min([2*site_i+1,2*site_j+1])+1,max([2*site_i+1,2*site_j+1])):
                if(IBITS(state_in, m)): x_dn *= -1
            #if(IBITS(state_in, 2*site_i)):
            #    x_dn *= -1
            x.append(x_dn)
        #x = ... # return a value of you could flip them, otherwise 0
    return (state_out,x) # after the operator's action, becomes stateout with coefficient x


def ch(vector,k):#kth subblock ranges from 1 to L, corresponding to [0,2pi)
    thetak = 2 * np.pi * (k - 1) / L
    #print("thetak: ",thetak)
    c = np.exp(1j * (vector) * thetak)
    #print("c: ",c)
    return(c)

def hami_mtx_k(k,L,Nup,Ndn,t,U):#wrong, don't use
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
                out_state, x = hami_real_space(new_x, term,t,U)
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


if __name__ == "__main__":
    #print("check fermion sign: ", fermion_sign([1, -3, 2, -12],228))#fermion_sign([3, -7, 2, -6]))
    #print("k= 0, ",hami_mtx_k(0,L,Nup,Ndn,t,U))
    #print("k= 1, ",hami_mtx_k(1,L,Nup,Ndn,t,U))
    #print("k= 2, ",hami_mtx_k(2,L,Nup,Ndn,t,U))
    #print("k= 3, ",hami_mtx_k(3,L,Nup,Ndn,t,U))





