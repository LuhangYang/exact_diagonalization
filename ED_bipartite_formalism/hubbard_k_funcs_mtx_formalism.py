import numpy as np
import matplotlib.pyplot as plt

Lx = 4
Ly = 1
L = Lx*Ly
#N = 8
#S = 6
dt = np.dtype(np.complex128)

t = -1.
U = 4.


class BoundaryCondition:
      OBC, PBC = range(2)

xBC = BoundaryCondition.OBC
#yBC = BoundaryCondition.PBC



def psi_dot_psi(psi1, psi2): #<psi1|psi2>
    x = 0.
    for i in range(psi1.shape[0]):
        for j in range(psi2.shape[1]):
            x += np.conjugate(psi1[i,j])*psi2[i,j]
    print(x)
    return x

def lanczos_customize(m, seed, maxiter, tol, use_seed = False, force_maxiter = False):
    x1 = seed
    x2 = seed
    gs = seed
    a = np.zeros(100,dtype=m.dt)
    b = np.zeros(100,dtype=m.dt)
    z = np.zeros((100,100),dtype=m.dt)
    lvectors = []
    control_max = maxiter
    e0 = 9999

    if(maxiter == -1):
        force_maxiter = False

    if(control_max == 0):
        gs = 1
        maxiter = 1
        return(e0,gs)
    
    x1[:,:] = 0
    x2[:,:] = 0
    gs[:,:] = 0
    a[:] = 0.0
    b[:] = 0.0
    if(use_seed):
        x1 = seed
    else:
        for i in range(x1.shape[0]):
            for j in range(x1.shape[1]):
                x1[i,j] = (2*np.random.random()-1.)

#    x1[:,:] = 1
    b[0] = psi_dot_psi(x1,x1)
    b[0] = np.sqrt(b[0])
    x1 = x1 / b[0]
    x2[:] = 0
    b[0] = 1.

    e0 = 9999
    nmax = min(99, maxiter)
    print("nmax: ", nmax)
    for iter in range(1,nmax+1):
        eini = e0
        if(b[iter - 1] != 0.):
            aux = x1
            x1 = -b[iter-1] * x2
            x2 = aux / b[iter-1]
        #print("state: ",x2)
        aux = m.product(x2)

        x1 = x1 + aux
        a[iter] = psi_dot_psi(x1,x2)
        x1 = x1 - x2*a[iter]

        b[iter] = psi_dot_psi(x1,x1)
        b[iter] = np.sqrt(b[iter])
        lvectors.append(x2)                                                  
#        print ("Iter =",iter,a[iter],b[iter])
        z.resize((iter,iter))
        z[:,:] = 0
        for i in range(0,iter-1):
            z[i,i+1] = b[i+1]
            z[i+1,i] = b[i+1]
            z[i,i] = a[i+1]
        z[iter-1,iter-1]=a[iter]
        d, v = np.linalg.eig(z)

        col = 0
        n = 0
        e0 = 9999
        for e in d:
            if(e < e0):
                e0 = e
                col = n
            n+=1
        e0 = d[col]
        
       
        print ("Iter = ",iter," Ener = ",e0) #check the convergence of the gs
        if((force_maxiter and iter >= control_max) or (iter >= gs.shape[0]*gs.shape[1] or iter == 99 or abs(b[iter]) < tol) or \
            ((not force_maxiter) and abs(eini-e0) <= tol)):
            # converged
            gs[:,:] = 0.
            for n in range(0,iter):
                gs += v[n,col]*lvectors[n]

            print ("E0 = ", e0)
            maxiter = iter
            return(e0,gs) # We return with ground states energy

    return(e0,gs)



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


def build_all_basis(param_L, Nup): # Using quantum number m to build bases
    L = param_L
    maxdim = 2 ** L
    basis = []
    dim = 0
    for state in range(int(maxdim)):
        n_ones = 0
        for bit in range(L):
            if ((state & (1 << bit)) != 0):  # the bit(th) bit of state is 1
                n_ones += 1
        if (n_ones != Nup):
            continue
        else:
            basis.append(state)
            dim += 1
    print("dim: ",dim)
    return dim,basis


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








