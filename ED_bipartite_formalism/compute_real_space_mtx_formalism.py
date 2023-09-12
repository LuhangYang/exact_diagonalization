import numpy as np
from hubbard_k_funcs_mtx_formalism import *


######################## build basis for each k sector

dim, basisup = build_all_basis(Lx*Ly,int(Lx*Ly/2)) 
basisdn = np.copy(basisup)
print(dim)
dt = float

###############3# calculate eigen energy, eigen vectors, and H matrix
S = BipartiteSystem(L,basisup,dt,U)

seed = np.eye(dim)
print(seed)
e0,gs = lanczos_customize(S, seed, 99, 1e-8, use_seed = False, force_maxiter = False) 
print("energys :", e0)
print("psi: ", gs)


################################ Compute correction vector in r-space ##############################

#gs = eigenmatrix_hami.transpose()[0]*1.0


#-1.9531014996351634 #t=-1, U=4
#-4.472135954999584  #t=-1, U=0



