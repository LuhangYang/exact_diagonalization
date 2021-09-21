import numpy as np
import matplotlib.pyplot as plt
from tj_ED_funcs import *


L=9
J=1.0
lamda = 0.0+J
t=-1.
dt = np.dtype(np.complex128)
alpha_all=[1.,2.0,3.,4.,20.]
hole = [0]
Nup = 4
Ndn = 4
print("Nup=",Nup," Ndn=",Ndn)
gs_all = []
gs_all0 = []
colors=['blue','red','gold','brown','green']
fig, ax1 = plt.subplots()
ic=0
basis,Nf,dim,period =build_basis_tj(L, Nup,Ndn)
print("length of basis: ",dim)
for alpha in alpha_all:
    for ns in range(1):
        f1= open("e0_new_L=%d_Nup=%d_Ndn=%d_J=%.1f_alpha=%.1f.txt" % (L, Nup,Ndn,J,alpha),"w+")
        #dimq,basisq = build_basis_useQN(L,Nup[ns],Ndn[ns])
        gs = []
        gs0 = 0.0
        for kx in range(L+1):
            #ek = hami_mtx_block_pbc(mx,L,basis_set,t,J,lamda,alpha)
            ek,vectos = hami_mtx_tj_pbc(basis,Nf,dim,period,kx+1,t,J,lamda,alpha,L)
            print("blokrized k=%d: " %(kx+1),ek)
            if(min(ek)<gs0): gs0 = min(ek)
            gs.append(ek)
            for my in range(len(ek)):
                f1.write("blokrized k=%d ek[%d]= %.8f\r\n" %(kx,my,ek[my]))
        #gs.append(hami_mtx_tj_pbc(basis,Nf,dim,period,1-L//2,t,J,lamda,alpha,L)[0])
        #np.savetxt("e0_L=%d_Nup=%d.txt" % (L, Nup[ns]),gs, fmt="%s")
        ks = np.linspace(0, 2 * (L) / (L), L+1)
        gs_plot = []
        #gs=np.transpose(gs)
        print(len(gs),len(gs[0]), len(ks))
        gs_all0.append(gs0)
        gs_all.append(gs)
        f1.close()
        
        for i in range(1):#(len(gs[0])):#(1):#:
            gs_plot = [gs[ii][i]-gs[0][0] for ii in range(len(gs))]#-gs0
            print(ks,gs_plot)
            ax1.scatter(ks, gs_plot,color=colors[ic])#,s=35)
            ax1.plot(ks, gs_plot,color=colors[ic],label=r'$\alpha=%.1f$' %(alpha))
        ic+=1
ax1.set_title(r'Spectrum of L=%d, Nup=%d, Ndn=%d, J=%.1f '%(L,Nup,Ndn,J),fontsize=16)#, x=0.5, y=0.86)

plt.legend()
    #plt.savefig(r'Spectrum_all_L=%d_Nup=%d_Ndn=%d_J=%.1f_alpha=%.1f.png' %(L,Nup[ns],Ndn[ns],J,alpha),fontsize=16)#, x=0.5, y=0.86)
plt.savefig(r'gs_energy_L=%d_Nup=%d_Ndn=%d_J=%.1f_alpha=%.1f.png' %(L,Nup,Ndn,J,alpha),fontsize=16)#, x=0.5, y=0.86)
plt.show()
    #plt.clf()







