import numpy as np
import matplotlib.pyplot as plt
from hubbard_k_funcs import *


const=0
cc_terms = []
list_n = np.zeros(int(L*2))
list_nn = np.zeros(shape=(int(L*2),int(L*2)))
list_nnn = np.zeros(shape=(int(L*2),int(L*2),int(L*2)))
list_nnnn = np.zeros(shape=(int(L*2),int(L*2),int(L*2),int(L*2)))
for i in range(int(L*2)):
    list_n[i] += t*2.0*np.cos(np.pi * 2. * int(i / 2) / L)
    #const += list_n[i]
    cc_terms.append(c_terms([i+1, -i-1],list_n[i]))
    print(" cc terms ",[i+1, -i-1],list_n[i])

k = np.zeros(L)
p = np.zeros(L)
q = np.zeros(L)

for i in range(L): #momentum
    k[i] = 2 * np.pi * i / L
    p[i] = 2 * np.pi * i / L
    q[i] = 2 * np.pi * i / L
#print("k",k)
cccc_num=0
nn_num=0
goal = [1,2,3,4]
cccc_terms=[]
nn_terms=[]
for k1 in range(L):
    for k2 in range(L):
        for k3 in range(L):
            k4=(k1-k2+k3+100*L)%L
            #print(2 * k1, 2 * k2, 2 * k3 + 1, 2 * k4 + 1)
            if(k1==k2):
                nn_terms.append([2 * k1+1, -(2 * k2+1), 2 * k3 + 2, -(2 * k4 + 2)])
                nn_num+=1
                list_nn[2 * k1][2 * k3 + 1] += U
                cccc_terms.append(c_terms([2 * k1+1, -(2 * k2+1), 2 * k3 + 2, -(2 * k4 + 2)],U/L))

            else:
                cccc_terms.append(c_terms([2 * k1+1, -(2 * k2+1), 2 * k3 + 2, -(2 * k4 + 2)],U/L))
                cccc_num+=1
#print(nn_terms)
for i in range(len(cccc_terms)):
    print(i,cccc_terms[i].lst,cccc_terms[i].coef)  ##odd is spin up on site (num-1)/2, even is spin down on site (num-2)/2
#print(nn_num,cccc_num)

all_term_list = []
for i in cc_terms:
    all_term_list.append(i)

for i in cccc_terms:
    all_term_list.append(i)


S = [2,4]

Nup = [int((N+S[i])/2.) for i in range(len(S))]
Ndn = [int((N-S[i])/2.) for i in range(len(S))]
print("Nup=",Nup," Ndn=",Ndn)
gs_all = []
gs_all0 = []
for ns in range(len(S)):
    f1= open("e0_L=%d_Nup=%d_U=%.1f.txt" % (L, Nup[ns],U),"w+")
    dimq,basisq = build_basis_useQN(L,Nup[ns],Ndn[ns])
    basis_set = [ [] for _ in range(L)]
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
   
    print(len(basis_set[0]))
    #print("eigenvalues: ", hami_mtx(build_basis_useQN(L,Nup,Ndn)[1],L,all_term_list))


    gs = []
    gs0 = 0.0
    for mx in range(L):
        ek = hami_mtx_block(mx,L,all_term_list,basis_set)
        print("blokrized k=%d: " %(mx),ek)
        if(min(ek)<gs0): gs0 = min(ek)
        gs.append(ek)
        for my in range(len(ek)):
            f1.write("blokrized k=%d ek[%d]= %.8f\r\n" %(mx,my,ek[my]))
    gs.append(hami_mtx_block(0,L,all_term_list,basis_set))
    #np.savetxt("e0_L=%d_Nup=%d.txt" % (L, Nup[ns]),gs, fmt="%s")
    ks = np.linspace(0, 2 * (L) / (L), L+1)
    gs_plot = []
    #gs=np.transpose(gs)
    print(len(gs),len(gs[0]), len(ks))
    gs_all0.append(gs0)
    gs_all.append(gs)
    f1.close()
    '''
    fig, ax1 = plt.subplots()
    for i in range(7):#(len(gs[0])):
        gs_plot = [gs[ii][i] for ii in range(len(gs))]#-gs0
        ax1.scatter(ks, gs_plot)#,s=35)
        ax1.plot(ks, gs_plot)
    ax1.set_title(r'Spectrum of L=%d, N=%d, S=%d' %(L,N,S[ns]),fontsize=16)#, x=0.5, y=0.86)

    plt.show()
    '''
ks = np.linspace(0, 2 * (L) / (L), L+1)
fig, ax1 = plt.subplots()

lengs1x = len(gs_all[1])
lengs1y = min([len(gs_all[1][i]) for i in range(L+1)])
lengs0x = len(gs_all[0])
lengs0y = min([len(gs_all[0][i]) for i in range(L+1)])

newgs=np.zeros(shape=(lengs0x,lengs0y))
for i in range(lengs0x):
    for j in range(lengs0y):
        newgs[i][j] = gs_all[0][i][j]



del_list=[]
for nx in range(lengs1x):
    for ny in range(lengs1y):
        for mx in range(lengs0x):
            for my in range(lengs0y):
                #print(nx,ny,mx,my)
                if(abs(gs_all[0][mx][my]-gs_all[1][nx][ny])< 1e-10):
                    if([mx,my] not in del_list): 
                        del_list.append([mx,my])
                        print("delete: " , [mx,my])
                    


for ay in range(lengs0y):
    ks_plot=[]
    plot_line=[]
    for ax in range(lengs0x):
        if([ax,ay] in del_list): continue
        else: 
            ks_plot.append(ks[ax])
            plot_line.append(gs_all[0][ax][ay]-(-10.6144071606))#-gs_all0[0])#-(-7.952325596992009))
    if(ay<100):
        ax1.scatter(ks_plot, plot_line)#,s=35)
        ax1.plot(ks_plot,plot_line)

ax1.set_title(r'Spectrum of L=%d, N=%d, S=%d - S=%d' %(L,N,int(S[0]/2),int(S[1]/2)),fontsize=16)#, x=0.5, y=0.86)
ax1.set_xlabel('$k/\pi$',fontsize=14)
ax1.set_ylabel('E',fontsize=14)
ax1.set_ylim([0,7])
plt.show()
print("gs:",gs_all0)

