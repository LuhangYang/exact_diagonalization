from spinless_fermion_ED_funcs import *
import numpy as np
import matplotlib.pyplot as plt

V4=8.#-500.0
mu=0.#1000./L
def measureSz(statein,k,r):
	sk = 0.
	for i in range(L):
		if(IBITS(statein,i)):
			sk += 0.5*np.exp(-1j*k*((i+10*L+r)%L))
		else:
			sk += -0.5*np.exp(-1j*k*((i+10*L+r)%L))
	return sk
def measureNiNj(basis,vect,L):
	nn = []
	for i in range(L):
		nin = []
		for j in range(L):
			ninj = 0.
			for bi in range(len(basis)):
				if(IBITS(basis[bi],j)==1 and IBITS(basis[bi],i)==1):
					ninj += (vect[bi]*np.conjugate(vect[bi])).real
			nin.append(ninj)
		nn.append(nin)
	return nn

def measureCdiCj(basis,vect,L):
	cdc = np.zeros(shape=(L, L),dtype=dt)
	for bi in range(len(basis)):
		for i in range(L):
			if(IBITS(basis[bi],i)==1):
				cdc[i][i] += vect[bi]*np.conjugate(vect[bi])
			for j in range(i+1,L):
				num = 0
				for num_i in range(i+1,L):
					if(IBITS(basis[bi],num_i)==1):num+=1
				for num_i in range(j+1,L):
					if(IBITS(basis[bi],num_i)==1):num+=1
				if(IBITS(basis[bi],i)==1 and IBITS(basis[bi],j)==0):#i<j CdjCi

					mask = (1 << i) | (1 << j)
					new_state = (basis[bi] ^ mask)
					bj = find_rep(new_state, basis)
					#if(i==0):
					#	print(i,j,bi,bj,num)
					cdc[i][j] += ((-1)**num) *vect[bi]*np.conjugate(vect[bj])
				if(IBITS(basis[bi],i)==0 and IBITS(basis[bi],j)==1):#i<j CdiCj
					mask = (1 << i) | (1 << j)
					new_state = (basis[bi] ^ mask)
					bj = find_rep(new_state, basis)
					#if(i==0):
					#	print(j,i,bi,bj,num+1)
					cdc[j][i] += ((-1)**(num+1)) *vect[bi]*np.conjugate(vect[bj])
				
	return cdc


def measureNk(basis,vect,L):
	cdc = measureCdiCj(basis,vect,L)
	nk = []
	kx = []
	for k in range(L):
		nxk=0.
		xk = -np.pi + k * 2.*np.pi/L;
		kx.append(xk)
		for i in range(L):
			for j in range(L):
				nxk += cdc[i][j]*np.exp(1j*xk*(j-i))
				
		if(abs(nxk.real/L)<1e-6 and abs(nxk.imag/L)<1e-6): nxk=0.
		if(abs(nxk.imag/L)<1e-6):
			nk.append(nxk.real/L)
		else:
			nk.append(nxk/L)
	return kx,nk

def measureNNk(basis,vect,L):
	ninj = measureNiNj(basis,vect,L)
	print("NiNj: ")
	print(ninj)
	nnk = []
	kx = []
	for k in range(L):
		nnxk=0.
		xk = -np.pi + k * 2.*np.pi/L;
		kx.append(xk)
		for i in range(L):
			for j in range(L):
				nnxk += (ninj[i][j] - ninj[i][i]*ninj[j][j]) *np.exp(1j*xk*(j-i))
		if(abs(nnxk.real/L)<1e-6): nnxk=0.
		nnk.append(nnxk.real/L)
	return kx,nnk



basis = build_all_basis(L, Nocu)

eigenval,eigenvect = spinless_r_space_V4_rise_k(basis,V4,mu)
print("E: ",eigenval)
eigenvect0 = np.transpose(eigenvect)[0]
print("E0_vec: ",eigenvect0)

print("Nk: ",measureNk(basis,eigenvect0,L))
print("NN(K): ",measureNNk(basis,eigenvect0,L))








