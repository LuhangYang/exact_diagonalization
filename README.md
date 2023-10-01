# Exact Diagonalization 

This library can be used to solve various problems in condensed matter physics using exact diagonalization method. This method consists of a series of different algorithms, which together can be a powerful tool to solve diverse problems in low energy physics. This library is in python. Most of the functions and classes defined are self-explanatory and is a good resource for the purpose of both teaching and conducting research.

#

## Features implemented:



  - **[Heisenberg 1/2 model](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_heis_RKKY)**

    1.  Both the original Heisenberg model and the one with RKKY-type of **long range interactions** are realized.
    2.  Both **real** and **momentum** space solutions can be obtained. The translational symmetry can be represented by the momentum, by using the good quantum number 'k'. We can get the spectrum as well as spectral weights for each momentum sector.
    3. The **Dynamical Structure Factor** or Green's function is calculated by the methods of both Lehmann representation and Correction Vectors.
    4. **Time evolution** is realized by Krylov formula.
    5. Other features such as basis truncation, correlation function measurements are also implmented.

  - **[Fermi-Hubbard model](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_hubbard)**

    1. Both in **real** space and **momentum** space. The implementation of momentum resolution is different with the Heisenberg model. In Hubbard model we divide the basis into different momentum sectors and compute directly with Hubbard Hamiltonian in the momentum space. For this purpose, I have a generic code for **any ranged hopping and interaction terms**.
    2. In order to fully utilize the parallellization capability of computers, I also implemented this part of ED in a **[bipartite formalism](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_bipartite_formalism)**
    3. **Dynamical Structure Factor** is obtained by solving the correction vectors.
    4. Since the Halmitonian terms are allowed to be in any distance, both **1D and 2D lattices** can be calculated.

- **[Spinless fermion model](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_spinless_fermion)**

     1. Tight binding model for spinless fermions with **long range interaction** terms.
     2. This code also includes the functions that measure several kinds of **correlations** in ground states.
    
- **[tJ model](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_tj)**

    1. This code can be used to calculate the tJ model ground states in **real space** and **momentum space**.
    2. I develope this code in order to understand the physics in the single-hole scenario.

- **[eta pairing states](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_SC_eta_pairing)**

    1. This code can be used to study the eta pairing states in extended Hubbard model. This code was developed in the purpose of studying some properties of superconductivity.

- **[Spinon dynamics using Villain approximation](https://github.com/LuhangYang/exact_diagonalization/tree/main/Spinon_time_evolving)**[1]

    1. Following the Villain's approximation we can study the spinons propagation and velocities.
    2. **Dynamic programming** technique is used to do the truncation of basis in a fast way.

- **[Three band Hubbard model](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_three_band_Hubbard)**
    1. Calculate the energy splitting between singlet and triplet states, and the energy difference between bonding and anti-bonding states for **Cu2O7** and **Cu2O8** clusters to get the effective t-t'-J-J' model coefficients.
    2. This code is developed to understand the low energy physics of three band Hubbard model.






    


[1] J.Villain,Propagative spin relaxation in the Ising-like antiferromagnetic linear chain,Phys- ica B+C 79,1(1975).
