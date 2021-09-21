# exact_diagonalization

This library can be used to calculate the following 1D models:

  - **[Heisenberg model with and without RKKY type long range interactions](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_heis_RKKY)**
    Both in real ann momentum space. The 1d Heisenberg model with PBC have translation symmetry, which can be represented by the momentum, by using the quantum number 'k', we can get the spectrum as well as spectral weights for every momentum sector.
    

  - **[Fermi-Hubbard model](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_hubbard)**

    Both in real space and momentum space. In momentum space every term in the Hamiltonian satisfies momentum conservation, which can be used to calculate the spectrum for every well-defined momentum sector.

- **[Spinless fermion model](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_spinless_fermion)**

     Tight binding model for spinless fermions with additional interacting terms.  This code also includes the functions that measure several kinds of correlations in ground states.
    
- **[tJ model](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_tj)**

    This code can be used to calculate the tJ model ground states in real space.

- **[eta pairing states](https://github.com/LuhangYang/exact_diagonalization/tree/main/ED_SC_eta_pairing)**

    This code can be used to study the eta pairing states in extended Hubbard model.

- **[Spinon dynamics using Villain approximation](https://github.com/LuhangYang/exact_diagonalization/tree/main/Spinon_time_evolving)**[1]

    Following the Villain's approximation we can study the spinons propagation and velocities.







    


[1] J.Villain,Propagative spin relaxation in the Ising-like antiferromagnetic linear chain,Phys- ica B+C 79,1(1975).
