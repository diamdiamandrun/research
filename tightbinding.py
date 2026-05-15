import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

k = sp.symbols('k', real=True)

def HMOnnn(A, plot=False, k_val=0.0, nnn=False, t_nnn=-0.2):
    """
    Solve a Hückel (tight-binding) Hamiltonian by direct diagonalisation.
    
    Parameters
    ----------
    A : list of lists (symbolic)
        Hückel matrix with symbolic k
    k_val : float
        Numerical value of k
    plot : bool
        Plot energy levels if True
    nnn : bool
        If True, introduces Next-Nearest-Neighbour (NNN) hopping.
    t_nnn : float
        The hopping amplitude for NNN interactions. (Default assumes beta ~ -1).
        
    Returns
    -------
    en : ndarray
        Sorted real eigenvalues (all levels, including degeneracies)
    """

    # convert to sympy matrix for symbolic manipulation
    M = sp.Matrix(A)
    N = M.shape[0]
    
    # nnn hopping; adds t_nnn to the second off-diagonal elements (i, i+2) and (i+2, i) with pbc
    if nnn and N > 2:
        added_bonds = set()
        for i in range(N):
            j = (i + 2) % N
            bond = tuple(sorted((i, j)))
            
            # prevent double counting bonds
            if bond not in added_bonds:
                # Add t_nnn to the symmetric off-diagonal elements
                M[bond[0], bond[1]] += t_nnn
                M[bond[1], bond[0]] += t_nnn
                added_bonds.add(bond)

    # sub k for numerical value
    M_num = M.subs(k, k_val)

    # to numpy for diagonalisation, dtype=complex to preserve complex phases
    M_np = np.array(M_num.tolist(), dtype=complex)

    # diagonalise and ensure eigenvalues are real
    en = np.linalg.eigvalsh(M_np)

    # energy levels
    if plot:
        plt.figure(figsize=(8,5))
        plt.axhline(0, color='r', linestyle='--', linewidth=1)  # fermi level
        for i, E in enumerate(en):
            plt.hlines(E, xmin=i - 0.3, xmax=i + 0.3,
                       colors='b', linewidth=2)
        plt.xticks([])
        plt.xlabel("Orbital index")
        plt.ylabel("Energy, ($-\\beta$)")
        
        # title 
        title = f"Hückel spectrum for [{N}]Annulene (k = {k_val/np.pi:.1f} $\pi$)"
        
        if nnn:
            title += f", NNN on ($t'$ = {t_nnn})"
        plt.title(title)
        
        plt.tight_layout()
        plt.show()
    
    return en


def build_pah(N, beta=-1.0, delta=0.0):
    """
    Builds an N-membered Hückel matrix with a topological phase 'k' 
    on the boundary bond, and optional bond alternation (delta).
    """
    M = sp.zeros(N, N) # zeroes sq matrix

    # nearest neighbour (alternating if delta != 0)
    for i in range(N - 1):
        current_beta = (beta - delta) if i % 2 == 0 else (beta + delta)
        M[i, i+1] = current_beta
        M[i+1, i] = current_beta
    
    # boundary condition (maintain alternating logic)
    closing_beta = (beta - delta) if (N - 1) % 2 == 0 else (beta + delta)
    M[0, N-1] = closing_beta * sp.exp(-sp.I * k)
    M[N-1, 0] = closing_beta * sp.exp(sp.I * k)

    return M.tolist()
