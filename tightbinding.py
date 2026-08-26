## To experiment with vertex perturbation for Td and Oh cages

import itertools
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist

# coordinate generator for graph, euclidean coordinates required because for symmfinder engine

def generate_tetrahedron_coords_topological(nodes_per_edge=4):
  v_base = np.array(
      [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float
  )

  edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
  coords = [pt for pt in v_base]
  is_vertex = [True, True, True, True]

  for e in edges:
    v1, v2 = v_base[e[0]], v_base[e[1]]
    for i in range(1, nodes_per_edge + 1):
      frac = i / (nodes_per_edge + 1)
      coords.append(v1 + frac * (v2 - v1))
      is_vertex.append(False)

  coords = np.array(coords)
  coords -= np.mean(coords, axis=0)
  return coords, np.array(is_vertex, dtype=bool)


def generate_diynecube_coords_topological(nodes_per_edge=4):
  L = nodes_per_edge + 1
  coords = []
  is_vertex = []
  for x in range(L + 1):
    for y in range(L + 1):
      for z in range(L + 1):
        if (x in (0, L)) + (y in (0, L)) + (z in (0, L)) >= 2:
          coords.append([x, y, z])
          # Identify true vertices of the cube (corners where 3 edges meet)
          is_vertex.append((x in (0, L)) and (y in (0, L)) and (z in (0, L)))

  coords = np.array(coords, dtype=float)
  coords -= np.mean(coords, axis=0)
  return coords, np.array(is_vertex, dtype=bool)


# symmfinder engine

def get_Td_symm(V, coords):
  N, dim = V.shape
  irrep_names = ["A1", "A2", "E", "T1", "T2"]
  char_table = np.array([
      [1, 1, 1, 1, 1],
      [1, 1, 1, -1, -1],
      [2, -1, 2, 0, 0],
      [3, 0, -1, 1, -1],
      [3, 0, -1, -1, 1],
  ])

  def get_class_index(M):
    det = int(round(np.linalg.det(M)))
    tr = int(round(np.trace(M)))
    if det == 1:
      if tr == 3:
        return 0
      if tr == 0:
        return 1
      if tr == -1:
        return 2
    elif det == -1:
      if tr == -1:
        return 3
      if tr == 1:
        return 4
    raise ValueError(f"Unknown Td class for det={det}, tr={tr}")

  td_vertices = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
  ops = []

  for p in itertools.permutations([0, 1, 2]):
    for signs in itertools.product([1, -1], repeat=3):
      M = np.zeros((3, 3))
      for i in range(3):
        M[i, p[i]] = signs[i]
      new_v = td_vertices @ M.T
      if np.all(np.min(cdist(new_v, td_vertices), axis=1) < 1e-4):
        ops.append(M)

  chi_V_ops = np.zeros(24)
  class_indices = np.zeros(24, dtype=int)

  for op_idx, M in enumerate(ops):
    class_indices[op_idx] = get_class_index(M)
    new_coords = coords @ M.T
    j_indices = np.argmin(cdist(new_coords, coords), axis=1)

    P_mat = np.zeros((N, N))
    P_mat[j_indices, np.arange(N)] = 1.0
    chi_V_ops[op_idx] = np.trace(V.T @ P_mat @ V)

  label_cell = []
  for i, irrep in enumerate(irrep_names):
    count = int(np.round(np.sum(char_table[i, class_indices] * chi_V_ops) / 24.0))
    if count > 0:
      label_cell.append(irrep if count == 1 else f"{count}{irrep}")

  return " $\oplus$ ".join(label_cell) if label_cell else "Unknown"


def get_Oh_symm_cube(V, coords):
  N, dim = V.shape
  irrep_names = [
      "A1g",
      "A2g",
      "Eg",
      "T1g",
      "T2g",
      "A1u",
      "A2u",
      "Eu",
      "T1u",
      "T2u",
  ]
  char_table = np.array([
      [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      [1, 1, -1, -1, 1, 1, -1, 1, 1, -1],
      [2, -1, 0, 0, 2, 2, 0, -1, 2, 0],
      [3, 0, -1, 1, -1, 3, 1, 0, -1, -1],
      [3, 0, 1, -1, -1, 3, -1, 0, -1, 1],
      [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
      [1, 1, -1, -1, 1, -1, 1, -1, -1, 1],
      [2, -1, 0, 0, 2, -2, 0, 1, -2, 0],
      [3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
      [3, 0, 1, -1, -1, -3, 1, 0, 1, -1],
  ])

  def get_class_index(M):
    det, tr = int(round(np.linalg.det(M))), int(round(np.trace(M)))
    zeros = int(np.sum(np.isclose(np.diag(M), 0)))
    if det == 1:
      if tr == 3:
        return 0
      if tr == 0:
        return 1
      if tr == -1 and zeros == 2:
        return 2
      if tr == 1:
        return 3
      if tr == -1 and zeros == 0:
        return 4
    elif det == -1:
      if tr == -3:
        return 5
      if tr == -1:
        return 6
      if tr == 0:
        return 7
      if tr == 1 and zeros == 0:
        return 8
      if tr == 1 and zeros == 2:
        return 9
    raise ValueError(f"Unknown O_h class for det={det}, tr={tr}")

  chi_V_ops = np.zeros(48)
  class_indices = np.zeros(48, dtype=int)

  op_idx = 0
  for p in itertools.permutations([0, 1, 2]):
    for signs in itertools.product([1, -1], repeat=3):
      M = np.zeros((3, 3))
      for i in range(3):
        M[i, p[i]] = signs[i]

      class_indices[op_idx] = get_class_index(M)
      j_indices = np.argmin(cdist(coords @ M.T, coords), axis=1)

      P_mat = np.zeros((N, N))
      P_mat[j_indices, np.arange(N)] = 1.0
      chi_V_ops[op_idx] = np.trace(V.T @ P_mat @ V)
      op_idx += 1

  label_cell = []
  for i, irrep in enumerate(irrep_names):
    count = int(
        np.round(np.sum(char_table[i, class_indices] * chi_V_ops) / 48.0)
    )
    if count > 0:
      label_cell.append(irrep if count == 1 else f"{count}{irrep}")

  return " $\oplus$ ".join(label_cell) if label_cell else "Unknown"


# plotter

def draw_spectrum(en, vecs, coords, symm_finder, title, font=None):
  plt.figure(figsize=(10, 6))
  plt.axhline(0, color="r", linestyle="--", linewidth=1)

  tol = 1e-4
  current_idx = 0

  while current_idx < len(en):
    mask = np.abs(en - en[current_idx]) < tol
    count = np.sum(mask)

    symm = symm_finder(vecs[:, mask], coords) if symm_finder else "Unknown"
    center_x = current_idx + (count - 1) / 2.0

    plt.text(
        center_x + 0.2,
        np.mean(en[mask]) - 0.2,
        symm,
        ha="center",
        va="top",
        fontsize=10 if font is None else font - 3,
        color="darkred",
    )
    current_idx += count

  for i, E in enumerate(en):
    plt.hlines(E, xmin=i - 0.3, xmax=i + 0.3, colors="b", linewidth=2)

  plt.xticks([])
  if font:
    plt.yticks(fontsize=font)
  plt.xlabel("Orbital index", fontsize=font)
  plt.ylabel("Energy, $n$ in $\\alpha-n\\beta$", fontsize=font)
  plt.title(title, fontsize=font)
  plt.tight_layout()
  plt.show()


# spherical topological engine

def spherical_topological_HMO(
    coords,
    is_vertex,
    k_val=0.0,
    symm_finder=None,
    plot=False,
    title="",
    big=False,
):
  """Solves the Pure Topological Hückel model with Vertex On-Site Energy Perturbation (k)."""
  N = len(coords)
  D = cdist(coords, coords)

  # NN bondlength
  D_no_diag = D.copy()
  np.fill_diagonal(D_no_diag, np.inf)
  bond_length = np.min(D_no_diag)

  # adjacency matrix, set to -\beta for NN
  H = np.zeros((N, N))
  H[np.abs(D - bond_length) < 1e-4] = -1.0

  # onsite perturbation for vertices
  for i in range(N):
    if is_vertex[i]:
      H[i, i] = k_val

  # diagonalise and sort eigenvalues and eigenvectors
  en, vecs = np.linalg.eigh(H)
  en = -en  # stacking convention
  idx = np.argsort(en)
  en, vecs = en[idx], vecs[:, idx]

  if plot:
    draw_spectrum(
        en, vecs, coords, symm_finder, title, font=15 if big else None
    )

  return en, vecs

# function definition

def HMO_perturbed(pointgroup="Td", n=3, k=0.5):
  """Solves and plots the spherical topological Hückel model.

  If k is non-zero, it plots both the unperturbed (k = 0.0) and
  vertex-perturbed states using the original title convention.
  """
  # Select configuration based on point group string
  if pointgroup == "Td":
    coords, is_vertex = generate_tetrahedron_coords_topological(
        nodes_per_edge=n
    )
    symm_func = get_Td_symm
    title_name = "Tetrahedron (Td)"
  elif pointgroup == "Oh":
    coords, is_vertex = generate_diynecube_coords_topological(nodes_per_edge=n)
    symm_func = get_Oh_symm_cube
    title_name = "Diynecube (Oh)"
  else:
    raise ValueError(f"Unsupported pointgroup: {pointgroup}")

  N = len(coords)
  D = cdist(coords, coords)

  # Find bond length from nearest neighbors
  D_no_diag = D.copy()
  np.fill_diagonal(D_no_diag, np.inf)
  bond_length = np.min(D_no_diag)

  # Build Adjacency Matrix (-1 for strict nearest neighbors)
  H_base = np.zeros((N, N))
  H_base[np.abs(D - bond_length) < 1e-4] = -1.0

  # --- 1. Unperturbed Solver & Plot ---
  en_base, vecs_base = np.linalg.eigh(H_base)
  en_base = -en_base
  idx_base = np.argsort(en_base)
  en_base, vecs_base = en_base[idx_base], vecs_base[:, idx_base]

  draw_spectrum(
      en_base,
      vecs_base,
      coords,
      symm_func,
      title=f"{title_name} Cage ($n={n}$), Unperturbed ($k = 0.0$)",
      font=15,
  )

  # --- 2. Perturbed Solver & Plot (if k is non-zero) ---
  if k != 0 and k is not False and k is not None:
    H_perturbed = H_base.copy()
    for i in range(N):
      if is_vertex[i]:
        H_perturbed[i, i] = k

    en_pert, vecs_pert = np.linalg.eigh(H_perturbed)
    en_pert = -en_pert
    idx_pert = np.argsort(en_pert)
    en_pert, vecs_pert = en_pert[idx_pert], vecs_pert[:, idx_pert]

    draw_spectrum(
        en_pert,
        vecs_pert,
        coords,
        symm_func,
        title=f"{title_name} Cage ($n={n}$), Vertex Perturbed ($k = {k}$)",
        font=15,
    )
    return (en_base, vecs_base), (en_pert, vecs_pert)

  return en_base, vecs_base


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################


# to explore NNN
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
