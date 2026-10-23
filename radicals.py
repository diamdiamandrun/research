import numpy as np

def yamaguchi(mode=None, e_ls=None, e_hs=None, s2_ls=None, s2_hs=None, s2_pure=0.0, n_hono=None, n_luno=None):
    "Input all energy values in Hartrees. The code automatically converts them to kcal/mol or eV."
    
    if mode is None:
        print('Please input mode! Available modes include "gap", "gap ev", and "y_n".')
        return
        
    if "gap" in mode.lower():
        if None not in (e_ls, e_hs, s2_ls, s2_hs):
            if s2_hs != s2_ls:

                # defining constants
                h = 6.62607015e-34  # J s
                a_0 = 5.29177210903e-11 # m
                hbar = h/(2*np.pi)  # J s
                e = 1.602176634e-19  # C
                m_e = 9.10938356e-31  # kg
                N_A = 6.02214076e23 # mol-1
                H = (hbar)**2/(m_e*a_0**2) # Hartrees in Joules

                j = (e_ls - e_hs) / (s2_hs - s2_ls) # Hartrees, exchange energy
                
                # gap
                if s2_pure == 0.0:
                    print('Be sure to check input for s2_pure. Set to 0.0 by default for broken singlets.')
                gap = j * (s2_hs - s2_pure) # Hartrees

                if "ev" in mode.lower():
                    HeV = H / e # Hartrees in eV
                    print(f'\u0394E_(HS-LS) = {HeV * gap:.3f} eV')
                else:
                    Hkcal = H * N_A / 4184 # Hartrees in kcal/mol
                    print(f'\u0394E_(HS-LS) = {Hkcal * gap:.3f} kcal mol-1')             
            else:
                print('Check your inputs! s2_hs cannot be equal to s2_ls!')
        else:
            print("Check your inputs! For 'gap' mode, e_ls, e_hs, s2_ls & s2_hs are required!")
            
    elif 'y_' in mode.lower():
        if None not in (n_hono, n_luno):
            T = 1/2 * (n_hono - n_luno) # "radical bond order"
            y_n = 1 - (2 * T) / (1 + T**2)
            print(f'{mode} = {y_n:.3f}')
        else:
            print(f"Check your inputs! For '{mode}' mode, n_hono & n_luno are required!")
    else:
        print('Unrecognized mode! Available modes include "gap", "gap ev", and "y_n".')
