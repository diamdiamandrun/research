def yamaguchi(mode=None,e_bs=None,e_trip=None,s2_bs=None,s2_trip=None,n_hono=None,n_luno=None):
    "Input all energy values in Hartrees. The code automatically converts them to kcal/mol or eV."
    if mode == None:
        print('Please input mode! Available modes include "gap" and "y_n", where n is arbitrary.')
        return
    elif "gap" in mode:
        if None not in (e_bs, e_trip, s2_bs, s2_trip):
            if s2_trip != s2_bs:

                import numpy as np
                # defining constants
                h = 6.62607015e-34  # J s
                a_0 = 5.29177210903e-11 # m
                hbar = h/(2*np.pi)  # J s
                e = 1.602176634e-19  # C
                m_e = 9.10938356e-31  # kg
                N_A = 6.02214076e23 # mol-1
                H = (hbar)**2/(m_e*a_0**2) # Hartrees

                j = (e_bs-e_trip)/(s2_trip-s2_bs) # Hartrees, exchange energy
                gap = (s2_trip)*j # Hartrees

                if "ev" in mode.lower():
                    HeV = H / e # Hartrees in eV
                    print(f'\u0394E_(S-T) = {HeV*gap:.3f} eV')
                else:
                    Hkcal = H*N_A/4184 # Hartrees in kcal/mol
                    print(f'\u0394E_(S-T) = {Hkcal*gap:.3f} kcal mol-1')             
            else:
                print('Check your inputs! s2_trip cannot be equal to s2_bs!')
        else:
            print("Check your inputs! For 'gap' mode, e_bs, e_trip, s2_bs & s2_trip are required!")
            return 
    elif 'y_' in mode:
        if None not in (n_hono, n_luno):
            T = 1/2 * (n_hono - n_luno) # "radical bond order"
            y_n = 1-(2*T)/(1+T**2)
            print(f'{mode} = {y_n:.3f}')
        else:
            print(f"Check your inputs! For '{mode}' mode, n_hono & n_luno are required!")
