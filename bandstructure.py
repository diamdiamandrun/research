from imports import *
from dataloading import *

# standard bandstructure plotting
def bandstructure(tdos=None, filename=None, high_symm_points=None, symm_points_name=None, material=None, colour=None, sigma=None, yscale=15, yspan=None, saveformat=None, titles=None):
    
    # Set default color if not provided
    if colour is None:
        colour = 'b'

    if plt.rcParams['figure.dpi'] == 100:
        size = 'small'
    elif plt.rcParams['figure.dpi'] > 100:
        size = 'big'
    else:
        size = None
    
    if size == 'big':
        font = 20
    else:
        font = None

    if tdos == 'On':
        try:
            data = loadexcel(f'{filename}_band')
            if data is None:
                raise FileNotFoundError
        except FileNotFoundError:
            data = copypaste()
            saveexcel(data, f'{filename}_band')
            print(f'"{filename}_band.csv" has been created. Process file in Excel to an Nx2 column and rerun the function.')
            return

        try:
            dosdata = loadexcel(f'{filename}_dos')
            if dosdata is None:
                raise FileNotFoundError
        except FileNotFoundError:
            dosdata = copypaste()
            saveexcel(dosdata, f'{filename}_dos')
            print(f'{filename}_dos.csv has been created. Process file in Excel by adding normalizing Fermi Energy to 0 eV and rerun the function.')
            return

        band = plotarray(data)
        dos = plotarray(dosdata)

        # reshape band.csv
        if len(band) == 2:
            df = pd.read_csv(f'{filename}_band.csv', header=None, names=['k', 'E-Ef'])
            pivot_df = df.pivot_table(index='k', columns=df.groupby('k').cumcount(), values='E-Ef', aggfunc='first')
            pivot_df.reset_index(inplace=True)
            saveexcel(pivot_df, f'{filename}_band')
            data = loadexcel(f'{filename}_band')
            band = plotarray(data)

        # create subplots
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios': [2, 1]})

        # band structure (left plot)
        axes[0].axhline(0, color='gray', linestyle='-', linewidth=0.5)

        if high_symm_points is not None and isinstance(high_symm_points, list):
            for point in high_symm_points:
                axes[0].axvline(point, color='gray', linestyle='-', linewidth=0.5)

        if symm_points_name is not None:
            axes[0].set_xticks(high_symm_points)
            axes[0].set_xticklabels(symm_points_name)

        for i in range(1, len(band)):
            # FIX: Prevent duplicate legend labels by only labelling the first line
            axes[0].plot(band[0], band[i], f'{colour}-', label=f'{material}' if i == 1 else "")
            
        axes[0].legend(loc='best', fontsize=.6*font if font else None)

        axes[0].set_xlabel('$k$-points', fontsize=font)
        axes[0].set_ylabel('$E - E_F$ (eV)', fontsize=font)
        axes[0].set_xlim(min(band[0]), max(band[0]))

        # total DOS (right plot)
        axes[1].tick_params(left=False)  # hide y-axis ticks
        axes[1].set_yticks([])  # remove y-axis labels
        axes[1].plot(dos[1], dos[0], f'{colour}-')  
        axes[1].set_xlabel('DOS (states eV$^{-1}$ cell$^{-1}$)', fontsize=font)
        axes[1].set_xlim(0)
        axes[1].axhline(0, color='gray', linestyle='-', linewidth=0.5)

        if type(yspan) != tuple:
            yspan = (-1*yscale, yscale)
        
        axes[0].set_ylim(yspan)
        axes[1].set_ylim(yspan)
        plt.tight_layout()
        
        if titles == 'On':
            axes[0].set_title(f'Band Structure of {material}', fontsize=font)
            axes[1].set_title("Total DOS", fontsize=font)

        # turn on for bandgap analysis
        if sigma is not None:
            homos, lumos, homosk, lumosk = [], [], [], []
            # eV, tolerance
            for i in range(1, len(band)):
                for j in range(len(band[i])):
                    if band[i][j] <= 0:
                        E, k = band[i][j], band[0][j]
                        if E >= np.max(band[i]) - sigma:
                            homos.append(E)
                            homosk.append(k)
                    elif band[i][j] >= 0:
                        E, k = band[i][j], band[0][j]
                        if E <= np.min(band[i]) + sigma:
                            lumos.append(E)
                            lumosk.append(k)
            
            if not homos or not lumos:
                print("Metallic system or no clear gap found within sigma tolerance.")
            else:
                # FIX: Add abs() to k-point distance check to prevent false positives for Dirac points
                if np.min(lumos) - np.max(homos) <= 10*sigma and abs(lumosk[np.argmin(lumos)] - homosk[np.argmax(homos)]) <= sigma:
                    axes[0].scatter(homosk[np.argmax(homos)], np.max(homos), marker='o', color='r')
                    axes[0].annotate(f"$k$ = {homosk[np.argmax(homos)]:.4f}", xy=[homosk[np.argmax(homos)] + 0.1, np.max(homos) + 0.4])
                    print(f'Dirac Point occurs at k = {homosk[np.argmax(homos)]}')

                else:
                    # encourage homo
                    dirEg, dirEgk, dirEgE = [], [], []
                    for i in range(len(homos)):
                        k1 = homosk[i]
                        for j in range(len(lumos)):
                            k2 = lumosk[j]
                            if k1 == k2 and homos[i] + sigma >= np.max(homos) and lumos[j] - sigma <= np.min(lumos):
                                dirEg.append(lumos[j] - homos[i])
                                dirEgk.append([homosk[i], lumosk[j]])
                                dirEgE.append([homos[i], lumos[j]])

                    # plot
                    if dirEg: # FIX: Pythonic empty list check
                        x, y = dirEgk[np.argmin(dirEg)], dirEgE[np.argmin(dirEg)]
                        axes[0].scatter(x, y, marker='o', color='r')
                        
                        if size is None or size == "small":
                            axes[0].annotate(f"({x[0]:.4f}, {y[0]:.4f})", xy=[x[0] + 0.01, y[0] + 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                            axes[0].annotate(f"({x[1]:.4f}, {y[1]:.4f})", xy=[x[1] + 0.01, y[1] - 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                        elif size == "big":
                            axes[0].annotate(f"({x[0]:.4f}, {y[0]:.4f})", xy=[x[0] + 0.01, y[0] + 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale], fontsize=15)
                            axes[0].annotate(f"({x[1]:.4f}, {y[1]:.4f})", xy=[x[1] + 0.01, y[1] - 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale], fontsize=15)
                        print(f'The direct bandgap, $E_g$ = {np.min(dirEg):.3f} eV')
                    else:
                        indirEgk = [homosk[np.argmax(homos)], lumosk[np.argmin(lumos)]]
                        indirEgE = [np.max(homos), np.min(lumos)]
                        axes[0].scatter(indirEgk, indirEgE, marker='o', color='r')
                        
                        if size is None or size == "small":
                            axes[0].annotate(f"({indirEgk[0]:.4f},{indirEgE[0]:.4f})", xy=[indirEgk[0] + 0.01, indirEgE[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                            axes[0].annotate(f"({indirEgk[1]:.4f},{indirEgE[1]:.4f})", xy=[indirEgk[1] + 0.01, indirEgE[1] - 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                        elif size == "big":
                            axes[0].annotate(f"({indirEgk[0]:.4f},{indirEgE[0]:.4f})", xy=[indirEgk[0] + 0.01, indirEgE[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale], fontsize=15)
                            axes[0].annotate(f"({indirEgk[1]:.4f},{indirEgE[1]:.4f})", xy=[indirEgk[1] + 0.01, indirEgE[1] - 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale], fontsize=15)
                        print(f'The indirect bandgap, $E_g$ = {indirEgE[1] - indirEgE[0]:.3f} eV')

    elif tdos is None or tdos == 'Off':
        try:
            data = loadexcel(f'{filename}_band')
            if data is None:
                raise FileNotFoundError
        except FileNotFoundError:
            data = copypaste()
            saveexcel(data, f'{filename}_band')
            print(f'"{filename}_band.csv" has been created. Process file in Excel to an Nx2 column and rerun the function.')
            return

        band = plotarray(data)

        if len(band) == 2:
            df = pd.read_csv(f'{filename}_band.csv', header=None, names=['k', 'E-Ef'])
            pivot_df = df.pivot_table(index='k', columns=df.groupby('k').cumcount(), values='E-Ef', aggfunc='first')
            pivot_df.reset_index(inplace=True)
            saveexcel(pivot_df, f'{filename}_band')
            data = loadexcel(f'{filename}_band')
            band = plotarray(data)

        plt.axhline(0, color='gray', linestyle='-', linewidth=0.5)

        if high_symm_points is not None and isinstance(high_symm_points, list):
            for point in high_symm_points:
                plt.axvline(point, color='gray', linestyle='-', linewidth=0.5)
        
        if symm_points_name is not None:
            plt.xticks(high_symm_points, symm_points_name)

        for i in range(1, len(band)):
            # FIX: Prevent duplicate legend labels
            plt.plot(band[0], band[i], f'{colour}-', label=f'{material}' if i == 1 else "")
            
        plt.legend(loc='best', fontsize=.6*font if font else None)

        if titles == 'On':
            plt.title(f'Band Structure of {material}', fontsize=font)
        plt.xlabel('$k$-points', fontsize=font)
        plt.ylabel('$E - E_F$ (eV)', fontsize=font)
        plt.xlim(min(band[0]), max(band[0]))
        plt.tight_layout()

        if type(yspan) != tuple:
            yspan = (-1*yscale, yscale)
        
        plt.ylim(yspan)

        # turn on for bandgap analysis
        if sigma is not None:
            homos, lumos, homosk, lumosk = [], [], [], []
            # eV, tolerance
            for i in range(1, len(band)):
                for j in range(len(band[i])):
                    if band[i][j] <= 0:
                        E, k = band[i][j], band[0][j]
                        if E >= np.max(band[i]) - sigma:
                            homos.append(E)
                            homosk.append(k)
                    elif band[i][j] >= 0:
                        E, k = band[i][j], band[0][j]
                        if E <= np.min(band[i]) + sigma:
                            lumos.append(E)
                            lumosk.append(k)
            
            if not homos or not lumos:
                print("Metallic system or no clear gap found within sigma tolerance.")
            else:
                if np.min(lumos) - np.max(homos) <= 10*sigma and abs(lumosk[np.argmin(lumos)] - homosk[np.argmax(homos)]) <= sigma:
                    plt.scatter(homosk[np.argmax(homos)], np.max(homos), marker='o', color='r')
                    plt.annotate(f"$k$ = {homosk[np.argmax(homos)]:.4f}", xy=[homosk[np.argmax(homos)] + 0.1, np.max(homos) + 0.4])
                    print(f'Dirac Point occurs at k = {homosk[np.argmax(homos)]}')

                else:
                    # encourage homo
                    dirEg, dirEgk, dirEgE = [], [], []
                    for i in range(len(homos)):
                        k1 = homosk[i]
                        for j in range(len(lumos)):
                            k2 = lumosk[j]
                            if k1 == k2 and homos[i] + sigma >= np.max(homos) and lumos[j] - sigma <= np.min(lumos):
                                dirEg.append(lumos[j] - homos[i])
                                dirEgk.append([homosk[i], lumosk[j]])
                                dirEgE.append([homos[i], lumos[j]])

                    # plot
                    if dirEg: # FIX: Pythonic empty list check
                        x, y = dirEgk[np.argmin(dirEg)], dirEgE[np.argmin(dirEg)]
                        plt.scatter(x, y, marker='o', color='r')
                        
                        if size is None or size == "small":
                            plt.annotate(f"({x[0]:.4f}, {y[0]:.4f})", xy=[x[0] + 0.01, y[0] + 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                            plt.annotate(f"({x[1]:.4f}, {y[1]:.4f})", xy=[x[1] + 0.01, y[1] - 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                        elif size == "big":
                            plt.annotate(f"({x[0]:.4f}, {y[0]:.4f})", xy=[x[0] + 0.01, y[0] + 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale], fontsize=15)
                            plt.annotate(f"({x[1]:.4f}, {y[1]:.4f})", xy=[x[1] + 0.01, y[1] - 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale], fontsize=15)
                        print(f'The direct bandgap, $E_g$ = {np.min(dirEg):.3f} eV')
                    else:
                        indirEgk = [homosk[np.argmax(homos)], lumosk[np.argmin(lumos)]]
                        indirEgE = [np.max(homos), np.min(lumos)]
                        plt.scatter(indirEgk, indirEgE, marker='o', color='r')
                        
                        if size is None or size == "small":
                            plt.annotate(f"({indirEgk[0]:.4f},{indirEgE[0]:.4f})", xy=[indirEgk[0] + 0.01, indirEgE[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                            plt.annotate(f"({indirEgk[1]:.4f},{indirEgE[1]:.4f})", xy=[indirEgk[1] + 0.01, indirEgE[1] - 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                        elif size == "big":
                            plt.annotate(f"({indirEgk[0]:.4f},{indirEgE[0]:.4f})", xy=[indirEgk[0] + 0.01, indirEgE[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale], fontsize=15)
                            plt.annotate(f"({indirEgk[1]:.4f},{indirEgE[1]:.4f})", xy=[indirEgk[1] + 0.01, indirEgE[1] - 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale], fontsize=15)
                        print(f'The indirect bandgap, $E_g$ = {indirEgE[1] - indirEgE[0]:.3f} eV')

    if saveformat is not None:
        if not os.path.exists('00_Plots'):
            os.makedirs('00_Plots')
        plt.savefig(f'00_Plots/{filename}.{saveformat}', bbox_inches='tight')
        print(f"Plot saved as {filename}.{saveformat} in 00_Plots")


# spin-split bandstructure plotting, collapses back to bandstructure() for bandstructure_ss(spin_split=False)
def bandstructure_ss(tdos=None, filename=None, high_symm_points=None, symm_points_name=None, material=None, colour=None, sigma=0.01, yscale=15, yspan=None, saveformat=None, titles=None, spin_split=False, colour_up='r', colour_down='b', ls_up='-', ls_down='--', elements=None, element_colours=None):
    
    # Prevent color override
    if colour is None:
        colour = 'b'
        
    if plt.rcParams['figure.dpi'] == 100:
        size='small'
    elif plt.rcParams['figure.dpi'] > 100:
        size='big'
    else:
        size=None
    
    if size=='big':
        font=20
    else:
        font=None

    # load band data
    try:
        data = loadexcel(f'{filename}_band')
        if data is None:
            raise FileNotFoundError
    except FileNotFoundError:
        data = copypaste()
        saveexcel(data, f'{filename}_band')
        print(f'"{filename}_band.csv" has been created. Process file in Excel to an Nx2 (or Nx3 for spin split) column and rerun the function.')
        return

    band = plotarray(data)

    # load dos data using the elements parameter when tdos is active
    dos_list = []
    pdos_labels = []
    pdos_colours = []
    
    if tdos == 'On' or tdos is True:
        if elements is not None:
            for idx, el in enumerate(elements):
                tag = f"{filename}_dos_{el}"
                try:
                    dosdata = loadexcel(tag)
                    if dosdata is None:
                        raise FileNotFoundError
                except FileNotFoundError:
                    dosdata = copypaste()
                    saveexcel(dosdata, tag)
                    print(f'"{tag}.csv" has been created. Process file and rerun the function.')
                    return
                dos_list.append(plotarray(dosdata))
                pdos_labels.append(el)
                
                # Assign custom color if provided, else fall back safely
                if element_colours is not None and idx < len(element_colours):
                    pdos_colours.append(element_colours[idx])
                else:
                    pdos_colours.append(colour)
        else:
            # Fallback legacy single file loading if elements list is omitted
            try:
                dosdata = loadexcel(f'{filename}_dos')
                if dosdata is None:
                    raise FileNotFoundError
            except FileNotFoundError:
                dosdata = copypaste()
                saveexcel(dosdata, f'{filename}_dos')
                print(f'{filename}_dos.csv has been created. Process file and rerun.')
                return
            dos_list.append(plotarray(dosdata))
            pdos_labels.append('Total')
            pdos_colours.append(colour)

    # reshape band data
    if spin_split and os.path.exists(f'{filename}_band_up.csv') and os.path.exists(f'{filename}_band_down.csv'):
        data_up = loadexcel(f'{filename}_band_up')
        band_up = plotarray(data_up)
        data_down = loadexcel(f'{filename}_band_down')
        band_down = plotarray(data_down)
    elif len(band) in [2, 3]:
#================== spin splitting processing ===============
        if spin_split:
            df = pd.read_csv(f'{filename}_band.csv', header=None, names=['k', 'E-Ef', 'Spin'], sep=None, engine='python')
                
            df_up = df[df['Spin'] == 1].copy().drop(columns=['Spin'])                
            df_down = df[df['Spin'] == 2].copy().drop(columns=['Spin'])

            def process_spin_df(sub_df):
                sub_df['band_id'] = (sub_df['k'] < sub_df['k'].shift(1)).cumsum()
                pivot_df = sub_df.pivot_table(index='k', columns='band_id', values='E-Ef', aggfunc='first')
                return pivot_df.reset_index()

            pivot_up = process_spin_df(df_up)
            saveexcel(pivot_up, f'{filename}_band_up')
            data_up = loadexcel(f'{filename}_band_up')
            band_up = plotarray(data_up)

            pivot_down = process_spin_df(df_down)
            saveexcel(pivot_down, f'{filename}_band_down')
            data_down = loadexcel(f'{filename}_band_down')
            band_down = plotarray(data_down)
#==============================================================
        else:
            df = pd.read_csv(f'{filename}_band.csv', header=None, usecols=[0, 1], names=['k', 'E-Ef'], sep=None, engine='python')
            pivot_df = df.pivot_table(index='k', columns=df.groupby('k').cumcount(), values='E-Ef', aggfunc='first').reset_index()
            saveexcel(pivot_df, f'{filename}_band')
            data = loadexcel(f'{filename}_band')
            band = plotarray(data)
    else:
        if spin_split:
            try:
                data_up = loadexcel(f'{filename}_band_up')
                band_up = plotarray(data_up)
                data_down = loadexcel(f'{filename}_band_down')
                band_down = plotarray(data_down)
            except FileNotFoundError:
                print(f'"{filename}_band_up.csv" OR DOWN NOT FOUND, UNABLE TO PLOT SPIN SPLIT WITHOUT RAW DATA.')
                return

    # initiate plots
    if tdos == 'On' or tdos is True:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios': [2, 1]})
        ax_band = axes[0]
        ax_dos = axes[1]
    else:
        fig, ax_band = plt.subplots(figsize=(7, 5))
        ax_dos = None

    # bandstructure key components
    ax_band.axhline(0, color='gray', linestyle='-', linewidth=0.5)

    if high_symm_points is not None and isinstance(high_symm_points, list):
        for point in high_symm_points:
            ax_band.axvline(point, color='gray', linestyle='-', linewidth=0.5)

    if symm_points_name is not None:
        ax_band.set_xticks(high_symm_points)
        ax_band.set_xticklabels(symm_points_name)

    # plot bandstructure
    if spin_split:
        for i in range(1, len(band_up)):
            ax_band.plot(band_up[0], band_up[i], color=colour_up, linestyle=ls_up, label=f'{material} (up)' if i == 1 else "")
        for i in range(1, len(band_down)):
            ax_band.plot(band_down[0], band_down[i], color=colour_down, linestyle=ls_down, label=f'{material} (down)' if i == 1 else "")
        ax_band.legend(loc='best', fontsize=.6*font if font else None)
        ax_band.set_xlim(min(band_up[0]), max(band_up[0]))
    else:
        for i in range(1, len(band)):
            ax_band.plot(band[0], band[i], f'{colour}-', label=f'{material}')
        ax_band.set_xlim(min(band[0]), max(band[0]))

    ax_band.set_ylabel('$E - E_F$ (eV)', fontsize=font)
    
    if type(yspan) != tuple:
        yspan = (-1*yscale, yscale)
    ax_band.set_ylim(yspan)

    if titles == 'On':
        ax_band.set_title(f'Band Structure of {material}', fontsize=font)
        ax_band.set_xlabel('$k$-points', fontsize=font)

    # dos key components
    if ax_dos is not None:
        ax_dos.tick_params(left=False)
        ax_dos.set_yticks([])
        
        for idx, dos in enumerate(dos_list):
            c = pdos_colours[idx]
            lbl = pdos_labels[idx]
            
            if spin_split and len(dos) >= 3:
                ax_dos.plot(dos[1], dos[0], color=c, linestyle=ls_up, label=f'{lbl} (up)')
                neg_dos_down = [-x for x in dos[2]]
                ax_dos.plot(neg_dos_down, dos[0], color=c, linestyle=ls_down, label=f'{lbl} (down)')
            else:
                ax_dos.plot(dos[1], dos[0], color=c, linestyle='-', label=lbl)
        
        ax_dos.legend(loc='best', fontsize=0.6*font if font else None)
        
        all_max = max([max(d[1]) for d in dos_list])
        if spin_split:
            all_min = min([min([-x for x in d[2]]) for d in dos_list if len(d) >= 3])
            ax_dos.set_xlim(left=all_min*1.1, right=all_max*1.1)
        else:
            ax_dos.set_xlim(0, all_max*1.1)
            
        ax_dos.set_xlabel('DOS (states eV$^{-1}$ cell$^{-1}$)', fontsize=font)
        ax_dos.axhline(0, color='gray', linestyle='-', linewidth=0.5)
        ax_dos.set_ylim(yspan)
        
        if titles == 'On':
            ax_dos.set_title("Projected DOS", fontsize=font)

    plt.tight_layout()

    # bandgap analysis
    if sigma is not None:
        def analyze_bandgap(bnd, label=""):
            homos, lumos, homosk, lumosk = [], [], [], []
            
            for i in range(1, len(bnd)):
                bnd_max = np.max(bnd[i])
                bnd_min = np.min(bnd[i])
                
                for j in range(len(bnd[i])):
                    E, k = bnd[i][j], bnd[0][j]
                    if E <= 0:
                        if E >= bnd_max - sigma:
                            homos.append(E)
                            homosk.append(k)
                    elif E >= 0:
                        if E <= bnd_min + sigma:
                            lumos.append(E)
                            lumosk.append(k)
            
            if not homos or not lumos:
                print(f"Metallic system or no clear gap found within sigma tolerance {label}.")
                return
            
            global_max_homo = np.max(homos)
            global_min_lumo = np.min(lumos)
            marker_color = colour_up if label == "(up)" else (colour_down if label == "(down)" else 'r')
            
            if global_min_lumo - global_max_homo <= 10*sigma and abs(lumosk[np.argmin(lumos)] - homosk[np.argmax(homos)]) <= sigma:
                ax_band.scatter(homosk[np.argmax(homos)], global_max_homo, marker='o', color=marker_color)
                ax_band.annotate(f"$k$ = {homosk[np.argmax(homos)]:.4f} {label}", xy=[homosk[np.argmax(homos)] + 0.1, global_max_homo + 0.4])
                print(f'Dirac Point occurs at k = {homosk[np.argmax(homos)]} {label}')
            else:
                dirEg, dirEgk, dirEgE = [], [], []
                
                valid_homos = [(E, k) for E, k in zip(homos, homosk) if E + sigma >= global_max_homo]
                valid_lumos = [(E, k) for E, k in zip(lumos, lumosk) if E - sigma <= global_min_lumo]
                
                for h_E, h_k in valid_homos:
                    for l_E, l_k in valid_lumos:
                        if abs(h_k - l_k) <= sigma:
                            dirEg.append(l_E - h_E)
                            dirEgk.append([h_k, l_k])
                            dirEgE.append([h_E, l_E])
                
                if dirEg:
                    min_idx = np.argmin(dirEg)
                    x, y = dirEgk[min_idx], dirEgE[min_idx]
                    ax_band.scatter(x, y, marker='o', color=marker_color)
                    
                    fs = 15 if size == "big" else None
                    ax_band.annotate(f"({x[0]:.2f}, {y[0]:.2f}) {label}", xy=[x[0] + 0.01, y[0] + 0.03 * (yspan[1] - yspan[0])/yscale], fontsize=fs)
                    ax_band.annotate(f"({x[1]:.2f}, {y[1]:.2f}) {label}", xy=[x[1] + 0.01, y[1] - 0.03 * (yspan[1] - yspan[0])/yscale], fontsize=fs)
                    print(f'The direct bandgap, $E_g$ = {np.min(dirEg):.3f} eV {label}')
                else:
                    indirEgk = [homosk[np.argmax(homos)], lumosk[np.argmin(lumos)]]
                    indirEgE = [global_max_homo, global_min_lumo]
                    ax_band.scatter(indirEgk, indirEgE, marker='o', color=marker_color)
                    
                    fs = 15 if size == "big" else None
                    ax_band.annotate(f"({indirEgk[0]:.2f},{indirEgE[0]:.2f}) {label}", xy=[indirEgk[0] + 0.01, indirEgE[0] + 0.3 * (yspan[1] - yspan[0])/yscale], fontsize=fs)
                    ax_band.annotate(f"({indirEgk[1]:.2f},{indirEgE[1]:.2f}) {label}", xy=[indirEgk[1] + 0.01, indirEgE[1] - 0.3 * (yspan[1] - yspan[0])/yscale], fontsize=fs)
                    print(f'The indirect bandgap, $E_g$ = {indirEgE[1] - indirEgE[0]:.3f} eV {label}')

        if spin_split:
            analyze_bandgap(band_up, label="(up)")
            analyze_bandgap(band_down, label="(down)")
        else:
            analyze_bandgap(band)

    if saveformat is not None:
        if not os.path.exists('00_Plots'):
            os.makedirs('00_Plots')
        plt.savefig(f'00_Plots/{filename}.{saveformat}', bbox_inches='tight')
        print(f"Plot saved as {filename}.{saveformat} in 00_Plots")
