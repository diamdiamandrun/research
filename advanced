def bandstructure(tdos=None, filename=None, high_symm_points=None, symm_points_name=None, material=None, colour=None, sigma=None, yscale=15, yspan=None, saveformat=None):
    
    # Set default color if not provided
    if colour is None:
        colour = 'b'
    
    if plt.rcParams['figure.dpi'] == 100:
        size='small'
    elif plt.rcParams['figure.dpi']== 1000:
        size='big'

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

        # band structure (left plot))
        axes[0].set_title(f'Band Structure of {material}')
        axes[0].axhline(0, color='gray', linestyle='-', linewidth=0.5)

        if high_symm_points is not None and isinstance(high_symm_points, list):
            for point in high_symm_points:
                axes[0].axvline(point, color='gray', linestyle='-', linewidth=0.5)

        if symm_points_name is not None:
            axes[0].set_xticks(high_symm_points)
            axes[0].set_xticklabels(symm_points_name)

        for i in range(1, len(band)):
            axes[0].plot(band[0], band[i], f'{colour}-', label=f'{material}')

        axes[0].set_xlabel('$k$-points')
        axes[0].set_ylabel('$E - E_F$ (eV)')
        axes[0].set_xlim(min(band[0]), max(band[0]))

        # total DOS (right plot)
        axes[1].set_title("Total DOS")
        axes[1].tick_params(left=False)  #hide y-axis ticks
        axes[1].set_yticks([])  #remove y-axis labels
        axes[1].plot(dos[1], dos[0], f'{colour}-')  
        axes[1].set_xlabel('DOS (states eV$^{-1}$ cell$^{-1}$)')
        axes[1].set_xlim(0)
        axes[1].axhline(0, color='gray', linestyle='-', linewidth=0.5)

        if type(yspan) != tuple:
            yspan = (-1*yscale,yscale)
        
        axes[0].set_ylim(yspan)
        axes[1].set_ylim(yspan)
        plt.tight_layout()

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
            
            #print(np.min(lumos) - np.max(homos),lumosk[np.argmin(lumos)]-homosk[np.argmax(homos)])
            
            # check for dirac points
            if np.min(lumos) - np.max(homos) <= 10*sigma and lumosk[np.argmin(lumos)]-homosk[np.argmax(homos)] <= sigma:
                axes[0].scatter(homosk[np.argmax(homos)],np.max(homos), marker = 'o', color = 'r')
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
                if dirEg != []:
                    x, y = dirEgk[np.argmin(dirEg)], dirEgE[np.argmin(dirEg)]
                    #print(x, y)
                    axes[0].scatter(x, y, marker='o', color='r')
                    
                    if size == None or size == "small":
                        axes[0].annotate(f"({x[0]:.4f}, {y[0]:.4f})", xy=[x[0] + 0.01, y[0] + 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                        axes[0].annotate(f"({x[1]:.4f}, {y[1]:.4f})", xy=[x[1] + 0.01, y[1] - 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                    elif size == "big":
                        axes[0].annotate(f"({x[0]:.4f}, {y[0]:.4f})", xy=[x[0] + 0.01, y[0] + 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale],fontsize=15)
                        axes[0].annotate(f"({x[1]:.4f}, {y[1]:.4f})", xy=[x[1] + 0.01, y[1] - 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale],fontsize=15)
                    print(f'The direct bandgap, $E_g$ = {np.min(dirEg):.3f} eV')
                else:
                    indirEgk = [homosk[np.argmax(homos)], lumosk[np.argmin(lumos)]]
                    indirEgE = [np.max(homos), np.min(lumos)]
                    axes[0].scatter(indirEgk, indirEgE, marker='o', color='r')
                    
                    if size is None or size == "small":
                        axes[0].annotate(f"({indirEgk[0]:.4f},{indirEgE[0]:.4f})", xy=[indirEgk[0] + 0.01, indirEgE[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                        axes[0].annotate(f"({indirEgk[1]:.4f},{indirEgE[1]:.4f})", xy=[indirEgk[1] + 0.01, indirEgE[1] - 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                    elif size == "big":
                        axes[0].annotate(f"({indirEgk[0]:.4f},{indirEgE[0]:.4f})", xy=[indirEgk[0] + 0.01, indirEgE[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale],fontsize=15)
                        axes[0].annotate(f"({indirEgk[1]:.4f},{indirEgE[1]:.4f})", xy=[indirEgk[1] + 0.01, indirEgE[1] - 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale],fontsize=15)
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
        
        plt.xticks(high_symm_points, symm_points_name)

        for i in range(1, len(band)):
            plt.plot(band[0], band[i], f'{colour}-', label=f'{material}')

        plt.xlabel('$k$-points')
        plt.ylabel('$E - E_F$ (eV)')
        plt.xlim(min(band[0]), max(band[0]))
        plt.tight_layout()

        if type(yspan) != tuple:
            yspan = (-1*yscale,yscale)
        
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
            
            #print(np.min(lumos) - np.max(homos),lumosk[np.argmin(lumos)]-homosk[np.argmax(homos)])
            
            # check for dirac points
            if np.min(lumos) - np.max(homos) <= 10*sigma and lumosk[np.argmin(lumos)]-homosk[np.argmax(homos)] <= sigma:
                plt.scatter(homosk[np.argmax(homos)],np.max(homos), marker = 'o', color = 'r')
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
                if dirEg != []:
                    x, y = dirEgk[np.argmin(dirEg)], dirEgE[np.argmin(dirEg)]
                    #print(x, y)
                    plt.scatter(x, y, marker='o', color='r')
                    
                    if size == None or size == "small":
                        plt.annotate(f"({x[0]:.4f}, {y[0]:.4f})", xy=[x[0] + 0.01, y[0] + 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                        plt.annotate(f"({x[1]:.4f}, {y[1]:.4f})", xy=[x[1] + 0.01, y[1] - 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                    elif size == "big":
                        plt.annotate(f"({x[0]:.4f}, {y[0]:.4f})", xy=[x[0] + 0.01, y[0] + 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale],fontsize=15)
                        plt.annotate(f"({x[1]:.4f}, {y[1]:.4f})", xy=[x[1] + 0.01, y[1] - 0.03 * (plt.ylim()[1] - plt.ylim()[0])/yscale],fontsize=15)
                    print(f'The direct bandgap, $E_g$ = {np.min(dirEg):.3f} eV')
                else:
                    indirEgk = [homosk[np.argmax(homos)], lumosk[np.argmin(lumos)]]
                    indirEgE = [np.max(homos), np.min(lumos)]
                    plt.scatter(indirEgk, indirEgE, marker='o', color='r')
                    
                    if size is None or size == "small":
                        plt.annotate(f"({indirEgk[0]:.4f},{indirEgE[0]:.4f})", xy=[indirEgk[0] + 0.01, indirEgE[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                        plt.annotate(f"({indirEgk[1]:.4f},{indirEgE[1]:.4f})", xy=[indirEgk[1] + 0.01, indirEgE[1] - 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale])
                    elif size == "big":
                        plt.annotate(f"({indirEgk[0]:.4f},{indirEgE[0]:.4f})", xy=[indirEgk[0] + 0.01, indirEgE[0] + 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale],fontsize=15)
                        plt.annotate(f"({indirEgk[1]:.4f},{indirEgE[1]:.4f})", xy=[indirEgk[1] + 0.01, indirEgE[1] - 0.3 * (plt.ylim()[1] - plt.ylim()[0])/yscale],fontsize=15)
                    print(f'The indirect bandgap, $E_g$ = {indirEgE[1] - indirEgE[0]:.3f} eV')

    if saveformat is not None:
        import os
        if not os.path.exists('00_Plots'):
            os.makedirs('00_Plots')
        plt.savefig(f'00_Plots/{filename}.{saveformat}', bbox_inches='tight')
        print(f"Plot saved as {filename}.{saveformat} in 00_Plots")
