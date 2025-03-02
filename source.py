# import packages
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

# Set all text elements to use 'Times New Roman' as the font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# Set math text to use 'Times New Roman'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['mathtext.it'] = 'serif:italic'
plt.rcParams['mathtext.bf'] = 'serif:bold'
plt.rcParams['mathtext.fontset'] = 'custom'


# high res
plt.rcParams['figure.dpi'] = 100 #default 100

# data loading functions 

def loadexcel(filename = None):
    if filename == None:
        filename = input("Insert .csv filename here: ")
        if filename.endswith(".csv"):
            df = pd.read_csv(f"{filename}")
        else:
            df = pd.read_csv(f"{filename}.csv")
    else:
        if filename.endswith(".csv"):
            df = pd.read_csv(f"{filename}")
        else:
            df = pd.read_csv(f"{filename}.csv")
        return df

def saveexcel(x, filename = None):
    if filename == None:
        filename = input("Insert .csv filename here: ")
        if filename.endswith(".csv"):
            intm = x.to_csv(f'{filename}', index = False)
        else:
            intm = x.to_csv(f'{filename}.csv', index = False)
    else:
        if filename.endswith(".csv"):
            intm = x.to_csv(f'{filename}', index = False)
        else:
            intm = x.to_csv(f'{filename}.csv', index = False)
    return intm

def copypaste():
    sep = [',', ';', '|', '/','\\', r'\s+']
    for i in sep:
        try:
            df = pd.read_clipboard(sep=i)
            if df is not None:
                break
        except:
            pass
    return df

def plotarray(filename):
    intm = np.array(filename)
    return [intm[:, i] for i in range(intm.shape[1])]

#advanced functions

def saveloadplot(filename):
    try:
        data = loadexcel(filename)
        if data is None:
            raise FileNotFoundError
    except FileNotFoundError:
        data = copypaste()

        if data is None or data.empty:
            print("Error: Clipboard empty. Aborting operation.")
            return
        
        saveexcel(data, filename)
        print(f'"{filename}.csv" has been created. Process the file if needed, then rerun the function.')

    return plotarray(data)

## advanced functions

def curvefit(xdata, ydata, fit_type, initial_guess=None):
    # Ensure xdata and ydata are numpy arrays
    xdata = np.asarray(xdata, dtype=float)
    ydata = np.asarray(ydata, dtype=float)

    # Define fitting functions
    def linear(x, a, b):
        return a * x + b
    def lin0(x,a):
        return a*x
    def norm_linear(x, a):
        return a * x + 1
    def zeroed_norm_linear(x, a):
        return a * x
    def exponential(x, a, b, c):
        return a * np.exp(c * x) + b
    def norm_exponential(x, a):
        return np.exp(a * x)
    def zeroed_norm_exponential(x, a):
        return np.exp(a * x) - 1
    def sinusoidal(x, a, b, c, d):
        return a * np.sin(b * x) + c * np.cos(d * x)
    def gaussian(x, a, b, c, d):
        return a * np.exp(-(x - b) ** 2 / (2 * c ** 2)) + d
    def polynomial(x, *coeffs):
        return np.polyval(coeffs, x)
    # dictionary mapping
    fit_functions = {'lin': (linear, 'y = {:.2f}x + {:.2f}'), 'lin0': (lin0, 'y={:.2f}x'), 'linnorm': (norm_linear, 'y = {:.2f}x + 1'), '0linnorm': (zeroed_norm_linear, 'y = {:.2f}x'), 'exp': (exponential, 'y = {:.2f}exp({:.2f}x) + {:.2f}'), 'expnorm': (norm_exponential, 'y = exp({:.2f}x)'), '0expnorm': (zeroed_norm_exponential, 'y = exp({:.2f}x) - 1'), 'sinusoidal': (sinusoidal, 'y = {:.2f}sin({:.2f}x) + {:.2f}cos({:.2f}x)'), 'gaussian': (gaussian, 'y = {:.2f} * exp(-(x - {:.2f})^2 / (2 * {:.2f}^2)) + {:.2f}')}

    # polynomial special treatment
    if fit_type == 'poly':
        o = int(input('Please insert order of polynomial fit: '))
        coeff = np.polyfit(xdata, ydata, o)
        fitted_func = polynomial
        eqn = 'y = ' + ' + '.join(f'{p:.2f}x^{i}' for i, p in enumerate(coeff[::-1]))
        # R^2
        y_fitted = np.polyval(coeff, xdata)
        SS_res = np.sum((ydata - y_fitted) ** 2)
        SS_tot = np.sum((ydata - np.mean(ydata)) ** 2)
        R_squared = 1 - (SS_res / SS_tot)
        
    elif fit_type in fit_functions:
        func, eqn_template = fit_functions[fit_type]
        coeff, covar = curve_fit(func, xdata, ydata, p0=initial_guess, maxfev=10000)
        fitted_func = func
        eqn = eqn_template.format(*coeff)
        # R^2 in general
        y_fitted = fitted_func(xdata, *coeff)
        SS_res = np.sum((ydata - y_fitted) ** 2)
        SS_tot = np.sum((ydata - np.mean(ydata)) ** 2)
        R_squared = 1 - (SS_res / SS_tot)
        
        if '0' in fit_type:
            xdata = np.append([0],xdata)
            y_fitted = np.append([0],y_fitted)
        if fit_type == 'gaussian':  # XRD signal fitting
            print(f'------------------XRD Maxima for \u03BC = {coeff[1]:.2f}\u00B0, \u03C3 = {coeff[2]:.2f}\u00B0------------------')
            print('Fitted parameters:', coeff)
            print('Fitted equation:', eqn)
    else:
        raise ValueError('Invalid fit_type. Expected one of: ' + ', '.join(fit_functions.keys()) + ', or "poly".')
    
    print('Fitted parameters:', coeff)
    print('Fitted equation:', eqn)
    print(f'R^2: {R_squared:.4f}')
    if fit_type != 'poly':
        print('Covariance Matrix:', covar)

    plt.plot(xdata, y_fitted, 'g-', label='Fitted Curve')
    plt.legend()


def xrdplot(debye_fit = None, fitspread = None, NP_name = None, filename = None, colour = None):
    if debye_fit == None and NP_name != None: # import
        xrddata = copypaste()
        if xrddata is not None:
            xrd = plotarray(xrddata)  
            print(xrddata)
            plt.plot(xrd[0], xrd[1],'b-',label=f'{NP_name}')
            plt.title(f'XRD plot of {NP_name}')
            plt.xlabel('Detector Angle, 2$\\theta$ ($^o$)')
            plt.ylabel('Counts (a.u.)')
            plt.legend()
            saveexcel(xrddata, f'{NP_name}_xrd')
        else:
            return None
    
    elif debye_fit != None: # standard display and/or gaussian fit
        xrddata = loadexcel(filename)
        xrd = plotarray(xrddata)
        
        if colour == None:
            colour = 'r'
        if isinstance(debye_fit,(int,float)) and debye_fit != 0:
            # first isolating the n maximas
            n = int(debye_fit)
            from scipy.signal import argrelextrema
            comparator = np.greater
            indices = argrelextrema(np.array(xrd[1]), comparator)[0] # this yields a list of extrema indices from left to right
            intm = [xrd[1][i] for i in indices] # get the intensities of the maxima from left to right
            extremaindices = indices[np.argsort(intm)[-n:]][::-1]
            # np.argsort sorts min to max, [-n:] splices the last n elements, [::-1] reverses the array
            
            # nth plot
            for i in range(len(extremaindices)):
                # now, defining range to curvefit()
                xrdzoom=[[],[]] 
                for j in range(extremaindices[i]-100*fitspread,extremaindices[i]+100*fitspread):
                    xrdzoom[0].append(xrd[0][j])
                    xrdzoom[1].append(xrd[1][j])
                    
                plt.figure(figsize=(6.4,4.8))
                #plt.subplot(1,2,2)
                tallest=[[],[]]
                tallest[0],tallest[1] = xrd[0][extremaindices[i]], xrd[1][extremaindices[i]]
                plt.plot(xrdzoom[0],xrdzoom[1],'b',label=f'{NP_name}_maxima#{i+1}_Zoomed')
                plt.xlabel('Detector Angle, 2$\\theta$ ($\circ$)')
                plt.ylabel('Intensity (a.u.)')
                plt.text(tallest[0], tallest[1], f'({tallest[0]:.2f}, {tallest[1]:.0f})', ha='right', va='bottom')
                a_init, b_init, c_init, d_init = max(xrdzoom[1]), xrdzoom[0][np.argmax(xrdzoom[1])], np.std(xrdzoom[0]), min(xrdzoom[1])
                p0 = [a_init, b_init, c_init, d_init] # to enable gaussian fit
                curvefit(xrdzoom[0], xrdzoom[1], 'gaussian',initial_guess = p0)
                plt.tight_layout()
                
            #plt.figure(figsize=(12.8,4.8))
            plt.figure(figsize=(6.4,4.8))
            # main/first plot
            plt.plot(xrd[0],xrd[1], '-', color = f"{colour}", label = f"{NP_name}")
            plt.xlabel('Detector Angle, 2$\\theta$ ($^o$)')
            plt.ylabel('Counts (a.u.)')
            plt.legend()
            plt.text(tallest[0], tallest[1], f'({tallest[0]:.2f}, {tallest[1]:.0f})', ha='right', va='bottom')
            plt.tight_layout()

        elif debye_fit == 0: # standard plot
            #plt.figure(figsize=(6.4,4.8)) # standard size in plt
            plt.plot(xrd[0],xrd[1], '-', color = f"{colour}", label = f"{NP_name}")
            plt.xlabel('Detector Angle, 2$\\theta$ ($^o$)')
            plt.ylabel('Counts (a.u.)')
            plt.legend()
            tallest=[[],[]]
            tallest[0],tallest[1] = xrd[0][np.argmax(xrd[1])], xrd[1][np.argmax(xrd[1])]
            plt.text(tallest[0], tallest[1], f'({tallest[0]:.2f}, {tallest[1]:.0f})', ha='right', va='bottom')
            plt.tight_layout()

    elif debye_fit == None and NP_name == None: # error
        print('Either "debye_fit" or "NP_name" argument is required!')
        return None

def debye(mu,dmu,sigma,dsigma,em):
    # first defining lambda
    if em == 'CuKa':
        Lambda = 1.541862 * 10**(-10) # Cu K\alpha
    
    # calculating the crystallite size
    if isinstance(mu and sigma and dmu and dsigma, (int,float)):
        dtheta,dbeta = (dmu/2)*np.pi/180,(dsigma*np.sqrt(8*np.log(2)))*np.pi/180
        beta = (np.pi/180) * np.sqrt(8*np.log(2)) * sigma
        theta = (np.pi/180) * 1/2 * mu
        D = (0.98 * Lambda) / beta * np.cos(theta)
        dD = np.sqrt(D**2 * (dbeta**2 + (dtheta*np.sin(theta))**2))
        #return [D,dD]
        print(f'----------Calculated Crystallite Size of {NP_name}----------')
        print(f'D = {D:} \u00B1 {dD}')

    elif isinstance(mu and sigma and dmu and dsigma, list):
        dtheta_intm,dbeta_intm = [i/2 * np.pi/180 for i in dmu],[i*np.sqrt(8*np.log(2))*np.pi/180 for i in dsigma]
        beta_intm = [(np.pi/180) * np.sqrt(8*np.log(2)) * i for i in sigma]
        theta_intm = [(np.pi/180) * i / 2 for i in mu]
        D_intm = [(0.98 * Lambda) / i * np.cos(j) for i,j in zip(beta_intm,theta_intm)]
        #dD = [np.sqrt(D_intm**2 * (dbeta**2 + (dtheta*np.sin(i))**2)) for i in theta]
        dD = [np.sqrt(i**2 * (j**2 + (k*np.sin(l))**2)) for i,j,k,l in zip(D_intm,dbeta_intm,dtheta_intm,theta_intm)]

        D,beta,theta=[],[],[]
        for i in range(len(D_intm) and len(dD)):
            D.append([])
            D[i].append(D_intm[i])
            D[i].append(dD[i])

            beta.append([])
            beta[i].append(beta_intm[i])
            beta[i].append(dbeta_intm[i])

            theta.append([])
            theta[i].append(theta_intm[i])
            theta[i].append(dtheta_intm[i])

        print(f'----------Calculated Crystallite Size of {NP_name}----------')
        for i in range(len(D)):
            print(f'Signal {i+1}: D = {D[i][0]*(10**9):.3f} \u00B1 {D[i][1]*(10**9):.3f} nm')
            print(f'\u03B8 = {theta[i][0]} \u00B1 {theta[i][1]} rad')
            print(f'\u03B2 = {beta[i][0]} \u00B1 {beta[i][1]} rad')

def ftirplot(filename):
    ftirdata = loadexcel(filename)
    ftir = plotarray(ftirdata)
    plt.plot(ftir[0],ftir[1], '-', color = f"{colour}", label = f"{NP_name}")
    #plt.title(f"FTIR Plot of {NP_name}")
    plt.xlabel('Wavenumber, (cm$^{-1}$)')
    plt.ylabel('Transmittance (%)')
    plt.gca().invert_xaxis()
    plt.ylim(0)
    #plt.legend()
    #plt.savefig(f'{filename}.png')

def uvvisplot(filename, fit = None,spectra = None):
    uvvisdata = loadexcel(filename)
    uvvis = plotarray(uvvisdata)
    if spectra == None:
        n = 1
    else:
        n = spectra
    plt.plot(uvvis[0],uvvis[n], '-', color = f"{colour}", label = f"{NP_name}")
    plt.title(f"UV-vis Spectra of {NP_name}")
    plt.xlabel('Wavelength, (nm)')
    plt.ylabel('Absorbance (a.u.)')
    if fit == True:
        curvefit(uvvis[0], uvvis[1], 'poly')
    else:
        return
    #plt.legend()
    #plt.savefig(f'{filename}.png')

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


def HOMA(folderpath=None, filename=None, analysis=None, correction=None):
    
    if analysis not in (99, "manual"):
        # Check if the folder exists
        if not os.path.exists(folderpath):
            print(f"The folder '{folderpath}' does not exist.")
            return

        filepath = os.path.join(folderpath, filename)

        try:
            ds = loadexcel(filepath)
            if ds is None:
                raise FileNotFoundError
        except FileNotFoundError:
            ds = copypaste()
            saveexcel(ds, filepath)
            print(f'"{filename}.csv" has been created in {filepath}. Process file in Excel and rerun the function.')
            return
        
        if analysis == 'local' or analysis == 0:
            # process all carbons with bond data
            carbonds = pd.DataFrame(ds[(ds['Symbol'] == 'C') & (ds['Highlight'] == 'Yes')])
            carbonds = carbonds[['Tag', 'X', 'Y', 'Z', 'Bond']].reset_index(drop=True)

            # identify carbons with NaN bond values
            nan_bonds = carbonds[carbonds['Bond'].isna()]
            inferred_bond_lengths = []

            # calculate bond length only for carbons with NaN bond values
            if correction == None:
                for index, row in nan_bonds.iterrows():
                    x1, y1, z1 = row[['X', 'Y', 'Z']]
                    for i, other_row in carbonds.iterrows():
                        if i != index:  # skip self-comparison
                            x2, y2, z2 = other_row[['X', 'Y', 'Z']]
                            distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

                            # include bond only if it falls within the bonding range
                            if 1.2 <= distance <= 1.6:
                                inferred_bond_lengths.append(distance)
            elif correction is not None:
                x1, y1, z1 = nan_bonds.iloc[0][['X', 'Y', 'Z']]
                connecting_carbon = carbonds[carbonds['Tag'] == correction]
                
                if connecting_carbon.empty:
                    raise ValueError(f"Connecting carbon with tag '{correction}' not found in carbonds.")

                x2, y2, z2 = connecting_carbon.iloc[0][['X', 'Y', 'Z']]
                distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
                #print(distance)
                inferred_bond_lengths.append(distance)
                
            # combine explicitly listed bonds with inferred ones
            explicit_bonds = carbonds['Bond'].dropna().to_numpy()
            all_bonds = np.unique(np.concatenate((explicit_bonds, np.array(inferred_bond_lengths))))

            # Perform HOMA calculation
            if 's1' in filename or 't1' in filename:
                alpha, Ropt, n = 950.74, 1.437, len(all_bonds)
                #print("All bond distances (explicit + inferred):", all_bonds)
                HOMERval = 1 - (alpha / n) * np.sum((all_bonds - Ropt) ** 2)
                print(f"HOMER value for excited-state {n}-\u03C0 conjugated {filename} is {HOMERval:.3f}")
            else:
                alpha, Ropt, n = 153.37, 1.392, len(all_bonds)
                #print("All bond distances (explicit + inferred):", all_bonds)
                HOMAval = 1 - (alpha / n) * np.sum((all_bonds - Ropt) ** 2)
                print(f"HOMAc value for groundstate {n}-\u03C0 conjugated {filename} is {HOMAval:.3f}")

        elif analysis == 'multiwfn' or analysis == 1:
            carbonds = pd.DataFrame(ds[(ds['Highlight'] == 'Yes') & (ds['Symbol'] == 'C')])
            carbonlist = carbonds['Tag'].to_numpy()
            output = ', '.join(map(str, carbonlist))
            print(output)

    elif analysis == 'manual' or analysis == 99:
        
        if correction is None:
            raise ValueError('Please provide bondlist in a list under argument "correction".')
        
        bondlist = np.array(correction, dtype=float)

        if 's1' in filename or 't1' in filename:
            alpha, Ropt, n = 950.74, 1.437, len(bondlist)
            HOMERval = 1 - (alpha / n) * np.sum((bondlist - Ropt) ** 2)
            print(f"HOMER value for excited-state {n}-\u03C0 conjugated {filename} is {HOMERval:.3f}")
        else:
            alpha, Ropt, n = 153.37, 1.392, len(bondlist)
            HOMAval = 1 - (alpha / n) * np.sum((bondlist - Ropt) ** 2)
            print(f"HOMAc value for groundstate {n}-\u03C0 conjugated {filename} is {HOMAval:.3f}")

    else:
        print(f'This mode of analysis is not supported as of {os.times()}.')

    print('The parameterised values used in calculation were obtained with gratitude from Enrique and Bo, Phys. Chem. Chem. Phys., 2023, 25, 16763–16771')

