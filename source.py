# import packages
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
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

# colour to integer mapping, try color=css[i] or for local functions, colour=css[colour]
css= ['red', 'Blue', 'green',
    'Brown', 'CadetBlue', 'Chocolate', 'Coral',
    'CornflowerBlue', 'Cornsilk', 'Crimson', 'Cyan', 'DarkBlue', 'DarkCyan',
    'DarkGoldenRod', 'DarkGray', 'DarkGreen', 'DarkKhaki', 'DarkMagenta',
    'DarkOliveGreen', 'DarkOrange', 'DarkOrchid', 'DarkRed', 'DarkSalmon',
    'DarkSeaGreen', 'DarkSlateBlue', 'DarkSlateGray', 'DarkTurquoise',
    'DarkViolet', 'DeepPink', 'DeepSkyBlue', 'DimGray', 'DodgerBlue',
    'FireBrick', 'FloralWhite', 'ForestGreen', 'Fuchsia', 'Gainsboro',
    'GhostWhite', 'Gold', 'GoldenRod', 'Gray', 'Green', 'GreenYellow',
    'HoneyDew', 'HotPink', 'IndianRed', 'Indigo', 'Ivory', 'Khaki',
    'Lavender', 'LavenderBlush', 'LawnGreen', 'LemonChiffon', 'LightBlue',
    'LightCoral', 'LightCyan', 'LightGoldenRodYellow', 'LightGreen',
    'LightGray', 'LightPink', 'LightSalmon', 'LightSeaGreen', 'LightSkyBlue',
    'LightSlateGray', 'LightSteelBlue', 'LightYellow', 'Lime', 'LimeGreen',
    'Linen', 'Magenta', 'Maroon', 'MediumAquaMarine', 'MediumBlue',
    'MediumOrchid', 'MediumPurple', 'MediumSeaGreen', 'MediumSlateBlue',
    'MediumSpringGreen', 'MediumTurquoise', 'MediumVioletRed', 'MidnightBlue',
    'MintCream', 'MistyRose', 'Moccasin', 'NavajoWhite', 'Navy', 'OldLace',
    'Olive', 'OliveDrab', 'Orange', 'OrangeRed', 'Orchid', 'PaleGoldenRod',
    'PaleGreen', 'PaleTurquoise', 'PaleVioletRed', 'PapayaWhip', 'PeachPuff',
    'Peru', 'Pink', 'Plum', 'PowderBlue', 'Purple', 'RebeccaPurple', 'Red',
    'RosyBrown', 'RoyalBlue', 'SaddleBrown', 'Salmon', 'SandyBrown', 'SeaGreen',
    'SeaShell', 'Sienna', 'Silver', 'SkyBlue', 'SlateBlue', 'SlateGray', 'Snow',
    'SpringGreen', 'SteelBlue', 'Tan', 'Teal', 'Thistle', 'Tomato', 'Turquoise',
    'Violet', 'YellowGreen'
]



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


import pandas as pd
import datetime
from pathlib import Path


def saveexcel(x, filename=None, header="auto"):
    """
    Export a pandas DataFrame to CSV with optional
    automatic header detection.

    Parameters
    ----------
    x : pandas.DataFrame
        DataFrame to export.

    filename : str or Path, optional
        Output filename. Automatically generated if None.

    header : "auto", True, or False
        "auto" : Automatically detect and correct headers.
        True   : Preserve existing DataFrame column names.
        False  : Export data without column names.

    Returns
    -------
    pathlib.Path
        Path to the exported CSV file.
    """

    # check inputs
    if not isinstance(x, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")

    if header not in ("auto", True, False):
        raise ValueError(
            "header must be 'auto', True, or False."
        )

    # generate filename
    if filename is None:
        timestamp = datetime.datetime.now().strftime(
            "%Y%m%d_%H%M%S"
        )
        filename = f"export_{timestamp}.csv"

    filename = Path(filename)

    if filename.suffix.lower() != ".csv":
        filename = filename.with_suffix(".csv")

   # protect original dataframe
    df = x.copy()

    # detect pandas structure
    has_default_columns = (
        isinstance(df.columns, pd.RangeIndex)
        and df.columns.equals(
            pd.RangeIndex(df.shape[1])
        )
    )

    #determine header behaviour
    write_header = True
    if header == "auto":
        if has_default_columns:
            # suppress default integer headers
            write_header = False

            # check first row for header
            if len(df) >= 2:
                row_0_numeric = pd.to_numeric(df.iloc[0], errors="coerce")
                row_1_numeric = pd.to_numeric(df.iloc[1], errors="coerce")
                row_0_is_text = (row_0_numeric.isna().all())
                row_1_has_numbers = (row_1_numeric.notna().any())

                if (row_0_is_text and row_1_has_numbers):

                    # promote first row to column name
                    df.columns = df.iloc[0].astype(str)
                    # remove promoted header
                    df = df.iloc[1:].reset_index(drop=True)
                    write_header = True
    else:
        # user input
        write_header = header
        
    df.to_csv(
        filename,
        index=False,
        header=write_header
    )

    print(f"Successfully exported: {filename}")

    return filename


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
def curvefit(xdata, ydata, fit_type, initial_guess=None, force_coord=None, colour=None, fitname=None, lgd=False, smooth=None):
    xdata = np.asarray(xdata, dtype=float)
    ydata = np.asarray(ydata, dtype=float)
    def linear(x, a, b):
        return a * x + b
    def lin0(x, a):
        return a * x
    def norm_linear(x, a):
        return a * x + 1
    def exponential(x, a, b, c):
        return a * np.exp(c * x) + b
    def norm_exponential(x, a):
        return np.exp(a * x)
    def zeroed_norm_exponential(x, a):
        return np.exp(a * x) - 1
    def log(x, a, b):
        return a * np.log(x) + b
    def sinusoidal(x, a, b, c, d):
        return a * np.sin(b * x) + c * np.cos(d * x)
    def gaussian(x, a, b, c, d):
        return a * np.exp(-((x - b) ** 2)/(2 * c ** 2)) + d
    def polynomial(x, *coeffs):
        return np.polyval(coeffs, x)
    #dictionary
    fit_functions = {
        'lin': (linear, 'y = {:.2f}x + {:.2f}'),
        'lin0': (lin0, 'y = {:.2f}x'),
        'linnorm': (norm_linear, 'y = {:.2f}x + 1'),
        'exp': (exponential, 'y = {:.2f}exp({:.2f}x) + {:.2f}'),
        'expnorm': (norm_exponential, 'y = exp({:.2f}x)'),
        'expnorm0': (zeroed_norm_exponential, 'y = exp({:.2f}x) - 1'),
        'log': (log, 'y = {:.2f} log(x) + {:.2f}'),
        'sinusoidal': (sinusoidal, 'y = {:.2f}sin({:.2f}x) + {:.2f}cos({:.2f}x)'),
        'gaussian': (gaussian, 'y = {:.2f} * exp(-(x - {:.2f})^2 / (2 * {:.2f}^2)) + {:.2f}')
    }
    if fit_type == 'smooth': #handle smooth
        sort_idx = np.argsort(xdata) #sort ascending
        xdata_sorted = xdata[sort_idx]
        ydata_sorted = ydata[sort_idx]
        spline = UnivariateSpline(xdata_sorted, ydata_sorted, s=smooth)
        x_smooth = np.linspace(xdata_sorted.min(), xdata_sorted.max(), 500)
        y_smooth = spline(x_smooth)
        y_fit = spline(xdata_sorted)
        SS_res = np.sum((ydata_sorted - y_fit)**2)
        SS_tot = np.sum((ydata_sorted - np.mean(ydata_sorted))**2)
        R_squared = 1 - (SS_res/SS_tot)
        print('Spline smoothing with smoothing factor:', smooth)
        print(f'R²: {R_squared:.4f}')
        if colour==None:
            colour=45
        plt.plot(x_smooth, y_smooth, '-', color=css[colour], label=f'{fitname}')
        plt.legend()

        return spline

    elif fit_type == 'polyforce': #to force coordinates
        if force_coord is None: 
            raise ValueError('For polyforce, you must supply force_coord=[force_x, force_y].') #error handling
        force_x = np.asarray(force_coord[0], dtype=float)
        force_y = np.asarray(force_coord[1], dtype=float)
        if len(force_x) != len(force_y):
            raise ValueError('force_x and force_y must have the same number of elements.')
        order = int(input('Insert order of polynomial fit (polyforce): '))
        E = np.vander(force_x, N=order+1, increasing=False) #build constraint matrix w vandermonde (descending coeff)
        f_vec = force_y
        c_p, residuals, rank, s = np.linalg.lstsq(E, f_vec, rcond=None) # solve for E*c=f_vec
        U, S, Vh = np.linalg.svd(E) #complete null space of E via SVD, null space dim = order + 1 - rank
        tol = np.max(E.shape) * np.amax(S) * np.finfo(S.dtype).eps
        ns_dim = (order+1) - rank
        if ns_dim > 0:
            Z = Vh[-ns_dim:].T
        else:
            Z = np.zeros((order+1, 0))
        V_data = np.vander(xdata, N=order+1, increasing=False)
        if ns_dim > 0:
            A = V_data @ Z
            b = ydata - V_data @ c_p
            d, residuals_d, rank_d, s_d = np.linalg.lstsq(A, b, rcond=None)
            coeff = c_p + Z @ d
        else:
            coeff = c_p
        def poly_func(x):
            return np.polyval(coeff, x)
        fitted_func = poly_func
        eqn_terms = []
        for i, c in enumerate(coeff):
            power = order - i
            if power == 0:
                term = f"{c:.2f}"
            elif power == 1:
                term = f"{c:.2f}x"
            else:
                term = f"{c:.2f}x^{power}"
            eqn_terms.append(term)
        eqn = "y = " + " + ".join(eqn_terms)
        y_fitted = fitted_func(xdata)
        SS_res = np.sum((ydata - y_fitted)**2)
        SS_tot = np.sum((ydata - np.mean(ydata))**2)
        R_squared = 1 - (SS_res/SS_tot)
        covar = None

    elif fit_type == 'poly':
        o = int(input('Insert order of polynomial fit: '))
        coeff = np.polyfit(xdata, ydata, o)
        fitted_func = polynomial
        eqn = 'y = ' + ' + '.join(f'{p:.2f}x^{i}' for i, p in enumerate(coeff[::-1]))
        y_fitted = np.polyval(coeff, xdata)
        SS_res = np.sum((ydata - y_fitted)**2)
        SS_tot = np.sum((ydata - np.mean(ydata))**2)
        R_squared = 1 - (SS_res/SS_tot)

    elif fit_type in fit_functions:
        func, eqn_template = fit_functions[fit_type]
        coeff, covar = curve_fit(func, xdata, ydata, p0=initial_guess, maxfev=10000)
        fitted_func = lambda x: func(x, *coeff)
        eqn = eqn_template.format(*coeff)
        y_fitted = fitted_func(xdata)
        SS_res = np.sum((ydata - y_fitted)**2)
        SS_tot = np.sum((ydata - np.mean(ydata))**2)
        R_squared = 1 - (SS_res/SS_tot)
        if '0' in fit_type:
            xdata = np.append([0], xdata)
            y_fitted = np.append([0], y_fitted)
        if fit_type == 'gaussian':
            print(f'------------------XRD Maxima for μ = {coeff[1]:.2f}°, σ = {coeff[2]:.2f}°------------------')
            print('Fitted parameters:', coeff)
            print('Fitted equation:', eqn)
    else:
        raise ValueError('Invalid fit_type. Expected one of: ' + ', '.join(list(fit_functions.keys()) + ['poly', 'polyforce', 'smooth']))
    print('Fitted parameters:', coeff)
    print('Fitted equation:', eqn)
    print(f'R²: {R_squared:.4f}')
    if fit_type not in ['poly', 'polyforce'] and fit_type in fit_functions:
        print('Covariance Matrix:', covar)
    if colour is None:
        colour = 'g'
    elif isinstance(colour, int):
        colour = css[colour % len(css)]

    #plotting
    if lgd:
        if fitname is None:
            fitname = 'Fitted Curve'
        plt.plot(xdata, y_fitted, '-', color=str(colour), label=str(fitname))
        plt.legend()
    else:
        plt.plot(xdata, y_fitted, '-', color=str(colour))
        plt.legend()

    if fit_type == 'poly':
        return coeff, R_squared
    elif fit_type in ['polyforce'] or fit_type in fit_functions:
        return coeff, covar, R_squared
    else:
        return coeff


def curvefitold(xdata, ydata, fit_type, initial_guess=None):
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

    def exponential(x, a, b, c):
        return a * np.exp(c * x) + b

    def norm_exponential(x, a):
        return np.exp(a * x)

    def zeroed_norm_exponential(x, a):
        return np.exp(a * x) - 1

    def log(x,a,b):
        return a * np.log(x) + b
        
    def sinusoidal(x, a, b, c, d):
        return a * np.sin(b * x) + c * np.cos(d * x)

    def gaussian(x, a, b, c, d):
        return a * np.exp(-(x - b) ** 2 / (2 * c ** 2)) + d

    def polynomial(x, *coeffs):
        return np.polyval(coeffs, x)

    # dictionary mapping fit types to functions
    fit_functions = {
        'lin': (linear, 'y = {:.2f}x + {:.2f}'),
        'lin0': (lin0, 'y={:.2f}x'),
        'linnorm': (norm_linear, 'y = {:.2f}x + 1'),
        'exp': (exponential, 'y = {:.2f}exp({:.2f}x) + {:.2f}'),
        'expnorm': (norm_exponential, 'y = exp({:.2f}x)'),
        'expnorm0': (zeroed_norm_exponential, 'y = exp({:.2f}x) - 1'),
        'log': (log, 'y = {:.2f} log(x) + {:.2f}'),
        'sinusoidal': (sinusoidal, 'y = {:.2f}sin({:.2f}x) + {:.2f}cos({:.2f}x)'),
        'gaussian': (gaussian, 'y = {:.2f} * exp(-(x - {:.2f})^2 / (2 * {:.2f}^2)) + {:.2f}')
    }

    # handle polynomial separately
    if fit_type == 'poly':
        o = int(input('Insert order of polynomial fit: '))
        coeff = np.polyfit(xdata, ydata, o)
        fitted_func = polynomial
        eqn = 'y = ' + ' + '.join(f'{p:.2f}x^{i}' for i, p in enumerate(coeff[::-1]))

        # Compute R^2 for polynomial
        y_fitted = np.polyval(coeff, xdata)
        SS_res = np.sum((ydata - y_fitted) ** 2)
        SS_tot = np.sum((ydata - np.mean(ydata)) ** 2)
        R_squared = 1 - (SS_res / SS_tot)

    elif fit_type in fit_functions:
        func, eqn_template = fit_functions[fit_type]
        coeff, covar = curve_fit(func, xdata, ydata, p0=initial_guess, maxfev=10000)
        fitted_func = func
        eqn = eqn_template.format(*coeff)

        # compute R^2 for all other fit types
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

    # Print results
    print('Fitted parameters:', coeff)
    print('Fitted equation:', eqn)
    print(f'R^2: {R_squared:.4f}')
    if fit_type != 'poly':
        print('Covariance Matrix:', covar)

    # Plotting
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

def debye(mu,dmu,sigma,dsigma,em,K=None,scale=None):
    # first defining lambda
    if em == 'CuKa':
        Lambda = 1.5418e-10 # Cu K\alpha
    elif em == None:
        raise ValueError('Insert blackbody emitter!')

    if K == 'None':
        K = 0.89
    elif K == 'sphere':
        K = 0.98
    
    # calculating the crystallite size
    if isinstance(mu and sigma and dmu and dsigma, (int,float)):
        dtheta,dbeta = (dmu/2)*np.pi/180,(dsigma*np.sqrt(8*np.log(2)))*np.pi/180
        beta = (np.pi/180) * np.sqrt(8*np.log(2)) * sigma
        theta = (np.pi/180)*(mu/2)
        D = (0.98 * Lambda) / beta * np.cos(theta)
        dD = np.sqrt(D**2 * ((dbeta/beta)**2 + (dtheta*np.tan(theta))**2))
        #return [D,dD]
        print(f'----------Calculated Crystallite Size of {NP_name}----------')
        print(f'D = {D:} \u00B1 {dD}')

    elif isinstance(mu and sigma and dmu and dsigma, list):
        dtheta_intm,dbeta_intm = [dmu[i]/2 * np.pi/180 for i in range(len(dmu))],[dsigma[i]*np.sqrt(8*np.log(2))*np.pi/180 for i in range(len(dsigma))]
        beta_intm = [(np.pi/180) * np.sqrt(8*np.log(2)) * sigma[i] for i in range(len(sigma))]
        theta_intm = [(np.pi/180) * mu[i] / 2 for i in range(len(mu))]
        D_intm = [(0.98 * Lambda) / beta_intm[i] * np.cos(theta_intm[i]) for i in range(len(beta_intm))] # or = len(theta_intm)
        #dD = [np.sqrt(D_intm**2 * (dbeta**2 + (dtheta*np.sin(i))**2)) for i in theta]
        dD = [np.sqrt(D_intm[i]**2 * ((dbeta_intm[i]/beta_intm[i])**2 + (dtheta_intm[i]*np.sin(theta_intm[i]))**2)) for i in range(len(D_intm))] # also = len(dbeta_intm,dtheta_intm,theta_intm)

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

        if scale in (None,'nm','nano'):
            print(f'----------Calculated Crystallite Size of {NP_name}----------')
            for i in range(len(D)):
                print(f'Signal {i+1}: D = {D[i][0]*(1e9):.3f} \u00B1 {D[i][1]*(1e9):.3f} nm')
                print(f'\u03B8 = {theta[i][0]:.4f} \u00B1 {theta[i][1]:.4f} rad')
                print(f'\u03B2 = {beta[i][0]:.4f} \u00B1 {beta[i][1]:.4f} rad')
        
            Dbar,dDbar=0,0
            for i in range(len(D)):
                Dbar += D[i][0]
                dDbar += D[i][1]  
            Dbar,dDbar=Dbar/len(D),dDbar/len(D)
            print(f'Average Crystallite Size = {Dbar*1e9:.3f} \u00B1 {dDbar*1e9:.3f} nm')
    
        elif scale in ('AA','Ang','ang','Angstrom','Angstrom'):
            print(f'----------Calculated Crystallite Size of {NP_name}----------')
            for i in range(len(D)):
                print(f'Signal {i+1}: D = {D[i][0]*(1e10):.3f} \u00B1 {D[i][1]*(1e10):.3f} Å')
                print(f'\u03B8 = {theta[i][0]:.4f} \u00B1 {theta[i][1]:.4f} rad')
                print(f'\u03B2 = {beta[i][0]:.4f} \u00B1 {beta[i][1]:.4f} rad')
        
            Dbar,dDbar=0,0
            for i in range(len(D)):
                Dbar += D[i][0]
                dDbar += D[i][1]  
            Dbar,dDbar=Dbar/len(D),dDbar/len(D)
            print('-----------------------------------------------------------------')
            print(f'---------Average Crystallite Size = {Dbar*1e10:.3f} \u00B1 {dDbar*1e10:.3f} Å-----------')
            print('-----------------------------------------------------------------')

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

def uvvisplot(filename, fit = None,spectra = None, name = None, colour = None):
    try:
        uvvisdata = loadexcel(filename)
        if uvvisdata is None:
            raise FileNotFoundError
    except FileNotFoundError:
        uvvisdata = copypaste()
        saveexcel(uvvisdata, f'{filename}')
        print(f'"{filename}.csv" has been created. Process file in Excel to an Nx2 column and rerun the function.')
        return
    
    uvvis = plotarray(uvvisdata)

    if colour == None:
        colour = 'b'
    
    if plt.rcParams['figure.dpi'] == 1000:
        size = 15 
    else:
        size = None
    
    if spectra == None:
        i = 1
    else:
        i = spectra
    plt.plot(uvvis[0],uvvis[i], '-', color = f"{colour}", label = f"{name}")
    #plt.title(f"UV-vis Spectra of {name}")
    plt.xlabel('Wavelength, $\lambda$ (nm)',fontsize=size)
    plt.ylabel('Absorbance (a.u.)',fontsize=size)
    
    if fit == True:
        curvefit(uvvis[0], uvvis[1], 'poly')
    else:
        return
    plt.legend(fontsize=size)


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
            all_bonds = np.concatenate((explicit_bonds, np.array(inferred_bond_lengths)))

            # Perform HOMA calculation
            if 's1' in filename or 't1' in filename:
                alpha, Ropt, n = 950.74, 1.437, len(all_bonds)
                #print("All bond distances (explicit + inferred):", all_bonds)
                HOMERval = 1 - (alpha / n) * np.sum((all_bonds - Ropt) ** 2)
                print(f"HOMER value for excited-state {n}-\u03C0 conjugated {filename} is {HOMERval:.3f}")
            else:
                alpha, Ropt, n = 257.70, 1.388, len(all_bonds)
                #print("All bond distances (explicit + inferred):", all_bonds)
                HOMAval = 1 - (alpha / n) * np.sum((all_bonds - Ropt) ** 2)
                print(f"HOMA value for groundstate {n}-\u03C0 conjugated {filename} is {HOMAval:.3f}")

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
            alpha, Ropt, n = 257.70, 1.388, len(bondlist)
            HOMAval = 1 - (alpha / n) * np.sum((bondlist - Ropt) ** 2)
            print(f"HOMA value for groundstate {n}-\u03C0 conjugated {filename} is {HOMAval:.3f}")

    else:
        print(f'This mode of analysis is not supported as of {os.times()}.')

    print('The parameterised values used in calculation were obtained with gratitude from Enrique and Bo, Phys. Chem. Chem. Phys., 2023, 25, 16763–16771')
