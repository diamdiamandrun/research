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


# high res, switch to 1000
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
    intm = pd.DataFrame(filename).to_numpy() #intm[row/no. of datapts][col/no. of var.]
    x,rows,cols = [], len(intm), len(np.transpose(intm))
    for i in range(cols):
        x.append([])
        for j in range(rows):
            x[i].append(intm[j][i]) # to create a list of all values in each col
    return x

def curvefit(xdata, ydata, fit_type, initial_guess=None):
    # Ensure xdata and ydata are numpy arrays
    xdata = np.asarray(xdata, dtype=float)
    ydata = np.asarray(ydata, dtype=float)

    # Define fitting functions
    def linear(x, a, b):
        return a * x + b

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

    # Dictionary mapping fit types to functions
    fit_functions = {
        'lin': (linear, 'y = {:.2f}x + {:.2f}'),
        'linnorm': (norm_linear, 'y = {:.2f}x + 1'),
        '0linnorm': (zeroed_norm_linear, 'y = {:.2f}x'),
        'exp': (exponential, 'y = {:.2f}exp({:.2f}x) + {:.2f}'),
        'expnorm': (norm_exponential, 'y = exp({:.2f}x)'),
        '0expnorm': (zeroed_norm_exponential, 'y = exp({:.2f}x) - 1'),
        'sinusoidal': (sinusoidal, 'y = {:.2f}sin({:.2f}x) + {:.2f}cos({:.2f}x)'),
        'gaussian': (gaussian, 'y = {:.2f} * exp(-(x - {:.2f})^2 / (2 * {:.2f}^2)) + {:.2f}')
    }

    # Handle polynomial separately
    if fit_type == 'poly':
        o = int(input('Please insert order of polynomial fit: '))
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

        # Compute R^2 for other fit types
        y_fitted = fitted_func(xdata, *coeff)
        SS_res = np.sum((ydata - y_fitted) ** 2)
        SS_tot = np.sum((ydata - np.mean(ydata)) ** 2)
        R_squared = 1 - (SS_res / SS_tot)

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
    plt.show()

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
            #plt.subplot(1,2,1)
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
