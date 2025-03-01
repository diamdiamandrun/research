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
