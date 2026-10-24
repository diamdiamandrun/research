# data loading functions 
from imports import *

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

def saveexcel(x, filename=None):
    if filename is None:
        filename = f"export_{datetime.datetime.now():%Y%m%d_%H%M%S}"
    if not filename.endswith(".csv"):
        filename += ".csv"
    return x.to_csv(filename, index=False)

def copypaste(header=None):
    "Remember to input header = 0 if first row are headers."
    sep = [',', ';', '|', '/', '\\', r'\s+']
    for i in sep:
        try:
            df = pd.read_clipboard(sep=i, header=header)
            if df is not None:
                break
        except:
            pass
    return df
    
def plotarray(filename):
    intm = np.array(filename)
    return [intm[:, i] for i in range(intm.shape[1])]
