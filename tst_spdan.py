import argparse
import sys
import os
import pandas as pd
import numpy as np
import math
from scipy import optimize
from numpy import exp, sqrt, mean, std, diag, inf
from scipy.optimize import curve_fit, fminbound, minimize
import scipy
from scipy.odr import Model, RealData, ODR
#from uncertainties import ufloat, correlated_values
# from scipy.stats import norm
import matplotlib.pyplot as plt
import afterpulse as aft
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})
import spd_analyze

directory = 'C:\\Users\\Андрей\\OneDrive - ООО КуРэйт\\Stand\\V1.2\\!SPDStandResults\\DEFNAME\\1'

SPD = spd_analyze.SPD(directory)