import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt
from uncertainties import ufloat
import uncertainties.umath as math
from uncertainties import unumpy as unp


def extraire_data(filepath, w_uncert=False):

    # np.array pour stocker le nombre de comptes
    data_array = np.empty(4096)

    # np.array d'abscisses (canaux)
    abscisses_array = np.arange(1, 4097, 1)

    # Initialisation des variables de temps
    live_time = 0
    real_time = 0

    # Ouverture du fichier
    f = open(filepath, "r")

    # Extraction des données
    for index, line in enumerate(f):
        if index == 7: # Live time
            ind = find_nth_occurrence(line, "-", 1)
            live_time = ufloat(float(line[ind+2:]), 0.0000005)

        if index == 8: # Real time
            ind = find_nth_occurrence(line, "-", 1)
            real_time = ufloat(float(line[ind+2:]), 0.0000005)

        if index > 4112: # Fin des données
            break

        if index > 16: # Début des données
            if not w_uncert:
                data_array[index-17] = int(line)
            if w_uncert:
                if index == 17:
                    data_array = unp.uarray(data_array, 0)
                data_array[index-17] = ufloat(int(line), np.sqrt(int(line)))

    f.close()

    return data_array, abscisses_array, live_time, real_time


def extraire_params(filename):

    # Indices des underscore pour titres
    index_first_under = find_nth_occurrence(filename, "_", 1)
    index_second_under = find_nth_occurrence(filename, "_", 2)

    # Extraction des valeurs de tension, du courant et des filtres
    tension = filename[:index_first_under]
    courant = filename[index_first_under+1:index_second_under]
    filtre = filename[index_second_under+1:]
    # print(f"Analyse du filtre [{filtre}] à {tension} V et {courant} A")

    return tension, courant, filtre

def fit_w_sigma(xdata, ydata, incert, func):

    popt, pcov = curve_fit(func, xdata, ydata, sigma=incert)

    return popt, pcov

# def linear_fit(xdata, ydata):
#     params = linregress(xdata, ydata)

#     return params.slope, params.intercept

def etalonnage(canaux, num):

    if num == 1:
        params = linregress([0, 921.8, 1173.19, 3895.6], [0, 13.95, 17.24, 59.54])
        # Ce qui nous donne un std error, nice
    elif num == 2:
        params = linregress([0, 931.54, 1173.3, 3896.69], [0, 13.95, 17.24, 59.54])

    return params.slope * canaux + params.intercept


def trouver_pic(data):
    indices = find_peaks_cwt(data, 4, min_snr=3)
    #indices = find_peaks(data, height=7, threshold=20, distance=50)[0]

    return indices


def gaussian(x, sigma, mu, H):
    return H + 1/(sigma * np.sqrt(2*np.pi))* np.exp(-0.5 * (x-mu)**2 / sigma**2)


def gaussian_fit(xdata, ydata, guess, yerr=0):

    popt, pcov = curve_fit(gaussian, xdata, ydata, p0=guess, sigma=yerr, absolute_sigma=True)


    return popt, pcov

def exponential(x, a, b):
    return a * np.exp(-b * x)
    # return np.exp(-b*x)

def linear(x, a, b):
    return a*x + b

def atau_Z(Z, c, n):
    return c * Z ** n

def exponential_fit(xdata, ydata, guess, func=exponential, yerr=0):

    if type(yerr) == int:
        yerr = np.zeros(len(xdata))
    # xdata = [float(i) for i in xdata]
    # ydata = [float(i) for i in ydata]
    # yerr = [float(i) for i in yerr]

    xdata = xdata.astype(dtype=np.float32)
    ydata = ydata.astype(dtype=np.float32)
    yerr = yerr.astype(dtype=np.float32)
    popt, pcov = curve_fit(func, xdata, ydata, p0=guess, sigma=yerr, absolute_sigma=True)

    return popt[0], popt[1], pcov

def linear_fit(xdata, ydata, guess, func=linear, yerr=0):

    popt, pcov = curve_fit(func, xdata, ydata, p0=guess, absolute_sigma=True)

    return popt[0], popt[1], pcov


def find_nearest(array, value):
    index = np.argmin(np.abs(array-value))
    return index

def a_tau(N0, Nt, A, rho, t):

    # print(f"A={A}, rho={rho}, t={t}")

    t_incert = 0.06375
    N0_val = N0.nominal_value
    N0_incert = N0.std_dev
    Nt_val = Nt.nominal_value
    Nt_incert = Nt.std_dev

    # print(f"N0={N0_val}±{N0_incert}"); print(f"Nt={Nt_val}±{Nt_incert}")

    N = Nt/N0
    # N_incert = N * np.sqrt((N0_incert/N0_val)**2 + (Nt_incert/Nt_val)**2)
    # print(f"N = {N}±{N_incert}")
    # ln_t = np.log(N)/t
    # incert_ln_t = ln_t * np.sqrt((N_incert/N)**2 + (t_incert/t)**2)
    # incert_ln_t = N_incert

    avo = 6.022e23
    
    return -math.log(N) * A / (rho * t * avo) - 0.2 * A / avo


def find_nth_occurrence(string, substring, n):
    # Finds the index of the nth occurrence of substring in string
    # Tiré de StackOverflow: https://stackoverflow.com/questions/1883980/find-the-nth-occurrence-of-substring-in-a-string
    start = string.find(substring)
    while start >= 0 and n > 1:
        start = string.find(substring, start+len(substring))
        n -= 1
    return start