import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt


def extraire_data(filepath):

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
            live_time = float(line[ind+2:])

        if index == 8: # Real time
            ind = find_nth_occurrence(line, "-", 1)
            real_time = float(line[ind+2:])

        if index > 4112: # Fin des données
            break

        if index > 16: # Début des données
            data_array[index-17] = int(line)

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

def linear_fit(xdata, ydata):
    params = linregress(xdata, ydata)

    return params.slope, params.intercept

def etalonnage(canaux):

    params = linregress([0, 921.8, 1173.19, 3895.6], [0, 13.95, 17.24, 59.54])
    # Ce qui nous donne un std error, nice

    return params.slope * canaux + params.intercept


def trouver_pic(data):
    indices = find_peaks_cwt(data, 5, min_snr=3)
    #indices = find_peaks(data, height=7, threshold=20, distance=50)[0]

    return indices


def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x-x0)**2 / (2 * sigma ** 2))


def gaussian_fit(xdata, ydata):
    return

def exponential(x, a, b):
    return a * np.exp(-b * x)

def linear(x, a, b):
    return a*x + b

def atau_Z(Z, c, n):
    return c * Z ** n

def exponential_fit(xdata, ydata, guess, func=exponential, yerr=0):
    
    popt, pcov = curve_fit(func, xdata, ydata, p0=guess, sigma=yerr, absolute_sigma=True)

    return popt[0], popt[1], pcov


def find_nearest(array, value):
    index = np.argmin(np.abs(array-value))
    return index

def a_tau(N0, Nt, A, rho, t):

    # print(f"A={A}, rho={rho}, t={t}")

    t_incert = 0.06375
    N0_val = N0[0]
    N0_incert = N0[1]
    Nt_val = Nt[0]
    Nt_incert = Nt[1]

    # print(f"N0={N0_val}±{N0_incert}"); print(f"Nt={Nt_val}±{Nt_incert}")

    N = Nt_val/N0_val
    N_incert = N * np.sqrt((N0_incert/N0_val)**2 + (Nt_incert/Nt_val)**2)
    # print(f"N = {N}±{N_incert}")
    ln_t = np.log(N)/t
    incert_ln_t = ln_t * np.sqrt((N_incert/N)**2 + (t_incert/t)**2)

    avo = 6.022e23
    
    return (-np.log(N) * A / (rho * t * avo) - 0.2 * A / avo, (incert_ln_t)*A/(rho*avo))


def find_nth_occurrence(string, substring, n):
    # Finds the index of the nth occurrence of substring in string
    # Tiré de StackOverflow: https://stackoverflow.com/questions/1883980/find-the-nth-occurrence-of-substring-in-a-string
    start = string.find(substring)
    while start >= 0 and n > 1:
        start = string.find(substring, start+len(substring))
        n -= 1
    return start