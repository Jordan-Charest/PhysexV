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
    print(f"Analyse du filtre [{filtre}] à {tension} V et {courant} A")

    return tension, courant, filtre

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

def atau_Z(Z, c, n):
    return c * Z ** n

def exponential_fit(xdata, ydata, guess, func=exponential):
    
    popt, pcov = curve_fit(func, xdata, ydata, p0=guess)

    return popt[0], popt[1]


def find_nearest(array, value):
    index = np.argmin(np.abs(array-value))
    return index

def a_tau(N, A, rho, t):
    print(f"N={N}, A={A}, rho={rho}, t={t}")
    avo = 6.022e23
    return -np.log(N) * A / (rho * t * avo) - 0.2 * A / avo


def find_nth_occurrence(string, substring, n):
    # Finds the index of the nth occurrence of substring in string
    # Tiré de StackOverflow: https://stackoverflow.com/questions/1883980/find-the-nth-occurrence-of-substring-in-a-string
    start = string.find(substring)
    while start >= 0 and n > 1:
        start = string.find(substring, start+len(substring))
        n -= 1
    return start