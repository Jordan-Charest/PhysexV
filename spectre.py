import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from fonctions import *
from uncertainties import ufloat
from uncertainties import unumpy as unp


### PARAMÈTRES ###

filename = "50_15_Ag&1" # Fichier à charger
multifiltre = True
filepath = f"./Data/seance3/{filename}.mca" # Nom du fichier à analyser
diviser_par_temps = True # Diviser le nb de comptes par le live time



### EXTRACTION DES DONNÉES ###

tension, courant, filtre = extraire_params(filename)
data_array, abscisses_array, live_time, real_time = extraire_data(filepath)
indices_pics = trouver_pic(data_array)
# print(indices_pics)


### TRAITEMENT DES DONNÉES ###

if diviser_par_temps: # Diviser par le live time
    data_array = data_array / live_time

# Étalonnage de l'axe des canaux (abscisses) en énergie
abscisses_array = etalonnage(abscisses_array, 2)



### AFFICHAGE DU GRAPHIQUE ###

data_array = unp.nominal_values(data_array)
fig = plt.plot(abscisses_array, data_array)
# plt.scatter([abscisses_array[indice] for indice in indices_pics], [data_array[indice] for indice in indices_pics])
# for i, indice in enumerate(indices_pics[1:-2]):
#     plt.annotate("%.2f" % abscisses_array[indice], (abscisses_array[indice], data_array[indice]))

plt.yscale("log")
plt.xlabel("Énergie [keV]")

if diviser_par_temps:
    plt.ylabel("Nombre de comptes par seconde [éch. log]")
else:
    plt.ylabel("Nombre de comptes total [log]")

# DEMANDER: nb de comptes, diviser par live time ou real time?

# plt.title(f"Nombre de comptes en fonction de l'énergie, {tension} V, {courant} A, filtres {filtre}")

# plt.axvline(x=10, color="r", label="10 keV")

plt.show()


