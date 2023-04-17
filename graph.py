import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from fonctions import *

### PARAMÈTRES ###

filename = "50_15_pas" # Fichier à charger
filepath = f"./Data/{filename}.mca" # Nom du fichier à analyser
diviser_par_temps = True # Diviser le nb de comptes par le live time



### EXTRACTION DES DONNÉES ###

tension, courant, filtre = extraire_params(filename)
data_array, abscisses_array, live_time, real_time = extraire_data(filepath)
# indices_pics = trouver_pic(data_array)
# print(indices_pics)



### TRAITEMENT DES DONNÉES ###

if diviser_par_temps: # Diviser par le live time
    data_array = data_array / live_time

# Étalonnage de l'axe des canaux (abscisses) en énergie
abscisses_array = etalonnage(abscisses_array)



### AFFICHAGE DU GRAPHIQUE ###

fig = plt.plot(abscisses_array, data_array)
# plt.scatter([abscisses_array[indice] for indice in indices_pics], [data_array[indice] for indice in indices_pics])
# for i, indice in enumerate(indices_pics):
#     plt.annotate("%.2f" % abscisses_array[indice], (abscisses_array[indice], data_array[indice]))

plt.yscale("log")
plt.xlabel("Énergie [keV]")

if diviser_par_temps:
    plt.ylabel("Nombre de comptes par seconde [log]")
else:
    plt.ylabel("Nombre de comptes total [log]")

# DEMANDER: nb de comptes, diviser par live time ou real time?

plt.title(f"Nombre de comptes en fonction de l'énergie, {tension} V, {courant} A, filtres {filtre}")

plt.show()


