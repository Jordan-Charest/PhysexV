from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees


filenames = sets_donnees.courant_set
diviser_par_temps = True # Diviser le nb de comptes par le live time
section = True # Analyser seulement la section précisée ci-bas
largeur = 20 # Largeur en nombre de canaux à analyser
centre = 20 # Énergie en keV du centre à analyser
filtre_nom = "Courant"

tension_array = np.zeros(len(filenames))
courant_array = np.zeros(len(filenames))
livetime_array = np.zeros(len(filenames))
realtime_array = np.zeros(len(filenames))
epaisseur_array = np.zeros(len(filenames))
moy_comptes_array = np.zeros(len(filenames))
moy0 = 0

for i, filename in enumerate(filenames):

    filepath = f"./Data/{filename}.mca" # Nom du fichier à analyser
    tension, courant, filtre = extraire_params(filename)
    data_array, abscisses_array, live_time, real_time = extraire_data(filepath)
    abscisses_array = etalonnage(abscisses_array)

    if filtre == "pas":
        epaisseur = 0
    else:
        epaisseur = int(filtre[2:])

    if section:
        indice_centre = find_nearest(abscisses_array, centre)
        indice_min = int(indice_centre-largeur/2)
        indice_max = int(indice_centre+largeur/2)

        data_reduit = data_array[indice_min:indice_max]
        moy = np.sum(data_reduit)
    else:
        moy = np.sum(data_array)

    tension_array[i] = tension.replace(",", ".")
    courant_array[i] = courant.replace(",", ".")
    livetime_array[i] = live_time
    realtime_array[i] = real_time
    epaisseur_array[i] = epaisseur
    moy_comptes_array[i] = moy

    if i == 0:
        moy0 = moy

moy_comptes_array = moy_comptes_array / moy0

### AFFICHAGE DU GRAPHIQUE ###


fig = plt.scatter(courant_array, moy_comptes_array)




# plt.yscale("log")
# plt.xscale("log")
plt.xlabel("Courant [uA]")

if diviser_par_temps:
    plt.ylabel("Somme du nombre de comptes par seconde")
else:
    plt.ylabel("Somme du nombre de comptes total")

# DEMANDER: nb de comptes, diviser par live time ou real time?

plt.title(f"Nombre de comptes en fonction du courant, {tension} kV")

plt.show()
