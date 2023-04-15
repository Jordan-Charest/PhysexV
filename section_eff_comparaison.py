from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees


filenames = sets_donnees.Al_set
diviser_par_temps = True # Diviser le nb de comptes par le live time
largeur = 50 # Largeur en nombre de canaux à analyser
centre = 31 # Énergie en keV du centre à analyser
filtre_nom = "Cu"
A = 63.55 # Cu
rho = 8.933 # Cu, g/cm3

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

    indice_centre = find_nearest(abscisses_array, centre)
    indice_min = int(indice_centre-largeur/2)
    indice_max = int(indice_centre+largeur/2)

    data_reduit = data_array[indice_min:indice_max]
    moy = np.mean(data_reduit)

    tension_array[i] = tension
    courant_array[i] = courant
    livetime_array[i] = live_time
    realtime_array[i] = real_time
    epaisseur_array[i] = epaisseur
    moy_comptes_array[i] = moy

    if i == 0:
        moy0 = moy

moy_comptes_array = moy_comptes_array / moy0

# array_a_tau = moy_comptes

### AFFICHAGE DU GRAPHIQUE ###


#### NOTE: mu est ici le coefficient d'atténuation pour le matériau

fig = plt.scatter(epaisseur_array, moy_comptes_array)

x_points = np.arange(epaisseur_array[0], epaisseur_array[-1], (epaisseur_array[-1]-epaisseur_array[0])/50)
plt.plot(x_points, exponential(x_points, a, mu))



# plt.yscale("log")
plt.xlabel("Épaisseur de filtre [mil]")

if diviser_par_temps:
    plt.ylabel("Rapport Nt/N0 moyen par seconde")
else:
    plt.ylabel("Rapport Nt/N0 moyen total")

# DEMANDER: nb de comptes, diviser par live time ou real time?

plt.title(f"Nombre de comptes en fonction de l'épaisseur de filtre, {tension} V, {courant} A, filtres {filtre_nom}")

plt.show()
