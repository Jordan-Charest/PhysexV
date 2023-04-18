from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees

def generer_graph(filenames, path, title, selected_range="all", uncertainties=(0,0)):

    uncx, uncy = uncertainties[0], uncertainties[1]

    # Si custom
    largeur = 30 # Largeur en nombre de canaux à analyser
    centre = 31 # Énergie en keV du centre à analyser
    diviser_par_temps = True # Diviser le nb de comptes par le live time

    tension_array = np.zeros(len(filenames))
    courant_array = np.zeros(len(filenames))
    livetime_array = np.zeros(len(filenames))
    realtime_array = np.zeros(len(filenames))
    epaisseur_array = np.zeros(len(filenames))
    somme_comptes_array = np.zeros(len(filenames))
    somme0 = 0

    for i, filename in enumerate(filenames):

        filepath = f"{path}{filename}.mca" # Nom du fichier à analyser

        try:
            tension, courant, filtre = extraire_params(filename)
            data_array, abscisses_array, live_time, real_time = extraire_data(filepath)
        except:
            continue

        abscisses_array = etalonnage(abscisses_array)
        print(f"live time: {live_time}")
        print(f"real time: {real_time}")

        if filtre == "pas" or filtre == "pos" or filtre == "pus":
            epaisseur = 0
        else:
            index_epp = find_nth_occurrence(filtre, "&", 1)
            epaisseur = int(filtre[index_epp+1:])

        # Selection de la bonne plage

        if selected_range == "all":
            indice_min = find_nearest(abscisses_array, 7.2)
            indice_max = find_nearest(abscisses_array, 50)

        elif selected_range == "31keV":
            indice_centre = find_nearest(abscisses_array, 31)
            indice_min = int(indice_centre-15)
            indice_max = int(indice_centre+15)

        elif selected_range == "12keV":
            indice_centre = find_nearest(abscisses_array, 12)
            indice_min = int(indice_centre-15)
            indice_max = int(indice_centre+15)

        elif selected_range == "custom":
            indice_centre = find_nearest(abscisses_array, centre)
            indice_min = int(indice_centre-largeur/2)
            indice_max = int(indice_centre+largeur/2)


        data_reduit = data_array[indice_min:indice_max]
        somme = np.sum(data_reduit)


        if diviser_par_temps:
            somme = somme/live_time

        tension_array[i] = tension
        courant_array[i] = courant
        livetime_array[i] = live_time
        realtime_array[i] = real_time
        epaisseur_array[i] = epaisseur
        somme_comptes_array[i] = somme

        if i == 0:
            somme0 = somme

    moy_comptes_array = somme_comptes_array / somme0

    ### AFFICHAGE DU GRAPHIQUE ###
    guess = (somme0, 0)
    a, mu = exponential_fit(epaisseur_array, moy_comptes_array, guess)
    print(f"a = {a}")
    print(f"mu = {mu}")

    #### NOTE: mu est ici le coefficient d'atténuation pour le matériau

    fig = plt.errorbar(epaisseur_array, moy_comptes_array, xerr = uncx, yerr = uncy, label=f"{selected_range}, mu={mu:.3f}", fmt='o', capsize=3)
    x_points = np.arange(epaisseur_array[0], epaisseur_array[-1], (epaisseur_array[-1]-epaisseur_array[0])/50)
    plt.plot(x_points, exponential(x_points, a, mu))



    # plt.yscale("log")
    # plt.xscale("log")
    plt.xlabel("Épaisseur de filtre [mil]")

    if diviser_par_temps:
        plt.ylabel("Rapport Nt/N0 moyen par seconde")
    else:
        plt.ylabel("Rapport Nt/N0 moyen total")

    plt.title(title)

    return fig


ranges = ["all", "31keV", "12keV", "custom"]

filenames = sets_donnees.Al2_set
path = "./Data/seance2/"

title = "Nb comptes en fct de l'épaisseur de filtre, 50 kV, 15 uA, Aluminium"

uncx = np.zeros(len(filenames)-1)
uncy = np.ones(len(filenames)-1)
fig1 = generer_graph(filenames[1:], path, title, selected_range="31keV", uncertainties=(uncx, uncy))

uncx = np.zeros(len(filenames)-1)
uncy = np.ones(len(filenames)-1)
fig2 = generer_graph(filenames[1:], path, title, selected_range="all", uncertainties=(uncx, uncy))

uncx = np.zeros(len(filenames)-1)
uncy = np.ones(len(filenames)-1)
fig3 = generer_graph(filenames[1:], path, title, selected_range="12keV", uncertainties=(uncx, uncy))

plt.legend()
plt.show()

plt.close()

filenames = sets_donnees.Cu_set
path = "./Data/"

title = "Nb comptes en fct de l'épaisseur de filtre, 50 kV, 15 uA, Cuivre"
fig1 = generer_graph(filenames, path, title, selected_range="31keV")
fig2 = generer_graph(filenames, path, title, selected_range="all")
fig3 = generer_graph(filenames, path, title, selected_range="12keV")

plt.legend()
plt.show()