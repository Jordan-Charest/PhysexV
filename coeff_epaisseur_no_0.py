from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees

def generer_graph(filenames, path, title, selected_range="all", uncertainties=(0,0), color="black"):

    uncx, uncy = uncertainties[0], uncertainties[1]

    # Si custom
    largeur = 25 # Largeur en nombre de canaux à analyser
    centre = 31 # Énergie en keV du centre à analyser
    diviser_par_temps = True # Diviser le nb de comptes par le live time

    tension_array = np.zeros(len(filenames))
    courant_array = np.zeros(len(filenames))
    livetime_array = np.zeros(len(filenames))
    realtime_array = np.zeros(len(filenames))
    epaisseur_array = np.zeros(len(filenames))
    incert_array = np.zeros(len(filenames))
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
            epaisseur = int(filtre[index_epp+1:]) * 0.00254

        # Selection de la bonne plage

        if selected_range == "all":
            indice_min = find_nearest(abscisses_array, 5)
            indice_max = find_nearest(abscisses_array, 50)
            label_str = "Tout le spectre"

        elif selected_range == "31keV":
            indice_centre = find_nearest(abscisses_array, 31)
            indice_min = int(indice_centre-largeur)
            indice_max = int(indice_centre+largeur)
            label_str = "31 keV"

        elif selected_range == "12keV":
            indice_centre = find_nearest(abscisses_array, 12)
            indice_min = int(indice_centre-largeur)
            indice_max = int(indice_centre+largeur)
            label_str = "12 keV"

        elif selected_range == "custom":
            indice_centre = find_nearest(abscisses_array, centre)
            indice_min = int(indice_centre-largeur)
            indice_max = int(indice_centre+largeur)

        elif type(selected_range) == int:
            indice_centre = find_nearest(abscisses_array, selected_range)
            indice_min = int(indice_centre-largeur)
            indice_max = int(indice_centre+largeur)
            label_str = f"{selected_range} keV"



        

        data_reduit = data_array[indice_min:indice_max]
        

        somme = np.sum(data_reduit)


        if diviser_par_temps:
            somme = somme/live_time

        somme_comptes_array_all = somme_comptes_array
        epaisseur_array_all = epaisseur_array
        incert_array_all = incert_array

        tension_array[i] = tension
        courant_array[i] = courant
        livetime_array[i] = live_time
        realtime_array[i] = real_time
        epaisseur_array[i] = epaisseur
        somme_comptes_array[i] = somme
        incert_array[i] = np.sqrt(somme)

    #     if i == 1:
    #         somme1 = somme
    #         somme1_incert = np.sqrt(somme)

        if i == 0:
            somme1 = somme
            somme1_incert = np.sqrt(somme)


    # tension_array = tension_array[1:]
    # courant_array = courant_array[1:]
    # livetime_array = livetime_array[1:]
    # realtime_array = realtime_array[1:]
    # epaisseur_array = epaisseur_array[1:]
    # somme_comptes_array = somme_comptes_array[1:]
    # incert_array = incert_array[1:]

    moy_comptes_array = somme_comptes_array / somme1
    incert_array = moy_comptes_array * np.sqrt((incert_array / somme_comptes_array)**2 + (somme1_incert / somme1)**2)

    moy_comptes_array_all = somme_comptes_array_all / somme1
    incert_array_all = moy_comptes_array_all * np.sqrt((incert_array_all / somme_comptes_array_all)**2 + (somme1_incert / somme1)**2)

    # pythasson_array = moy_comptes_array * np.sqrt((pythasson_array/somme_comptes_array)**2 + (pythasson_array[0]/somme0)**2)

    ### AFFICHAGE DU GRAPHIQUE ###
    guess = (somme1, 0)
    a, mu, pcov = exponential_fit(epaisseur_array, moy_comptes_array, guess, yerr=incert_array)
    print(f"a = {a}")
    print(f"mu = {mu}")

    err_mu = np.sqrt(np.diag(pcov))[1]

    #### NOTE: mu est ici le coefficient d'atténuation pour le matériau

    fig = plt.errorbar(epaisseur_array_all, moy_comptes_array_all, xerr = 0, yerr = incert_array_all, label=f"{label_str}, mu={mu:.2f}±{err_mu:.2f}", fmt='o', color=color, capsize=5, markersize=3)
    # fig = plt.scatter(epaisseur_array, moy_comptes_array, label=f"{label_str}, mu={mu:.3f}", color=color)
    x_points = np.arange(epaisseur_array[0], epaisseur_array[-1], (epaisseur_array[-1]-epaisseur_array[0])/50)
    plt.plot(x_points, exponential(x_points, a, mu), color=color)



    # plt.yscale("log")
    # plt.xscale("log")
    plt.xlabel("Épaisseur de filtre [cm]")

    if diviser_par_temps:
        plt.ylabel("Rapport Nt/N0 par seconde")
    else:
        plt.ylabel("Rapport Nt/N0 moyen total")

    plt.title(title)

    return fig


ranges = ["all", "31keV", "12keV", "custom"]

filenames = sets_donnees.Al2_set
path = "./Data/seance2/"

title = "Nb comptes en fct de l'épaisseur de filtre, 50 kV, 15 uA, Aluminium"

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig2 = generer_graph(filenames, path, title, selected_range="all", uncertainties=(uncx, uncy), color="cyan")

# uncx = np.zeros(len(filenames))
# uncy = np.ones(len(filenames))
# fig3 = generer_graph(filenames, path, title, selected_range=10, uncertainties=(uncx, uncy), color="blue")

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig4 = generer_graph(filenames, path, title, selected_range=15, uncertainties=(uncx, uncy), color="green")

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig5 = generer_graph(filenames, path, title, selected_range=20, uncertainties=(uncx, uncy), color="orange")

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig6 = generer_graph(filenames, path, title, selected_range=30, uncertainties=(uncx, uncy), color="red")

# uncx = np.zeros(len(filenames))
# uncy = np.ones(len(filenames))
# fig1 = generer_graph(filenames, path, title, selected_range="31keV", uncertainties=(uncx, uncy), color="purple")

plt.legend()
plt.show()

plt.close()

filenames = sets_donnees.Cu_set
path = "./Data/"

title = "Nb comptes en fct de l'épaisseur de filtre, 50 kV, 15 uA, Cuivre"
uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig2 = generer_graph(filenames, path, title, selected_range="all", uncertainties=(uncx, uncy), color="cyan")

# uncx = np.zeros(len(filenames))
# uncy = np.ones(len(filenames))
# fig3 = generer_graph(filenames, path, title, selected_range=10, uncertainties=(uncx, uncy), color="blue")

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig4 = generer_graph(filenames, path, title, selected_range=15, uncertainties=(uncx, uncy), color="green")

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig5 = generer_graph(filenames, path, title, selected_range=20, uncertainties=(uncx, uncy), color="orange")

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig6 = generer_graph(filenames, path, title, selected_range=30, uncertainties=(uncx, uncy), color="red")

plt.legend()
plt.show()