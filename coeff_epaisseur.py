from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees

def generer_graph(filenames, path, title, selected_range="all", uncertainties=(0,0), color="C0"):

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
    pythasson_array = np.zeros(len(filenames))
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
            indice_min = int(indice_centre-15)
            indice_max = int(indice_centre+15)
            label_str = "31 keV"

        elif selected_range == "12keV":
            indice_centre = find_nearest(abscisses_array, 12)
            indice_min = int(indice_centre-15)
            indice_max = int(indice_centre+15)
            label_str = "12 keV"

        elif selected_range == "custom":
            indice_centre = find_nearest(abscisses_array, centre)
            indice_min = int(indice_centre-largeur/2)
            indice_max = int(indice_centre+largeur/2)

        elif type(selected_range) == int:
            indice_centre = find_nearest(abscisses_array, selected_range)
            indice_min = int(indice_centre-largeur/2)
            indice_max = int(indice_centre+largeur/2)
            label_str = f"{selected_range} keV"



        

        data_reduit = data_array[indice_min:indice_max]

        array_poisson = np.sqrt(data_reduit)
        pythasson = 0
        for j in range(len(array_poisson)):
            pythasson += array_poisson[j]**2
        pythasson = np.sqrt(pythasson)/live_time

        somme = np.sum(data_reduit)


        if diviser_par_temps:
            somme = somme/live_time

        tension_array[i] = tension
        courant_array[i] = courant
        livetime_array[i] = live_time
        realtime_array[i] = real_time
        epaisseur_array[i] = epaisseur
        somme_comptes_array[i] = somme
        pythasson_array[i] = pythasson
        incert_array[i] = np.sqrt(somme)

        if i == 0:
            somme0 = somme
            somme0_incert = np.sqrt(somme)

    moy_comptes_array = somme_comptes_array / somme0

    incert_array = moy_comptes_array * np.sqrt((incert_array / somme_comptes_array)**2 + (somme0_incert / somme0)**2)

    ### AFFICHAGE DU GRAPHIQUE ###
    guess = (somme0, 0)
    a, mu, pcov = exponential_fit(epaisseur_array[1:], moy_comptes_array[1:], guess, yerr=incert_array[1:])
    print(f"a = {a}")
    print(f"mu = {mu}")

    err_mu = np.sqrt(np.diag(pcov))[1]

    #### NOTE: mu est ici le coefficient d'atténuation pour le matériau

    fig = plt.errorbar(epaisseur_array, moy_comptes_array, xerr = uncx, yerr = pythasson_array, fmt='o', color=color, capsize=5, markersize=3)
    x_points = np.arange(epaisseur_array[1], epaisseur_array[-1]+0.003, (epaisseur_array[-1]-epaisseur_array[0])/50)
    plt.plot(x_points, exponential(x_points, a, mu), color="C1", label=f"Lissage exponentiel avec incertitude\nNt/N0={a:.2f}*exp(-{mu:.2f}x)")

    # plt.text(epaisseur_array[int(len(epaisseur_array)/2)-2], moy_comptes_array[int(len(moy_comptes_array)/2)+3],f'{a:.2f}*exp(-{mu:.2f}x)', horizontalalignment='center',
    #  verticalalignment='center')



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
fig2 = generer_graph(filenames, path, title, selected_range="all", uncertainties=(uncx, uncy))
plt.legend()
plt.show()

plt.close()

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig3 = generer_graph(filenames, path, title, selected_range=15, uncertainties=(uncx, uncy))
plt.legend()
plt.show()

plt.close()

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig4 = generer_graph(filenames, path, title, selected_range=20, uncertainties=(uncx, uncy))
plt.legend()
plt.show()

plt.close()

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig5 = generer_graph(filenames, path, title, selected_range=30, uncertainties=(uncx, uncy))
plt.legend()
plt.show()

plt.close()

# uncx = np.zeros(len(filenames))
# uncy = np.ones(len(filenames))
# fig6 = generer_graph(filenames, path, title, selected_range=28, uncertainties=(uncx, uncy), color="red")
# plt.legend()
# plt.show()

# plt.close()

# uncx = np.zeros(len(filenames))
# uncy = np.ones(len(filenames))
# fig1 = generer_graph(filenames, path, title, selected_range="31keV", uncertainties=(uncx, uncy), color="purple")

# plt.legend()
# plt.show()

# plt.close()

filenames = sets_donnees.Cu_set
path = "./Data/"

title = "Nb comptes en fct de l'épaisseur de filtre, 50 kV, 15 uA, Cuivre"
uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig2 = generer_graph(filenames, path, title, selected_range="all", uncertainties=(uncx, uncy))
plt.legend()
plt.show()

plt.close()

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig3 = generer_graph(filenames, path, title, selected_range=15, uncertainties=(uncx, uncy))
plt.legend()
plt.show()

plt.close()

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig4 = generer_graph(filenames, path, title, selected_range=20, uncertainties=(uncx, uncy))
plt.legend()
plt.show()

plt.close()

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig5 = generer_graph(filenames, path, title, selected_range=30, uncertainties=(uncx, uncy))
plt.legend()
plt.show()

plt.close()