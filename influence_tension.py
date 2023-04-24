from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees



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
    pythasson_array = np.zeros(len(filenames))
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
            epaisseur = int(filtre[index_epp+1:].replace(",", "."))

        # Selection de la bonne plage

        if selected_range == "all":
            # indice_min = find_nearest(abscisses_array, 5)
            indice_min = 0
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

        array_poisson = np.sqrt(data_reduit)
        pythasson = 0
        for j in range(len(array_poisson)):
            pythasson += array_poisson[j]**2
        pythasson = np.sqrt(pythasson)/live_time

        somme = np.sum(data_reduit)


        if diviser_par_temps:
            somme = somme/live_time

        tension_array[i] = float(tension.replace(",", "."))
        courant_array[i] = float(courant.replace(",", "."))
        livetime_array[i] = live_time
        realtime_array[i] = real_time
        epaisseur_array[i] = epaisseur
        somme_comptes_array[i] = somme
        pythasson_array[i] = pythasson

        if i == 0:
            somme0 = somme

    # moy_comptes_array = somme_comptes_array / somme0
    moy_comptes_array = somme_comptes_array

    pythasson_array = moy_comptes_array * np.sqrt((pythasson_array[i]/somme_comptes_array[i])**2 + (pythasson_array[0]/somme0)**2)

    ### AFFICHAGE DU GRAPHIQUE ###
    # guess = (somme0, 0)
    # a, mu = exponential_fit(epaisseur_array, moy_comptes_array, guess)
    # print(f"a = {a}")
    # print(f"mu = {mu}")

    #### NOTE: mu est ici le coefficient d'atténuation pour le matériau

    fig = plt.scatter(tension_array, moy_comptes_array, label=f"{selected_range}")

    # plt.yscale("log")
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.xscale("log")
    plt.xlabel("Tension [kV]")

    plt.ylabel("Somme du nombre de comptes total")

# DEMANDER: nb de comptes, diviser par live time ou real time?

    plt.title(f"Nombre de comptes en fonction de la tension, {courant} uA")

    return fig


ranges = ["all", "31keV", "12keV", "custom"]

filenames = sets_donnees.tension_set
path = "./Data/"

title = "Nb comptes en fct de l'épaisseur de filtre, 50 kV, 15 uA, Aluminium"

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig1 = generer_graph(filenames, path, title, selected_range="31keV", uncertainties=(uncx, uncy))

plt.legend()
plt.show()

plt.close()

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig2 = generer_graph(filenames, path, title, selected_range="all", uncertainties=(uncx, uncy))

plt.legend()
plt.show()

plt.close()

uncx = np.zeros(len(filenames))
uncy = np.ones(len(filenames))
fig3 = generer_graph(filenames, path, title, selected_range="12keV", uncertainties=(uncx, uncy))

plt.legend()
plt.show()

plt.close()


###############################################################################
