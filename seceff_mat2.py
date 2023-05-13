from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees
from decimal import Decimal

def generer_graph(filenames, path, title, selected_range="all", uncertainties=(0,0)):
    
    uncx, uncy = uncertainties[0], uncertainties[1]

    array_dict = [sets_donnees.aluminium, sets_donnees.cuivre, sets_donnees.molybdene, sets_donnees.argent, sets_donnees.tungstene]
    # array_dict = [sets_donnees.aluminium, sets_donnees.cuivre, sets_donnees.molybdene, sets_donnees.tungstene]
    Z_array = [dico["Z"] for dico in array_dict]
    A_array = [dico["A"] for dico in array_dict]
    rho_array = [dico["rho"] for dico in array_dict]

    # Si custom
    largeur = 30 # Largeur en nombre de canaux à analyser
    centre = 31 # Énergie en keV du centre à analyser
    diviser_par_temps = False # Diviser le nb de comptes par le live time

    tension_array = np.zeros(len(filenames))
    courant_array = np.zeros(len(filenames))
    livetime_array = np.zeros(len(filenames))
    realtime_array = np.zeros(len(filenames))
    epaisseur_array = np.zeros(len(filenames))
    somme_comptes_array = np.zeros(len(filenames))
    incert_y_array = np.zeros(len(filenames))
    somme0 = 0

    for i, filename in enumerate(filenames):

        filepath = f"{path}{filename}.mca" # Nom du fichier à analyser

        try:
            tension, courant, filtre = extraire_params(filename)
            data_array, abscisses_array, live_time, real_time = extraire_data(filepath)
        except:
            print("Erreur de chargement")
            continue

        abscisses_array = etalonnage(abscisses_array)
        # print(f"live time: {live_time}")
        # print(f"real time: {real_time}")

        if filtre == "pas":
            epaisseur = 0
        else:
            index_epp = find_nth_occurrence(filtre, "&", 1)
            epaisseur = int(filtre[index_epp+1:])

        # Selection de la bonne plage

        if selected_range == "all":
            indice_min = find_nearest(abscisses_array, 5)
            indice_max = find_nearest(abscisses_array, 50)

        elif selected_range == "31keV":
            plage = "31±0,23 keV"
            indice_centre = find_nearest(abscisses_array, 31)
            indice_min = int(indice_centre-15)
            indice_max = int(indice_centre+15)

        elif selected_range == "10keV":
            plage = "10±0,23 keV"
            indice_centre = find_nearest(abscisses_array, 10)
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

        tension_array[i] = tension.replace(",", ".")
        courant_array[i] = courant.replace(",", ".")
        livetime_array[i] = live_time
        realtime_array[i] = real_time
        epaisseur_array[i] = epaisseur
        somme_comptes_array[i] = somme

        print(f"N{i}={somme_comptes_array[i]}")

        if i == 0:
            somme0 = somme

    # print(f"N0={somme0}")

    incert_y_array = np.sqrt(somme_comptes_array)
    # print(f"ARRAY Y: {somme_comptes_array}")
    # print(f"INCERTITUDE SUR ARRAY Y: {incert_y_array}")
    N0_tuple = (somme_comptes_array[0], incert_y_array[0])

    # somme_comptes_array = somme_comptes_array / somme0

    tension_array = tension_array[1:]
    courant_array = courant_array[1:]
    livetime_array = livetime_array[1:]
    realtime_array = realtime_array[1:]
    epaisseur_array = epaisseur_array[1:]
    somme_comptes_array = somme_comptes_array[1:]
    incert_y_array = incert_y_array[1:]
    

    # ### AFFICHAGE DU GRAPHIQUE ###
    # guess = (somme0, 0)
    # a, mu = exponential_fit(epaisseur_array, somme_comptes_array, guess)
    # print(f"a = {a}")
    # print(f"mu = {mu}")

    atau_array = np.zeros(len(array_dict))
    atau_incert_array = np.zeros(len(array_dict))
    for i in range(len(array_dict)):
        # print(f"Matériau: Z={Z_array[i]}")
        (atau, atau_incert) = a_tau(N0_tuple, (somme_comptes_array[i], incert_y_array[i]), A_array[i], rho_array[i], epaisseur_array[i]*.002540)
        # print(atau_incert)
        atau_array[i] = atau
        atau_incert_array[i] = atau_incert

    print(f"Range: {selected_range}")
    print(f"ARRAY ATAU: {atau_array}")

    # On met les données en log-log pour faire un fit linéaire
    Z_array_log = np.log10(Z_array)
    atau_array_log = np.log10(atau_array)

    slope, intercept, pcov = linear_fit(Z_array_log, atau_array_log, (0, 0))
    start = Z_array_log[0]
    stop = Z_array_log[-1]
    step = (stop-start)/100
    xrange = np.arange(start, stop, step)
    yrange = xrange*slope + intercept


    # ### AFFICHAGE DU GRAPHIQUE ###




    fig = plt.scatter(Z_array_log, atau_array_log)
    plt.plot(xrange, yrange)


    plt.title(title)

    # plt.yscale("log")
    # plt.xscale("log")
    plt.xlabel("Z")
    plt.ylabel("atau")

    return fig

#####################################

    # if selected_range=="all":
    #     fig = plt.scatter(Z_array_log, atau_array_log, label=f"a_tau, tout le spectre")
    #     # fig = plt.scatter(Z_array, atau_array, label=f"a_tau, tout le spectre")
    # else:
    #     fig = plt.scatter(Z_array_log, atau_array_log, label=f"a_tau, plage de {plage}")
    #     # fig = plt.scatter(Z_array, atau_array, label=f"a_tau, plage de {plage}")
    # x_steps = np.arange(np.floor(np.log10(10)), np.ceil(np.log10(75)), 1)
    # # plt.plot(x_steps, atau_Z(x_steps, c, n), label=f"Lissage exponentiel, a_tau={'%.2E' % Decimal(c)}*Z^{n:.2f}")
    # # fig = plt.errorbar(Z_array, atau_array, xerr = 0, yerr = atau_incert_array, label=f"{selected_range}", fmt='o', capsize=3, markersize=3)
    # # fig = plt.scatter(Z_array, atau_array, label=f"{selected_range}")
    # # for i, Zlog in enumerate(Z_array_log):
    # #     plt.annotate(f"{10**Zlog:.0f}", (Zlog, atau_array_log[i]) )
    # plt.plot(x_steps, (x_steps*slope + intercept))
    


ranges = ["all", "31keV", "12keV", "custom"]

filenames = sets_donnees.materiaux_set2
path = "./Data/Filtres/"

# filenames = sets_donnees.materiaux_set_pas_ag
# path = "./Data/seance1/"

title = "a_tau en fct du Z du matériau, 50 kV, 15 uA"

uncx = np.zeros(len(filenames)-1)
uncy = np.ones(len(filenames)-1)
fig3 = generer_graph(filenames, path, title, selected_range="10keV", uncertainties=(uncx, uncy))

# uncx = np.zeros(len(filenames)-1)
# uncy = np.zeros(len(filenames)-1)
# fig1 = generer_graph(filenames, path, title, selected_range="31keV", uncertainties=(uncx, uncy))

# uncx = np.zeros(len(filenames)-1)
# uncy = np.ones(len(filenames)-1)
# fig2 = generer_graph(filenames, path, title, selected_range="all", uncertainties=(uncx, uncy))



plt.legend()
plt.show()

# ==========================================================



# tension_array = np.zeros(len(filenames))
# courant_array = np.zeros(len(filenames))
# livetime_array = np.zeros(len(filenames))
# realtime_array = np.zeros(len(filenames))
# epaisseur_array = np.zeros(len(filenames))
# moy_comptes_array = np.zeros(len(filenames))
# moy0 = 0

# for i, filename in enumerate(filenames):

#     filepath = f"./Data/{filename}.mca" # Nom du fichier à analyser
#     tension, courant, filtre = extraire_params(filename)

#     data_array, abscisses_array, live_time, real_time = extraire_data(filepath)
#     abscisses_array = etalonnage(abscisses_array)

#     if filtre == "pas":
#         epaisseur = 0
#     else:
#         index_epp = find_nth_occurrence(filtre, "&", 1)
#         epaisseur = int(filtre[index_epp+1:])

#     if section:
#         indice_centre = find_nearest(abscisses_array, centre)
#         indice_min = int(indice_centre-largeur/2)
#         indice_max = int(indice_centre+largeur/2)

#         data_reduit = data_array[indice_min:indice_max]
#         moy = np.sum(data_reduit)
#     else:
#         moy = np.sum(data_array)

#     tension_array[i] = tension.replace(",", ".")
#     courant_array[i] = courant.replace(",", ".")
#     livetime_array[i] = live_time
#     realtime_array[i] = real_time
#     epaisseur_array[i] = epaisseur
#     moy_comptes_array[i] = moy

#     if i == 0:
#         moy0 = moy

# moy_comptes_array = moy_comptes_array / moy0

# tension_array = tension_array[1:]
# courant_array = courant_array[1:]
# livetime_array = livetime_array[1:]
# realtime_array = realtime_array[1:]
# epaisseur_array = epaisseur_array[1:]
# moy_comptes_array = moy_comptes_array[1:]

# atau_array = np.zeros(len(array_dict))
# for i in range(len(array_dict)):
#     print(f"Matériau: Z={Z_array[i]}")
#     atau = a_tau(moy_comptes_array[i], A_array[i], rho_array[i], epaisseur_array[i]*0.00254)
#     atau_array[i] = atau

# Z_array_log = np.log10(Z_array)
# atau_array_log = np.log10(atau_array)

# slope, intercept = linear_fit(Z_array_log, atau_array_log)
# print(slope, intercept)

# # guess = (0,4)
# # c, n = exponential_fit(Z_array, atau_array, guess, func=atau_Z)
# # print(f"c={c}, n={n}")

# ### AFFICHAGE DU GRAPHIQUE ###





# # fig = plt.scatter(Z_array, atau_array)
# x_steps = np.arange(1, 74, 1)
# # plt.plot(x_steps, atau_Z(x_steps, c, n))
# fig = plt.scatter(Z_array_log, atau_array_log)
# plt.plot(Z_array_log, Z_array_log*slope + intercept)




# # plt.yscale("log")
# # plt.xscale("log")
# plt.xlabel("log10(Z)")
# plt.ylabel("log10(atau)")

