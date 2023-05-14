from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees
from uncertainties import ufloat
from uncertainties import unumpy as unp
from decimal import Decimal

def generer_graph(filenames, path, title, selected_range="all", uncertainties=(0,0), color="C0", mat="Al"):

    uncx, uncy = uncertainties[0], uncertainties[1]

    # Si custom
    largeur = 30 # Largeur en nombre de canaux à analyser
    centre = 31 # Énergie en keV du centre à analyser
    diviser_par_temps = True # Diviser le nb de comptes par le live time

    tension_array = np.zeros(len(filenames))
    tension_array = unp.uarray(tension_array, 0)

    courant_array = np.zeros(len(filenames))
    courant_array = unp.uarray(courant_array, 0)

    livetime_array = np.zeros(len(filenames))
    livetime_array = unp.uarray(livetime_array, 0)

    realtime_array = np.zeros(len(filenames))
    realtime_array = unp.uarray(realtime_array, 0)

    epaisseur_array = np.zeros(len(filenames))
    epaisseur_array = unp.uarray(epaisseur_array, 0)

    # print(epaisseur_array)

    # pythasson_array = np.zeros(len(filenames))
    # incert_array = np.zeros(len(filenames))

    somme_comptes_array = np.zeros(len(filenames))
    somme_comptes_array = unp.uarray(somme_comptes_array, 0)

    for i, filename in enumerate(filenames):
        # print(f"i={i}")

        filepath = f"{path}{filename}.mca" # Nom du fichier à analyser

        try:
            tension, courant, filtre = extraire_params(filename)
            data_array, abscisses_array, live_time, real_time = extraire_data(filepath, w_uncert=True)
        except:
            raise ValueError("Ayoye")
            continue

        if mat == "Cu":
            abscisses_array = etalonnage(abscisses_array, 1)
        elif mat == "Al":
            abscisses_array = etalonnage(abscisses_array, 2)
        # print(f"live time: {live_time}")
        # print(f"real time: {real_time}")

        if filtre == "pas" or filtre == "pos" or filtre == "pus":
            epaisseur = ufloat(0, 0)
        else:
            index_epp = find_nth_occurrence(filtre, "&", 1)
            ep = int(filtre[index_epp+1:]) * 0.00254
            epaisseur = ufloat(ep, 0.06375*0.00254)

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

        # array_poisson = np.sqrt(data_reduit)
        # pythasson = 0
        # for j in range(len(array_poisson)):
        #     pythasson += array_poisson[j]**2
        # pythasson = np.sqrt(pythasson)/live_time

        somme = data_reduit.sum()
        # incert_somme = np.sqrt(somme)

        somme = ufloat(somme.nominal_value, somme.std_dev) / live_time


        # if diviser_par_temps:
        #     somme = somme/live_time

        tension_array[i] = tension
        courant_array[i] = courant
        livetime_array[i] = live_time
        realtime_array[i] = real_time
        epaisseur_array[i] = epaisseur
        somme_comptes_array[i] = somme
        # pythasson_array[i] = pythasson
        # incert_array[i] = np.sqrt(somme)
        # print(f"i={i}")
        if i == 0:
            somme0 = somme

    # comptes_array_normalise = somme_comptes_array

    # for i in comptes_array_normalise:
    #     i = ufloat(i.nominal_value/somme0.nominal_value, i.nominal_value * np.sqrt((somme0.std_dev/somme0.nominal_value)**2+(i.std_dev/i.nominal_value)**2))

    comptes_array_normalise = somme_comptes_array / somme0

    # print(comptes_array_normalise)

    ### AFFICHAGE DU GRAPHIQUE ###

    comptes_array_normalise[0] = ufloat(comptes_array_normalise[0].nominal_value, np.sqrt(somme_comptes_array[0].nominal_value)/somme0.nominal_value)


    epaisseur_vals = np.array(unp.nominal_values(epaisseur_array))
    comptes_vals = np.array(unp.nominal_values(comptes_array_normalise))
    print(comptes_vals)
    comptes_incert = np.array(unp.std_devs(comptes_array_normalise))
    print(comptes_incert)
    guess = (10)
    # a, mu, pcov = exponential_fit(epaisseur_vals, comptes_vals, guess, yerr=comptes_incert)
    mu, pcov = exponential_fit(epaisseur_vals, comptes_vals, guess, yerr=comptes_incert)
    # print(f"a = {a}")
    mu = float(mu[0])
    print(f"mu = {mu}")
    # err_mu = np.sqrt(np.diag(pcov))[1]
    err_mu = np.sqrt(pcov)[0][0]
    print(f"err mu = {err_mu}")
    # err_a = np.sqrt(pcov[0][0])
    # err_mu = np.sqrt(pcov[1][1])

    x_points = np.arange(epaisseur_vals[0], epaisseur_vals[-1]+0.003, (epaisseur_vals[-1]-epaisseur_vals[0])/50)
    y_points = exponential(x_points, mu)

    # indice_CDA1 = find_nearest(y_points, 0.5)
    # indice_CDA2 = find_nearest(y_points, 0.25)

    # a_ufloat = ufloat(a, err_a)
    mu_ufloat = ufloat(mu, err_mu)
    print(mu_ufloat)

    # Calcul de la CDA avec incertitudes
    CDA1 = -1/mu_ufloat * math.log(0.5)
    CDA2 = -1/mu_ufloat * math.log(0.25)

    print(CDA1)
    print(CDA2)



    # CDA1 = x_points[indice_CDA1]
    # CDA1_incert = CDA1 * np.sqrt((0.06375/CDA1)**2+(err_mu/mu)**2)
    # CDA2 = x_points[indice_CDA2]
    # CDA2_incert = CDA2 * np.sqrt((0.06375*0.25/CDA1)**2+(err_mu/mu)**2)


    comptes_array_normalise[0] = ufloat(comptes_array_normalise[0].nominal_value, 0)

    actual = unp.nominal_values(comptes_array_normalise)
    expected = exponential(unp.nominal_values(epaisseur_array), mu)

    MSE = mse(actual, expected)

    #### NOTE: mu est ici le coefficient d'atténuation pour le matériau

    fig = plt.errorbar(unp.nominal_values(epaisseur_array), unp.nominal_values(comptes_array_normalise), xerr = unp.std_devs(epaisseur_array), yerr = unp.std_devs(comptes_array_normalise), fmt='o', color=color, capsize=3, markersize=3)
    # fig = plt.scatter(epaisseur_array, comptes_array_normalise, color=color)
    
    plt.plot(x_points, y_points, color="C1", label=f"Lissage exponentiel avec incertitude\nNt/N0=exp(-ux)\n u={mu:.2f}±{err_mu:.5f}\nEQM={'%.2E' % Decimal(MSE)}\nCDA1 = {CDA1.nominal_value:.5f}±{CDA1.std_dev:.5f}\nCDA2 = {CDA2.nominal_value:.5f}±{CDA2.std_dev:.5f}")

    # plt.text(epaisseur_array[int(len(epaisseur_array)/2)-2], moy_comptes_array[int(len(moy_comptes_array)/2)+3],f'{a:.2f}*exp(-{mu:.2f}x)', horizontalalignment='center',
    #  verticalalignment='center')



    # plt.yscale("log")
    # plt.xscale("log")
    plt.xlabel("Épaisseur de filtre [cm]")

    if diviser_par_temps:
        plt.ylabel("Rapport Nt/N0 par seconde")
    else:
        plt.ylabel("Rapport Nt/N0 moyen total")

    plt.title(title+f", {selected_range}")

    return fig


ranges = ["all", "31keV", "12keV", "custom"]

filenames = sets_donnees.Al3_set
path = "./Data/seance3/"

title = "50 kV, 15 uA, Al"

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

############################

filenames = sets_donnees.Cu_set
path = "./Data/seance1/"

title = "50 kV, 15 uA, Cu"
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