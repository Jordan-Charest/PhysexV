from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees

def generer_mu(filenames, path, title, selected_range="all", uncertainties=(0,0), set="Al", color="black"):

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

        

        tension_array[i] = tension
        courant_array[i] = courant
        livetime_array[i] = live_time
        realtime_array[i] = real_time
        epaisseur_array[i] = epaisseur
        somme_comptes_array[i] = somme
        incert_array[i] = np.sqrt(somme)

        if i == 1 and set == "Al":
            somme1 = somme
            somme1_incert = np.sqrt(somme)

        if i == 0 and set == "Cu":
            somme1 = somme
            somme1_incert = np.sqrt(somme)


    tension_array = tension_array[1:]
    courant_array = courant_array[1:]
    livetime_array = livetime_array[1:]
    realtime_array = realtime_array[1:]
    epaisseur_array = epaisseur_array[1:]
    somme_comptes_array = somme_comptes_array[1:]
    incert_array = incert_array[1:]


    moy_comptes_array = somme_comptes_array
    # moy_comptes_array = somme_comptes_array / somme1
    # incert_array = moy_comptes_array * np.sqrt((incert_array / somme_comptes_array)**2 + (somme1_incert / somme1)**2)

    # pythasson_array = moy_comptes_array * np.sqrt((pythasson_array/somme_comptes_array)**2 + (pythasson_array[0]/somme0)**2)

    ### AFFICHAGE DU GRAPHIQUE ###
    guess = (somme1, 0)
    a, mu, pcov = exponential_fit(epaisseur_array, moy_comptes_array, guess, yerr=incert_array)
    print(f"a = {a}")
    print(f"mu = {mu}")

    err_mu = np.sqrt(np.diag(pcov))[1]

    return mu, err_mu


###########

filenames = sets_donnees.Al2_set
path = "./Data/seance2/"

title = "Nb comptes en fct de l'épaisseur de filtre, 50 kV, 15 uA, Aluminium"

# uncx = np.zeros(len(filenames))
# uncy = np.ones(len(filenames))
# fig3 = generer_graph(filenames, path, title, selected_range=10, uncertainties=(uncx, uncy), color="blue")

al_mu_all, al_err_mu_all = generer_mu(filenames, path, title, selected_range="all", set="Al", color="C0")

al_mu15, al_err_mu15 = generer_mu(filenames, path, title, selected_range=15, set="Al", color="C0")

al_mu20, al_err_mu20 = generer_mu(filenames, path, title, selected_range=20, set="Al", color="C0")

al_mu30, al_err_mu30 = generer_mu(filenames, path, title, selected_range=30, set="Al", color="C0")

list_al_data = [al_mu15, al_mu20, al_mu30]
list_al_err = [al_err_mu15, al_err_mu20, al_err_mu30]

print(f"AL: {list_al_data}, {list_al_err}")

list_al_data_nist = [21.479, 9.291, 3.046]




# uncx = np.zeros(len(filenames))
# uncy = np.ones(len(filenames))
# fig1 = generer_graph(filenames, path, title, selected_range="31keV", uncertainties=(uncx, uncy), color="purple")

filenames = sets_donnees.Cu_set
path = "./Data/"

title = "Nb comptes en fct de l'épaisseur de filtre, 50 kV, 15 uA, Cuivre"
cu_mu_all, cu_err_mu_all = generer_mu(filenames, path, title, selected_range="all", set="Cu", color="C0")

cu_mu15, cu_err_mu15 = generer_mu(filenames, path, title, selected_range=15, set="Cu", color="C0")

cu_mu20, cu_err_mu20 = generer_mu(filenames, path, title, selected_range=20, set="Cu", color="C0")

cu_mu30, cu_err_mu30 = generer_mu(filenames, path, title, selected_range=30, set="Cu", color="C0")

list_cu_data = [cu_mu15, cu_mu20, cu_mu30]
list_cu_err = [cu_err_mu15, cu_err_mu20, cu_err_mu30]
print(f"CU: {list_cu_data}, {list_cu_err}")

list_cu_data_nist = [661.489, 301.846, 97.548]

fig, (ax1, ax2) = plt.subplots(2)

fig.supxlabel("Énergie [keV]")
fig.supylabel("Coefficient d'atténuation [cm^-1]")

xdata = [1, 4, 7]
xdata_nist = [2, 5, 8]

for i in range(3):
    ax1.errorbar(xdata[i], list_al_data[i], yerr=list_al_err[i], fmt='o', capsize=4, color="C0", label="Données expérimentales")
    ax1.errorbar(xdata_nist[i], list_al_data_nist[i], yerr=0, fmt='o', capsize=4, color="C1", label="Valeurs de référence du NIST")
    ax2.errorbar(xdata[i], list_cu_data[i], yerr=list_cu_err[i], fmt='o', capsize=4, color="C0")
    ax2.errorbar(xdata_nist[i], list_cu_data_nist[i], yerr=0, fmt='o', capsize=4, color="C1")

ax1.get_xaxis().set_visible(False)
ax2.get_xaxis().set_visible(False)

ax1.set_ylabel("Al")
ax2.set_ylabel("Cu")

ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, ncol=2)

plt.legend()
plt.show()