from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees
from decimal import Decimal
from uncertainties import ufloat
from uncertainties import unumpy as unp
# from uncertainties.umath import *

largeur = 20

def generer_graph(filenames, path):

    # Arrays de données pour le calcul de a_tau
    array_dict = [sets_donnees.aluminium, sets_donnees.cuivre, sets_donnees.molybdene, sets_donnees.argent, sets_donnees.tungstene]
    Z_array = [dico["Z"] for dico in array_dict]
    A_array = [dico["A"] for dico in array_dict]
    rho_array = [dico["rho"] for dico in array_dict]

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

    somme_comptes_array = np.zeros(len(filenames))
    somme_comptes_array = unp.uarray(somme_comptes_array, 0)

    for i, filename in enumerate(filenames):

        filepath = f"{path}{filename}.mca" # Nom du fichier à analyser

        # Extraction des données
        try:
            tension, courant, filtre = extraire_params(filename)
            data_array, abscisses_array, live_time, real_time = extraire_data(filepath)
        except:
            print("Erreur de chargement")
            continue

        # Étalonnage des abscisses
        abscisses_array = etalonnage(abscisses_array, 2)


        # Extraction de l'épaisseur du filtre en mm
        if filtre == "pas":
            epaisseur = ufloat(0, 0)
        else:
            index_epp = find_nth_occurrence(filtre, "&", 1)
            epaisseur = ufloat(int(filtre[index_epp+1:])*0.002540, 0.06375*0.002540)


        #Sélection de la bonne plage (31 ± 0,23 keV)
        indice_centre = find_nearest(abscisses_array, 31)
        indice_min = int(indice_centre-largeur)
        indice_max = int(indice_centre+largeur)
        energie_max = abscisses_array[indice_max]
        plage = energie_max - 31

        #Réduction des données en y à la plage réduite
        data_reduit = data_array[indice_min:indice_max]

        # Somme des comptes sur la plage normalisée par le live time
        somme = ufloat(np.sum(data_reduit), np.sqrt(np.sum(data_reduit))) / live_time

        # Extraction des données dans l'array approprié
        tension_array[i] = ufloat(tension.replace(",", "."), 0)
        courant_array[i] = ufloat(courant.replace(",", "."), 0)
        livetime_array[i] = live_time
        realtime_array[i] = real_time
        epaisseur_array[i] = epaisseur
        somme_comptes_array[i] = somme

        # Si i=0: sans filtre, donc la somme p/r à laquelle on va normaliser
        if i == 0:
            somme0 = somme

    ### FIN DE LA BOUCLE SUR LES FICHIERS

    # Retire les données sans filtre
    tension_array = tension_array[1:]
    courant_array = courant_array[1:]
    livetime_array = livetime_array[1:]
    realtime_array = realtime_array[1:]
    epaisseur_array = epaisseur_array[1:]
    somme_comptes_array = somme_comptes_array[1:]
    # incert_y_array = incert_y_array[1:]

    # On crée un array pour les sections efficaces
    atau_array = np.zeros(len(array_dict))
    atau_array = unp.uarray(atau_array, 0)

    # On calcule a_tau pour chaque point
    for i in range(len(array_dict)):
        N0_tuple = somme0 # Valeur, incert
        Nt_tuple = somme_comptes_array[i]
        atau = a_tau(N0_tuple, Nt_tuple, A_array[i], rho_array[i], epaisseur_array[i])
        atau_array[i] = atau



    ### DONNÉES FINALES

    donnees_atau = [unp.nominal_values(atau_array[:-1]), unp.std_devs(atau_array[:-1])]
    donnees_Z = Z_array[:-1]

    # donnees_atau[0] = donnees_atau[0][:-1]
    # donnees_Z = donnees_Z[:-1]

    print(donnees_atau[0])
    print(donnees_atau[1])
    print(donnees_Z)

    # On effectue un fit linéaire en log-log

    atau_log10 = np.log10(donnees_atau[0])
    Z_log10 = np.log10(donnees_Z)

    a, b, pcov = linear_fit(Z_log10, atau_log10, (4, 0), yerr=unp.std_devs(atau_array[:-1])) # Changer yerr éventuellement
    err_a = np.sqrt(pcov)[0][0]
    popt = [a, b]

    # Calcul R2:
    residuals = atau_log10 - linear(Z_log10, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((atau_log10-np.mean(atau_log10))**2)
    R2 = 1 - (ss_res/ss_tot)


    # Produit le graph

    x_steps = np.arange(Z_log10[0], Z_log10[-1]+0.02, (Z_log10[-1]-Z_log10[0])/50)

    fig = plt.errorbar(donnees_Z, donnees_atau[0], yerr=donnees_atau[1], fmt='o', capsize=4, label=f"Section efficace a_tau")

    plt.plot(10**x_steps, 10**(x_steps*a + b), color="C1", label=f"Lissage linéaire sur les données en log-log\n pente={a:.2f}±{err_a:.2f}, R^2={R2:.3f}")

    title = f"Section efficace a_tau en fonction du Z du matériau, 50 kV, 15 uA, plage de 31±{plage:.2f} keV"

    plt.title(title)

    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Z [éch. log]")
    plt.ylabel("atau [éch. log]")

    return fig


filenames = sets_donnees.materiaux_set2
path = "./Data/Filtres/"

fig1 = generer_graph(filenames, path)

plt.legend()
plt.show()