from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees


def print_somme_nb(filename, path, selected_range="all"):

    filepath = f"{path}{filename}.mca" # Nom du fichier à analyser

    tension, courant, filtre = extraire_params(filename)
    data_array, abscisses_array, live_time, real_time = extraire_data(filepath, w_uncert=True)

    abscisses_array = etalonnage(abscisses_array, 1)
    # print(f"live time: {live_time}")
    # print(f"real time: {real_time}")

    # if filtre == "pas" or filtre == "pos" or filtre == "pus":
    #     epaisseur = 0
    # else:
    #     index_epp = find_nth_occurrence(filtre, "&", 1)
    #     epaisseur = int(filtre[index_epp+1:].replace(",", ".")) * 0.00254

    # Selection de la bonne plage

    if selected_range == "all":
        indice_min = find_nearest(abscisses_array, 5)
        indice_max = find_nearest(abscisses_array, 50)


    energie = abscisses_array[indice_min: indice_max]
    comptes = data_array[indice_min: indice_max]/live_time

    moy = ufloat(0, 0)
    nb = ufloat(0, 0)

    for i in range(len(energie)):
        moy += energie[i]*comptes[i]
        nb += comptes[i]

    # print(f"Pour {filename}: moy={(moy/nb):.2f}, {nb:.2f} comptes\n")
    print(f"Pour {filename}: reduc de {( (ufloat(11239, 17.87)-nb)/ufloat(11239, 17.87) * 100):.2f}%")

filename = "50_15_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "50_15_Al&10"
print_somme_nb(filename, path)

filename = "50_15_Cu&1"
print_somme_nb(filename, path)

filename = "50_15_W&1"
print_somme_nb(filename, path)

filename = "50_15_Mo&1"
print_somme_nb(filename, path)

filename = "50_15_Ag&1_W&1"
print_somme_nb(filename, path)

filename = "50_15_W&1_Ag&1"
print_somme_nb(filename, path)

filename = "50_15_Al&10_W&1"
print_somme_nb(filename, path)

filename = "50_15_W&1_Al&10"
print_somme_nb(filename, path)

filename = "50_15_Al&10#Mo&1"
print_somme_nb(filename, path)

filename = "50_15_Mo&1_Al&10"
print_somme_nb(filename, path)

filename = "50_15_Mo&1_Cu&1"
print_somme_nb(filename, path)

filename = "50_15_Cu&1_Mo&1"
print_somme_nb(filename, path)

filename = "50_15_Cu&1_W&1"
print_somme_nb(filename, path)

filename = "50_15_W&1_Cu&1"
print_somme_nb(filename, path)