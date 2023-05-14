from fonctions import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sets_donnees

def print_somme_nb(filename, path, selected_range="all"):

    filepath = f"{path}{filename}.mca" # Nom du fichier à analyser

    tension, courant, filtre = extraire_params(filename)
    data_array, abscisses_array, live_time, real_time = extraire_data(filepath, w_uncert=True)

    abscisses_array, err = etalonnage_w_stdev(abscisses_array, 1)
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

    j = 0
    en_max = 0
    for i, nb in enumerate(comptes):
        if i < 50:
            continue

        if nb == 0:
            # print(f"Zéro comptes trouvés à i={i}")
            j += 1
        elif nb != 0:
            j = 0

        if j > 3:
            # print(f"Énergie max à {energie[i]}")
            en_max = energie[i]
            break


    moy = ufloat(0, 0)
    nb = ufloat(0, 0)

    for i in range(len(energie)):
        moy += energie[i]*comptes[i]
        nb += comptes[i]

    print(f"\n######### {filename} #########")
    print(f"Erreur sur l'énergie max: {err} keV")
    print(f"Énergie max={en_max:.2f} keV")
    print(f"Énergie moy={(moy/nb):.2f}, {nb:.2f} comptes total")
    print(f"reduc de {( (ufloat(11239, 17.87)-nb)/ufloat(11239, 17.87) * 100):.2f}%")

# filename = "50_15_pas"
# path = "./Data/seance1/"
# print_somme_nb(filename, path)

# filename = "50_15_Al&10"
# print_somme_nb(filename, path)

# filename = "50_15_Cu&1"
# print_somme_nb(filename, path)

# filename = "50_15_W&1"
# print_somme_nb(filename, path)

# filename = "50_15_Mo&1"
# print_somme_nb(filename, path)

# filename = "50_15_Ag&1_W&1"
# print_somme_nb(filename, path)

# filename = "50_15_W&1_Ag&1"
# print_somme_nb(filename, path)

# filename = "50_15_Al&10_W&1"
# print_somme_nb(filename, path)

# filename = "50_15_W&1_Al&10"
# print_somme_nb(filename, path)

# filename = "50_15_Al&10#Mo&1"
# print_somme_nb(filename, path)

# filename = "50_15_Mo&1_Al&10"
# print_somme_nb(filename, path)

# filename = "50_15_Mo&1_Cu&1"
# print_somme_nb(filename, path)

# filename = "50_15_Cu&1_Mo&1"
# print_somme_nb(filename, path)

# filename = "50_15_Cu&1_W&1"
# print_somme_nb(filename, path)

# filename = "50_15_W&1_Cu&1"
# print_somme_nb(filename, path)

######################################################

print("========================= ARGENT =========================")
filename = "50_15_Ag&1"
path = "./Data/seance3/"
print_somme_nb(filename, path)

print("========================= TENSION =========================")
filename = "10_10_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "15_10_pas"
print_somme_nb(filename, path)

filename = "20_10_pas"
print_somme_nb(filename, path)

filename = "25_10_pas"
print_somme_nb(filename, path)

filename = "30_10_pas"
print_somme_nb(filename, path)

filename = "35_10_pas"
print_somme_nb(filename, path)

filename = "40_10_pas"
print_somme_nb(filename, path)

filename = "45_10_pas"
print_somme_nb(filename, path)

filename = "50_10_pas"
print_somme_nb(filename, path)

print("========================= COURANT =========================")
filename = "30_5_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_6,5_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_8_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_9,5_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_11_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_12,5_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_14_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_15,5_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_17_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_18,5_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_20_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

filename = "30_5_pas"
path = "./Data/seance1/"
print_somme_nb(filename, path)

