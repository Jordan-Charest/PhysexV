import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt


### PARAMÈTRES ###

filename = "50_15_Al10" # Fichier à charger
diviser_par_temps = True # Diviser le nb de comptes par le live time



### EXTRACTION DES DONNÉES ###

filepath = f"./Data/{filename}.mca" # Nom du fichier à analyser

# Indices des underscore pour titres
index_first_under = filename.index("_")
index_second_under = filename.index("_",index_first_under+1,len(filename)-1)

# Extraction des valeurs de tension, du courant et des filtres
tension = filename[:index_first_under]
courant = filename[index_first_under+1:index_second_under]
filtre = filename[index_second_under+1:]
print(f"Analyse du filtre [{filtre}] à {tension} V et {courant} A")

# Ouverture du fichier
f = open(filepath, "r")

# np.array pour stocker le nombre de comptes
data_array = np.empty(4096)

# np.array d'abscisses (canaux)
abscisses_array = np.arange(1, 4097, 1)

# Initialisation des variables de temps
live_time = 0
real_time = 0

# Extraction des données
for index, line in enumerate(f):
    if index == 7: # Live time
        live_time = float(line[-10:])

    if index == 8: # Real time
        real_time = float(line[-10:])

    if index > 4112: # Fin des données
        break

    if index > 16: # Début des données
        data_array[index-17] = int(line)

f.close()



### TRAITEMENT DES DONNÉES ###

if diviser_par_temps: # Diviser par le live time
    data_array = data_array / live_time

# Étalonnage de l'axe des canaux (abscisses) en énergie

def etalonnage(canaux):

    params = linregress([0, 921.8, 1173.19, 3895.6], [0, 13.95, 17.24, 59.54])
    # Ce qui nous donne un std error, nice

    return params.slope * canaux + params.intercept

abscisses_array = etalonnage(abscisses_array)



### AFFICHAGE DU GRAPHIQUE ###

fig = plt.plot(abscisses_array, data_array)
plt.yscale("log")
plt.xlabel("Énergie [keV]")

if diviser_par_temps:
    plt.ylabel("Nombre de comptes par seconde [log]")
else:
    plt.ylabel("Nombre de comptes total [log]")

plt.title(f"Nombre de comptes en fonction de l'énergie, {tension} V, {courant} A, filtres {filtre}")

plt.show()


