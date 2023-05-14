import numpy as np
from sets_donnees import *
import matplotlib.pyplot as plt
from fonctions import *


def graph(x, y, yerr, title, xlabel, ylabel):

    fig = plt.errorbar(x, y, yerr = yerr, fmt='o', capsize=3, markersize=3)

    if title == "emax en fct de la tension":

        a, b, pcov = linear_fit(x, y, (1, 0), yerr=[0.0002 for i in x]) # Changer yerr éventuellement
        popt = [a, b]

        x = np.array(x)
        y = np.array(y)
        # Calcul R2:
        residuals = y - linear(x, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y-np.mean(y))**2)
        R2 = 1 - (ss_res/ss_tot)

        plt.plot(x, linear(x, *popt), label=f"Lissage linéaire y={a:.3f}x+{b:.3f}\nR^2={R2:.4f}")

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)



    plt.title(title)

    return fig

x = valeurs_courant
y = e_moy_courant
yerr = [0.02 for i in valeurs_courant]
title = "emoy en fct du courant"
xlabel = "Courant [uA]"
ylabel = "Énergie moyenne [keV]"

fig1 = graph(x, y, yerr, title, xlabel, ylabel)

plt.show()
plt.close()

x = valeurs_courant
y = e_max_courant
yerr = [0.0002 for i in valeurs_courant]
title = "emax en fct du courant"
xlabel = "Courant [uA]"
ylabel = "Énergie maximale [keV]"

fig2 = graph(x, y, yerr, title, xlabel, ylabel)

plt.show()
plt.close()

x = valeurs_tension
y = e_moy_tension
yerr = [0.02 for i in valeurs_tension]
title = "emoy en fct de la tension"
xlabel = "Tension [kV]"
ylabel = "Énergie moyenne [keV]"

fig3 = graph(x, y, yerr, title, xlabel, ylabel)

plt.show()
plt.close()

x = valeurs_tension
y = e_max_tension
yerr = [0.0002 for i in valeurs_tension]
title = "emax en fct de la tension"
xlabel = "Tension [kV]"
ylabel = "Énergie maximale [keV]"

fig4 = graph(x, y, yerr, title, xlabel, ylabel)

plt.legend()
plt.show()
plt.close()