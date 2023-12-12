import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg


def growth_rate(x, amp, sigma, a):
    return amp * np.exp(sigma * (x)) + a


def topo_space(x, amp, lamb, Lambda, shift=0):
    return amp * np.sin(2 * np.pi * (x + shift) / lamb) * np.exp((x + shift) / Lambda)


matplotlib.rcParams["font.size"] = 8  # default 12


path_general = "/media/cyril/Backup_dune/DUNE/PhD_Parts/Part3_Field_instability_Ewing/Paper/Figures/Figure_terrain"

path_transect = path_general + "/" + "Data/Data_2115.npy"
Data_transect = np.load(path_transect, allow_pickle=True).item()
Dates = sorted(Data_transect.keys())[:4]
date = "2007"
pad_end = 0

Starting_point = np.array(
    [Data_transect[date]["xlim"][0] + np.isnan(Data_transect[date]["transect"][:2500]).sum() for date in Dates]
).min()  # General coordinates
Ending_point = np.array(
    [
        Data_transect[date]["xlim"][1] + np.isnan(Data_transect[date]["transect"][:2500]).sum() + pad_end
        for date in Dates
    ]
).max()  # General coordinates
Profiles = {
    date: Data_transect[date]["filtered"][
        Starting_point - np.isnan(Data_transect[date]["transect"][:2500]).sum() : Ending_point
        - np.isnan(Data_transect[date]["transect"][:2500]).sum()
    ]
    if (Starting_point - np.isnan(Data_transect[date]["transect"][:2500]).sum()) > 0
    else np.nan
    for date in Dates
}

x_profiles = np.linspace(0, Profiles[date].size, Profiles[date].size)

x_fit_top = Data_transect[date]["Abscisses"][1::2] - (
    Starting_point - np.isnan(Data_transect[date]["transect"][:2500]).sum()
)
y_fit_top = Data_transect[date]["Ordonnees"][1::2]
x_fit_bot = Data_transect[date]["Abscisses"][0::2] - (
    Starting_point - np.isnan(Data_transect[date]["transect"][:2500]).sum()
)
y_fit_bot = -Data_transect[date]["Ordonnees"][0::2]
x_fit_top = np.insert(x_fit_top, 0, x_fit_top[0] - (x_fit_top[1] - x_fit_top[0]))  # padding zeros at the beginning
x_fit_top = np.insert(x_fit_top, 0, x_fit_top[0] - (x_fit_top[1] - x_fit_top[0]))  # padding zeros at the beginning
x_fit_bot = np.insert(x_fit_bot, 0, x_fit_bot[0] - (x_fit_bot[1] - x_fit_bot[0]))  # padding zeros at the beginning
x_fit_bot = np.insert(x_fit_bot, 0, x_fit_bot[0] - (x_fit_bot[1] - x_fit_bot[0]))  # padding zeros at the beginning
y_fit_top = np.insert(y_fit_top, 0, 0)
y_fit_top = np.insert(y_fit_top, 0, 0)
y_fit_bot = np.insert(y_fit_bot, 0, 0)
y_fit_bot = np.insert(y_fit_bot, 0, 0)

xtest = np.arange(Profiles[date].size)
p_test = [0.8, 120, 170, -450]


# #### Figure
figsize = (0.7 * quarto.regular_fig_width, 0.45 * quarto.regular_fig_height)
fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=figsize)

ax.plot(xtest + 160, topo_space(xtest, *p_test), label=r"theory")
ax.plot(x_profiles + 160, Profiles[date], label=r"data")
ax.set_ylabel("Bed elevation [m]")
ax.set_xlabel("Distance along profile [m]")
ax.set_xlim([x_profiles.min() + 160, x_profiles.max() + 160])
ax.set_yticks([-4, 0, 4])
ax.legend(loc="lower left")

# %% saving figure
figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=1200)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
