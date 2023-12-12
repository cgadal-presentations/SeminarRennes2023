import os
import sys

import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from matplotlib import pyplot as plt
from uncertainties import unumpy as unp

# #### parameters
manip_density = [1, 2, 4]
manip_particle = [5, 6, 8]

rho_c = {1: 1078.5, 2: 1024.4, 4: 1049.4, 5: 1047.5, 6: 1076.6, 8: 1017.9, 9: 1081.7}

rho_a = {1: 998.4, 2: 1007.5, 4: 1000.4, 5: 1000.4, 6: 1000.4, 8: 1000.1, 9: 1030.5}

phi = {1: 0, 2: 0, 3: 0, 5: 3, 6: 3, 8: 3, 9: 3}

zref = {
    1: 0.64,
    2: 0.57,
    4: 0.71,
    #
    5: 0.07,
    6: 0.05,
    8: 0.4,
    9: 0.4,
}

rho_p = 1050  # kg/m3
d_p = 500e-6  # m
Q = 1.2  # L/s
W = 21.5  # cm

labels = ["settling", "neutral", "buoyant"]

# %%
path_data_phi = "/media/cyril/IMFT data/202211_VISITE_JEAN/cyril/analysis_images/compute_average_profiles"

path_data_u = "/media/cyril/IMFT data/202211_VISITE_JEAN/cyril/analysis_UDV/runs"

manip_number = 5
data_phi = np.load(os.path.join(path_data_phi, f"manip{manip_number:0d}", "vertical_profile.npy"), allow_pickle=True)
data_u = np.load(os.path.join(path_data_u, f"manip{manip_number:0d}", "velocity_av.npy"), allow_pickle=True).item()

z_u, u = data_u["z_interp"], data_u["U_av"][0]
z_phi, phi = unp.nominal_values(data_phi)

u[z_u > 7.5] = np.nan
phi[z_phi > 7.5] = np.nan


figsize = (0.4 * quarto.regular_fig_width, 0.37 * quarto.regular_fig_height)
fig, axarr = plt.subplots(1, 2, layout="constrained", sharey=True, figsize=figsize)

axarr[0].plot(-u, z_u)
axarr[1].plot(phi, z_phi)

axarr[0].set_ylim(0, 11)
axarr[0].set_xlim(left=0)
axarr[1].set_xlim(0, 3)

axarr[0].set_ylabel("height [cm]")
axarr[0].set_xlabel("$u$ [cm/s]")
axarr[1].set_xlabel(r"$\phi~[\%]$")

axarr[0].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
axarr[0].xaxis.set_label_position("top")
axarr[1].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
axarr[1].xaxis.set_label_position("top")


figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=600)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
