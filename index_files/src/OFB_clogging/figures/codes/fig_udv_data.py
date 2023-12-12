import os
import sys

import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from matplotlib import pyplot as plt

sys.path.append("/media/cyril/IMFTPostDoc/OFB_Colmatage/processing/UDV")

path_data = "/media/cyril/IMFTPostDoc/OFB_Colmatage/processing/UDV/release01"

data = np.load(os.path.join(path_data, "data.npy"), allow_pickle=True).item()

z_interp = data["z_interp"]
u, v = data["U_av"]

probe = 13
z, phi = data[probe]["z"], data[probe]["phi_inferred"][1]

figsize = (0.5 * quarto.regular_fig_width, 0.5 * quarto.regular_fig_height)
fig, axarr = plt.subplots(1, 2, layout="constrained", figsize=figsize, sharey=True)

axarr[0].plot(-u, z_interp - 0.73)
# axarr[0].plot(v, z_interp - 0.73, label="v")

axarr[1].plot(phi, z - 0.73)

axarr[0].set(xlabel=r"Velocity [cm/s]", xlim=(0.03, 0.055))
axarr[0].set(ylabel=r"Height [cm]", ylim=(0, 6))
axarr[1].set(xlabel=r"$\phi~[\%]$", xlim=(0, 0.2))


# %% saving figure
figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=1200)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
