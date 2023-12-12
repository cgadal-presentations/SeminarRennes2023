import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from netCDF4 import Dataset

# %% Load data
path_data = "/home/cyril/Documents/Work/Research/Projects/Turbity_Currents/palagram_monograph/data/output_data"
list_runs = sorted(glob.glob(os.path.join(path_data, "*.nc")))
list_fitresults = sorted(glob.glob(os.path.join(path_data, "fitresult*")))

datasets = np.array([Dataset(run) for run in list_runs])

# %% mask data
particles = np.array([d.particle_type for d in datasets]).T
mask = particles != "saline water"

zorder_setups = {
    "Cyril Gadal": -9,
    "Marie Rastello": -10,
    "Jean Schneider": -6,
    "Julien Chauchat": -7,
    "Marie Rastello - Cyril Gadal": -10,
}

# %% figure

figsize = (0.6 * quarto.regular_fig_width, 0.81 * quarto.regular_fig_height)
fig, ax = plt.subplots(1, 1, figsize=figsize, layout="constrained")
# All runs

# 90, 0, 109, 118, 21
# selected run, non-dimensional
# i_runs = [0, 158, 94, 197, 210, 202]
i_runs = [90, 0, 109, 118, 94, 210, 202]

for i in i_runs:
    d = datasets[list_runs.index(os.path.join(path_data, f"run_{i:03d}.nc"))]
    # print(d.author)
    #
    t_ad = d.variables["t0"][:].data
    x_ad = d.variables["L0"][:].data
    #
    x_axis = d.variables["t"][:].data / t_ad
    y_axis = d.variables["x_front"][:].data / x_ad
    ax.plot(x_axis, y_axis)

# annotations
ax.annotate(
    r"",
    xy=(19.3, 10.3),
    xycoords="data",
    xytext=(9, 12),
    textcoords="data",
    arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0, connectionstyle="angle3,angleA=0,angleB=120"),
)
ax.text(2, 13, r"global slope $\nearrow$")

ax.annotate(
    r"",
    xy=(29.9, 13.45),
    xycoords="data",
    xytext=(36.6, 8.16),
    textcoords="data",
    arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0, connectionstyle="angle3,angleA=-90,angleB=-20"),
)
ax.text(35.2, 8.6, r"settling $\nearrow$", ha="right", va="center")

ax.set_ylabel(r"Front position, $x_{\rm f}/l_{0}$")
ax.set_xlabel(r"Time, $t/t_{0}$")
ax.set_xlim(-1, 40)
ax.set_ylim(-0.5, 17)

# axarr[0].axis("off")
# leg = axarr[0].legend(handles=tp.legend_datasets, ncol=2, title="Datasets", loc="upper center", borderaxespad=0)

fig.align_labels()

# %% saving figure
figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=600)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
