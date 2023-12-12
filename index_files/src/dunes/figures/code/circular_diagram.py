import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from color_template_dune import color_dune, color_flux, flow_color
from matplotlib.colors import to_rgba
from matplotlib.sankey import Sankey


def counter_color(limit, color, alpha):
    return to_rgba(color, alpha) if i < limit else color


def add_horizontal_text(ax, color, pos, text):
    offset = 0.1
    ax.text(pos[0], pos[1], text, weight="bold", ha="center", color=color)
    return


color_text = "k"
pos = [2, 1]
alpha = 0.3

limit_topo = 1
limit_wind = 2
limit_flux = 3

for i in range(4):
    figsize = 0.4 * np.array(quarto.regular_fig_size)
    fig, ax = plt.subplots(1, 1, layout="constrained", figsize=figsize)
    ax = plt.gca()
    plt.axis("off")
    sankey = Sankey(ax=ax, gap=0.5, scale=2.0 / pos[0], margin=0.1, shoulder=0.15)

    sankey.add(
        facecolor=counter_color(limit_topo, color_dune, alpha),
        flows=[pos[1], -pos[1]],
        labels=[None, None],
        pathlengths=[0.5, 0.25],
        orientations=[-1, -1],
        connect=(1, 0),
        edgecolor=counter_color(limit_topo, "k", alpha),
    )
    sankey.add(
        facecolor=counter_color(limit_wind, flow_color, alpha),
        flows=[pos[1], -pos[1]],
        labels=[None, None],
        pathlengths=[0.5, 0.25],
        orientations=[0, -1],
        prior=0,
        connect=(1, 0),
        edgecolor=counter_color(limit_wind, "k", alpha),
    )
    sankey.add(
        facecolor=counter_color(limit_flux, color_flux, alpha),
        flows=[pos[1], -pos[1]],
        labels=[None, None],
        pathlengths=[0.12, 0.75],
        orientations=[0, -1],
        prior=1,
        connect=(1, 0),
        edgecolor=counter_color(limit_flux, "k", alpha),
    )

    add_horizontal_text(ax, counter_color(limit_flux, "k", alpha), [-1, -1.7], "Flux")
    add_horizontal_text(ax, counter_color(limit_flux, "k", alpha), [-1, -2.05], "$Q$")
    add_horizontal_text(ax, counter_color(limit_wind, "k", alpha), [1, -1.7], "Wind")
    add_horizontal_text(ax, counter_color(limit_wind, "k", alpha), [1, -2.05], r"$\tau$")
    add_horizontal_text(ax, counter_color(limit_topo, "k", alpha), [0, 0], "Topography")
    add_horizontal_text(ax, counter_color(limit_topo, "k", alpha), [0, -0.35], r"$h$")

    diagrams = sankey.finish()

    # %% saving figure
    figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), i)
    fig.savefig(figname, dpi=1200)

    fig = sg.fromfile(figname)
    newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

    fig.set_size(newsize)
    fig.save(figname)
