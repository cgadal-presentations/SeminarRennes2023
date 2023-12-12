import os
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from color_template_dune import alpha_dune, color_dune, flow_color
from pydune.math import tand
from pydune.physics.dune.bedinstability_1D import temporal_growth_rate, temporal_velocity

A0 = 3.5
B0 = 2
mu = tand(35)
r = 1.26
Q = 38
Lsat = 2.8

linewidth_small = 0.6

k = np.linspace(0, 0.4, 1000)
sigma_t = temporal_growth_rate(k, A0, B0, mu, r)
celerity_t = temporal_velocity(k, A0, B0, mu, r)
kc = (B0 - (1 / mu) * (1 / r**2)) / A0


def topo(x, amp, kg):
    h0 = 0.1
    return amp * (np.sin(kg * x - np.pi / 2) + 1) + h0


############################ Figure sketch
for i in range(5):
    figsize = (0.75 * quarto.regular_fig_width, 0.4 * quarto.regular_fig_height)
    fig, ax = plt.subplots(1, 1, layout="constrained", figsize=figsize)

    x = np.linspace(0, 1, 200)
    amp = 0.2 if i > 0 else 0
    kg = 2 * np.pi / 1
    yg = topo(x, amp, kg)

    dTau_coeff = 0.25
    dLsat = 0.15
    hzone = 0.85

    plt.setp(ax.get_yticklabels(), visible=False)
    ax.set_yticks([])
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_xticks([])
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.plot(x, yg, color=color_dune)
    plt.fill_between(x, 0, yg, color=color_dune, alpha=alpha_dune)

    if i > 0:
        ax.vlines(0.5, 0, 0.59, transform=ax.transAxes, linewidth=linewidth_small, linestyle="--", color="k")
        ax.text(0.5 - 0.03, 0.625, r"$h_{\text{max}}$")
        ax.vlines(
            0.5,
            0.68,
            hzone if i > 3 else 1,
            transform=ax.transAxes,
            linewidth=linewidth_small,
            linestyle="--",
            color="k",
        )
        #
        xy = np.array([0, 0.05])
        xytext = np.array([0.63, 0.05])
        ax.annotate("", xy=xy, xytext=xytext, arrowprops=dict(arrowstyle="-|>", color="black", shrinkA=0, shrinkB=0))
        xy = np.array([0.67, 0.05])
        xytext = np.array([1, 0.05])
        ax.annotate("", xy=xy, xytext=xytext, arrowprops=dict(arrowstyle="<|-", color="black", shrinkA=0, shrinkB=0))
        ax.text(0.65, 0.035, r"$\lambda$", ha="center")

        if i > 1:
            xTau = 0.5 - dTau_coeff * 2 * np.pi / kg
            ax.vlines(xTau, 0.2, 0.59, transform=ax.transAxes, linewidth=linewidth_small, linestyle=":", color="k")
            ax.text(xTau - 0.04, 0.625, r"$\tau_{\text{max}}$")
            ax.vlines(
                xTau,
                0.68,
                hzone if i > 3 else 1,
                transform=ax.transAxes,
                linewidth=linewidth_small,
                linestyle=":",
                color="k",
            )

            xy = np.array([0.5, 0.2])
            xytext = np.array([xTau, 0.2])
            ax.annotate(
                "", xy=xy, xytext=xytext, arrowprops=dict(arrowstyle="<|-", color="black", shrinkA=0, shrinkB=0)
            )
            ax.text(xTau + 0.4 * (0.5 - xTau), 0.235, r"$\propto \lambda$")

            if i > 2:
                xQ = xTau + dLsat
                ax.vlines(
                    xQ,
                    topo(xQ, amp, kg),
                    0.715,
                    transform=ax.transAxes,
                    linewidth=linewidth_small,
                    color=("tab:red" if i > 3 else "k"),
                )
                ax.text(xQ - 0.04, 0.74, r"$Q_{\text{max}}$")
                ax.vlines(
                    xQ,
                    0.78,
                    1,
                    transform=ax.transAxes,
                    linewidth=linewidth_small,
                    color=("tab:red" if i > 3 else "k"),
                )

                xy = np.array([xTau, topo(xQ, amp, kg)])
                xytext = np.array([xQ, topo(xQ, amp, kg)])
                ax.annotate(
                    "", xy=xy, xytext=xytext, arrowprops=dict(arrowstyle="<|-", color="black", shrinkA=0, shrinkB=0)
                )
                ax.text(xTau + 0.33 * (xQ - xTau), 0.51, r"$L_{\text{sat}}$")

                if i > 3:
                    ax.hlines(hzone, 0, 1, transform=ax.transAxes, linewidth=0.8, color="k")
                    ax.text(0.15, 0.885, "Erosion")
                    ax.text(0.62, 0.89, "Deposition")

    arrow = mpatches.FancyArrowPatch((0.01, 0.6), (0.13, 0.6), mutation_scale=20, facecolor=flow_color)
    ax.add_patch(arrow)
    ax.text(0.013, 0.67, r"Wind")

    # %% saving figure
    figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), i)
    fig.savefig(figname, dpi=1200)

    fig = sg.fromfile(figname)
    newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

    fig.set_size(newsize)
    fig.save(figname)
