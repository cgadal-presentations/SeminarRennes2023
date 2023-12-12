import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from color_template_dune import flow_color
from pydune.math import cosd, sind

path_gen = (
    "/media/cyril/Backup_dune/DUNE/PhD_Parts/Part1_multidirectional/IPGP_intership/Experimental/Lit_plat/SerieV240"
)

linewidth_barscale = 2
figsize = (1.2 * quarto.regular_fig_width, 1.2 * 0.325 * quarto.regular_fig_height)

# Figure all patterns
theta = [0, 50, 70, 130, 180]
im_numbers = [50, 148, 151, 151, 151]
fig, axarr = plt.subplots(1, 5, figsize=figsize, layout="constrained")

for i, (ax, th, im_number) in enumerate(zip(axarr.flatten(), theta, im_numbers)):
    ax.set_xticks([])
    ax.set_yticks([])
    #
    try:
        image = os.path.join(path_gen, "Theta" + str(th), "Img/image" + f"{im_number:05d}" + ".BMP")
        im = np.array(plt.imread(image))
    except FileNotFoundError:
        image = os.path.join(path_gen, "Theta" + str(th), "Img2/image" + f"{im_number:05d}" + ".BMP")
        im = np.array(plt.imread(image))
    #
    im = im[::-1, ::-1]
    p1 = ax.imshow(im[278:890, 570:1200], animated=True)
    (p2,) = ax.plot([380, 600], [580, 580], color="k", linewidth=linewidth_barscale)
    p3 = ax.text((380 + 600) / 2, 560, r"$10$ cm", color="k", ha="center")
    # vents
    length = 200
    xy = np.array([50, 200])
    #
    th_arrow = 90 if th == 70 else th
    dx = int(round(cosd(th_arrow / 2) * length))
    dy = int(round(sind(th_arrow / 2) * length))
    xytext = xy + np.array([dx, dy])
    p4 = ax.annotate(
        "", xy=xy, xytext=xytext, arrowprops=dict(arrowstyle="<|-", color=flow_color, shrinkA=0, shrinkB=0)
    )
    #
    dx = int(round(cosd(-th / 2) * length))
    dy = int(round(sind(-th / 2) * length))
    xytext = xy + np.array([dx, dy])
    # plt.annotate("", xy=xy, xytext=xytext,arrowprops=dict(arrowstyle="<|-", color = 'white', shrinkA = 0, shrinkB = 0))
    p5 = ax.annotate(
        "", xy=xy, xytext=xytext, arrowprops=dict(arrowstyle="<|-", color=flow_color, shrinkA=0, shrinkB=0)
    )

# %% saving figure
figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=1200)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
