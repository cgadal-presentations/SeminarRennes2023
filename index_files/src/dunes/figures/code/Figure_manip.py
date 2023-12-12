import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg

dashes = (2, 1)


def length_arrow(
    points,
    bar_shift,
    label,
    color="k",
    va="center",
    ha="center",
    relative_arrow_pos=0.75,
    relative_text_pos=1.25,
    arrowstyle="<|-|>",
):
    points_arrow = points + relative_arrow_pos * bar_shift
    plt.annotate(
        "",
        xy=points_arrow[0],
        xycoords="data",
        xytext=points_arrow[1],
        textcoords="data",
        arrowprops=dict(arrowstyle=arrowstyle, shrinkA=0, shrinkB=0, color=color),
    )

    # for tip in points :
    #     plt.plot( [ tip, tip + bar_shift ], color = 'k', linestyle = '--', dashes = dashes )

    textpos = np.mean(points, axis=0) + relative_text_pos * bar_shift
    plt.text(textpos[0], textpos[1], label, color=color, va=va, ha=ha)


img = np.array(plt.imread("src/Photo_manip.png"))

figsize = (0.725 * quarto.regular_fig_width, 0.6 * quarto.regular_fig_height)
fig, ax = plt.subplots(1, 1, figsize=figsize)

plt.subplots_adjust(0, 0, 1, 1)

ax.set_axis_off()
ax.imshow(img)


############## Annotation
color_annotation = "crimson"
color_legend = "w"

xy = np.array([488, 1089])
xytext = np.array([2174, 829])
points = np.array([xy, xytext])
perp = -np.array([xytext[1] - xy[1], xy[0] - xytext[1]]) / np.linalg.norm(points)
length_arrow(points, perp, r"$2$ m", color=color_annotation, relative_arrow_pos=0, relative_text_pos=300)
#
xy = np.array([418, 941])
xytext = np.array([391, 363])
points = np.array([xy, xytext])
perp = -np.array([xytext[1] - xy[1], xy[0] - xytext[1]]) / np.linalg.norm(points)
length_arrow(points, perp, r"$30$ cm", color=color_annotation, relative_arrow_pos=0, relative_text_pos=-200)
#
xy = np.array([1480, 730])
xytext = np.array([1143, 650])
points = np.array([xy, xytext])
perp = -np.array([xytext[1] - xy[1], xy[0] - xytext[1]]) / np.linalg.norm(points)
length_arrow(points, perp, r"$1$ m", color=color_annotation, relative_arrow_pos=0, relative_text_pos=-120)
#
xy = np.array([1465, 639])
xytext = np.array([1562, 61])
ax.annotate(
    r"Rotating plate",
    xy=xy,
    xytext=xytext,
    color=color_legend,
    arrowprops=dict(arrowstyle="->", color=color_legend, connectionstyle="angle, angleA = 180, angleB = 90, rad = 0"),
)
#
xy = np.array([1042, 770])
xytext = np.array([1902, 1011])
ax.annotate(
    r"Water tank", xy=xy, xytext=xytext, color=color_legend, arrowprops=dict(arrowstyle="->", color=color_legend)
)

# %% saving figure
figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=1200)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
