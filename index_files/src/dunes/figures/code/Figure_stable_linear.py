import glob

import matplotlib.gridspec as gridspec
import matplotlib.patches as ptch
import matplotlib.pyplot as plt
import numpy as np
from General.Math import Rotation_matrix, cosd, sind
from Instability.Functions import Resultant_flux
from PIL import Image
from pydune.math import cosd, sind
from Wind_data.Wind_treatment import PDF_flux, Wind_to_flux, flux_rose


def North_arrow(fig, ax, center, length, length_small, width, radius, theta=0, color="k"):
    y_start = radius + length - length_small
    arrow = np.array([[0, y_start], [width / 2, radius], [0, radius + length], [-width / 2, radius], [0, y_start]])
    barycentre = np.sum(arrow, axis=0) / arrow.shape[0]
    # arrow = np.dot(Rotation_matrix(theta), (arrow-barycentre).T).T + barycentre
    arrow = np.dot(Rotation_matrix(theta + 180), arrow.T).T
    arrow = arrow + np.array(center)
    ax.add_patch(ptch.Polygon(arrow, color=color))
    t = ax.text(center[0], center[1], "N")
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    bb = t.get_window_extent(renderer=r).transformed(inv)
    t.set_visible(False)
    width_t = bb.width
    height_t = bb.height
    t = ax.text(center[0], center[1] - height_t / 2, r"\textbf{N}", color=color, ha="center")


fig = plt.figure(figsize=(0.8 * Beamer.fig_width, 0.8 * 0.37 * Beamer.fig_width))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 1])
gs.update(left=0.002, right=0.996, bottom=0.002, top=0.998, hspace=0.05)

ax0 = plt.subplot(gs[0])
ax0.set_xticks([])
ax0.set_yticks([])
path_gen = "/home/gadal/Documents/Work/Research/DUNE/PhD_Parts/Thèse/Chapitres/chapter4/figures/Figure_stable_linear"
path = path_gen + "/Niger"

# ####### image
img = np.array(Image.open(glob.glob(path + "/*.png")[0]))
plt.imshow(img[440:900, 160:])

# ####### échelle
backgrnd = ptch.Rectangle((0, 0), width=0.12, height=0.2, transform=ax0.transAxes, color="w", alpha=0.6)
ax0.add_patch(backgrnd)
plt.plot([30, 30 + 205], [445, 445], linewidth=Beamer.linewidth_barscale, color="k")
name = r"$" + str(250) + r"~\textrm{m}$"
plt.text((205 + 2 * 30) / 2, 430, name, color="k", ha="center")

# ####### Flux roses

data = np.loadtxt(path + "/wind_data.txt")
direction = data[:, 0] + 190
speed = data[:, 1]
qs, direction = Wind_to_flux(direction, speed, 180e-6)
pdfQ, Angles = PDF_flux(direction, qs)

# ######## Wind rose
#
RDD = Resultant_flux(Angles, pdfQ)[0]
#
anchor = [0.69, 0.05, 0.5, 0.5]
subax = ax0.inset_axes(bounds=anchor, transform=ax0.transAxes)
a = flux_rose(
    Angles,
    pdfQ,
    place=subax.get_position(),
    fig=fig,
    opening=1,
    nsector=20,
    edgecolor="k",
    linewidth=0.5,
    cmap=Beamer.flux_cmap,
)
a.annotate(
    "",
    (RDD * np.pi / 180, 0),
    (RDD * np.pi / 180, 0.85 * a.get_rmax()),
    arrowprops=dict(arrowstyle="<|-", shrinkA=0, shrinkB=0, color=Beamer.color_resultant_flux),
)
a.grid(linewidth=0.4, color=Beamer.wind_grid_color, linestyle="--")
a.set_axisbelow(True)
a.patch.set_alpha(0.4)
a.set_xticklabels([])
subax.remove()

ax1 = plt.subplot(gs[1])
ax1.set_xticks([])
ax1.set_yticks([])
img = np.array(Image.open(path_gen + "/DUN00143_t0.png").convert("RGBA"))
plt.imshow(img[4:200, 13:1450])

# ###########
radius = 100
xy = np.array([1310, 130])
xytext = xy + radius * np.array([cosd(-40), sind(-40)])
xytext2 = xy + radius * np.array([cosd(80), sind(80)]) / 2
xytext3 = xy + ((xytext - xy) + (xytext2 - xy))


ax1.annotate(
    "",
    xy,
    xytext,
    arrowprops=dict(arrowstyle="<|-", shrinkA=0, shrinkB=0, color=Beamer.flux_cmap(np.linspace(0, 1, 1)[0])),
)
ax1.annotate(
    "",
    xy,
    xytext2,
    arrowprops=dict(arrowstyle="<|-", shrinkA=0, shrinkB=0, color=Beamer.flux_cmap(np.linspace(0, 1, 1)[0])),
)
ax1.annotate(
    "", xy, xytext3, arrowprops=dict(arrowstyle="<|-", shrinkA=0, shrinkB=0, color=Beamer.color_resultant_flux)
)

plt.plot([107], [100], ".")

# ####### échelle
backgrnd = ptch.Rectangle((0, 0), width=0.13, height=0.33, transform=ax1.transAxes, color="w", alpha=0.6)
ax1.add_patch(backgrnd)
plt.plot([15, 15 + 150], [185, 185], linewidth=Beamer.linewidth_barscale, color="k")
name = r"$" + str(150) + r"~l_{0}$"
plt.text((150 + 2 * 15) / 2, 170, name, color="k", ha="center")


plt.savefig("../Figures/Figure_stable_linear.pdf", dpi=600)
