import glob
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from color_template_dune import color_BI, color_F, flux_cmap, linewidth_barscale, lw_arrow
from matplotlib.patches import FancyArrowPatch, PathPatch, Polygon, Rectangle
from pydune.data_processing import plot_flux_rose, velocity_to_shear
from pydune.math import (
    cosd,
    make_angular_PDF,
    sind,
    vector_average,
)
from pydune.physics import quartic_transport_law
from pydune.physics.dune.courrechdupont2014 import (
    MGBNT_orientation,
    elongation_direction,
)


def Rotation_matrix(theta):
    return np.array([[cosd(theta), -sind(theta)], [sind(theta), cosd(theta)]])


def north_arrow(fig, ax, center, length, length_small, width, radius, theta=0, color="k"):
    y_start = radius + length - length_small
    arrow = np.array([[0, y_start], [width / 2, radius], [0, radius + length], [-width / 2, radius], [0, y_start]])
    # barycentre = np.sum(arrow, axis=0) / arrow.shape[0]
    # arrow = np.dot(Rotation_matrix(theta), (arrow-barycentre).T).T + barycentre
    arrow = np.dot(Rotation_matrix(theta + 180), arrow.T).T
    arrow = arrow + np.array(center)
    ax.add_patch(Polygon(arrow, color=color))
    t = ax.text(center[0], center[1], "N")
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    bb = t.get_window_extent(renderer=r).transformed(inv)
    t.set_visible(False)
    # width_t = bb.width
    height_t = bb.height
    t = ax.text(center[0], center[1] - height_t / 2, r"N", weight="bold", color=color, ha="center")


def plot_arrow(ax, point, length, angle, type, arrowprops):
    dx = cosd(angle) * length
    dy = sind(-angle) * length
    # dx = int(round(cosd(angle)*length))
    # dy = int(round(sind(-angle)*length))
    if type == "centered":
        xy = point - np.array([dx, dy]) / 2
        xytext = point + np.array([dx, dy]) / 2
    else:
        xy = point
        xytext = xy + np.array([dx, dy])
    arrow = FancyArrowPatch(xytext, xy, **arrowprops)
    ax.add_patch(arrow)
    if arrow.get_linestyle() != "-":
        # Tail
        v1 = arrow.get_path().vertices[0:3, :]
        c1 = arrow.get_path().codes[0:3]
        p1 = matplotlib.path.Path(v1, c1)
        pp1 = PathPatch(
            p1, color=arrow.get_facecolor(), lw=arrow.get_linewidth(), linestyle=arrow.get_linestyle(), fill=False
        )
        ax.add_patch(pp1)
        # Heads ====> partie qui ne marche pas
        v2 = arrow.get_path().vertices[3:, :]
        c2 = arrow.get_path().codes[3:]
        c2[0] = 1
        p2 = matplotlib.path.Path(v2, c2)
        pp2 = PathPatch(p2, color=arrow.get_facecolor(), lw=arrow.get_linewidth(), linestyle="-")
        ax.add_patch(pp2)
        arrow.remove()


path_data = "/media/cyril/Backup_dune/DUNE/PhD_Parts/Thèse/Chapitres/chapter5/figures/Figure_comparison_shape/Data"
Sites = ["Transverse", "Linear"]

####### Winds
orientation = np.load(os.path.join(path_data, "Uorientation_extracted_points_2.npy"), allow_pickle=True)
velocity = np.load(os.path.join(path_data, "Ustrength_extracted_points_2.npy"), allow_pickle=True)
points_wind = list(np.load(os.path.join(path_data, "points_2.npy"), allow_pickle=True))

###### parameters
left_point = [950, 700, 400]
# list_numb = ['A', 'B', 'C']
list_numb = ["(a)", "(b)", "(c)"]
# list_numb = ['(a)', '(b)', '(c)']
scale = [1.1, 1.7, 1.4]
xy_arrows_BI = [np.array([615, 460]), np.array([600, 460])]
xy_arrows_F = [np.array([np.nan, np.nan]), np.array([np.nan, np.nan])]
# color = ['w', 'w', 'w']
color = ["k", "k", "k"]
North = [0, 0, 0]
bins = [20, 20, 20]
position = [0, 1, 2]
Wind_points = [898, 852, 858]
#
unit = ["km", "km", "km"]
Data = {}
path_gen = os.getcwd()

##### Constants
grain_diameter = 180e-6  # grain diameter [m]
rho_s = 2.65e3  # grain density [kg/m3]
rho_f = 1.2  # fluid density   [kg/m3]
z_era = 10  # heght of wind data [m]
z0_era = 1e-3  # Aerodynamic roughness
secTOyear = 365.25 * 24 * 3600  # number of second in one year
gamma = 1.7  # constant for courrech model
g = 9.81
bed_porosity = 0.6

# %%
################################################### Figure
figsize = (0.5 * quarto.regular_fig_width, 0.28 * quarto.regular_fig_height)

fig, axarr = plt.subplots(1, 2, layout="constrained", figsize=figsize)
for i, ax in enumerate(axarr.flatten()):
    print(Sites[i])

    ax.set_xticks([])
    ax.set_yticks([])

    path = os.path.join(path_data, Sites[i])
    os.chdir(path)
    img = np.array(plt.imread(glob.glob("*.png")[0]))[300:1000, left_point[i] : left_point[i] + 1120, :]
    os.chdir(path_gen)
    ax.imshow(img)

    ######## échelle
    backgrnd = Rectangle((0, 0), width=0.25, height=0.2, transform=ax.transAxes, color="w", alpha=0.6)
    ax.add_patch(backgrnd)
    ax.plot([30, 30 + 230], [665, 665], linewidth=linewidth_barscale, color=color[i])
    ax.text(30 + 230 / 2, 645, f"{scale[i]:.0f} {unit[i]}", color=color[i], ha="center", fontsize=8)

    ########## Numbering
    length = 90
    length_small = 0.8 * length
    width = 55
    radius = 35
    center = np.array([1065, 625])
    #
    backgrnd = Rectangle((0.9, 0), width=0.1, height=0.325, transform=ax.transAxes, color="w", alpha=0.6)
    ax.add_patch(backgrnd)
    north_arrow(fig, ax, center, length, length_small, width, radius, theta=North[i], color=color[i])
    #
    ############################ flux rose
    ori = orientation[points_wind.index(Wind_points[i])]
    strength = velocity[points_wind.index(Wind_points[i])]
    #
    shear_velocity = velocity_to_shear(strength, z_era)
    #
    Q = np.sqrt((rho_s - rho_f) * g * grain_diameter / rho_f) * grain_diameter  # characteristic flux [m2/s]
    shield_th_quartic = 0.0035  # threshold shield numbers for the quartic

    # shield number
    shield = (rho_f / ((rho_s - rho_f) * g * grain_diameter)) * shear_velocity**2
    # dimensional sand flux, [m2/day]
    sand_flux = (1 / bed_porosity) * Q * quartic_transport_law(shield, shield_th_quartic) * secTOyear
    # angular distribution
    angular_PDF, angles = make_angular_PDF(ori, sand_flux)
    #
    DP = np.mean(sand_flux)  # Drift potential, [m2/day]
    # Resultant drift direction [deg.] / Resultant drift potential, [m2/day]
    RDD, RDP = vector_average(ori, sand_flux)
    #
    Alpha_E = elongation_direction(angles, angular_PDF)
    Alpha_BI = MGBNT_orientation(angles, angular_PDF)
    #
    if i == 0:
        anchor = [0.5, 0.45, 0.55, 0.55]
    elif i == 1:
        anchor = [0.725, 0.45, 0.55, 0.55]
    subax = ax.inset_axes(bounds=anchor, transform=ax.transAxes)
    a = plot_flux_rose(
        angles,
        angular_PDF,
        ax=subax,
        fig=fig,
        edgecolor="k",
        linewidth=0.5,
        cmap=flux_cmap,
    )
    # #
    # a.annotate(
    #     "",
    #     (RDD * np.pi / 180, 0),
    #     (RDD * np.pi / 180, 0.85 * a.get_rmax()),
    #     arrowprops=dict(arrowstyle="<|-", shrinkA=0, shrinkB=0, color="saddlebrown"),
    # )
    a.grid(linewidth=0.4, color="k", linestyle="--")
    a.set_axisbelow(True)
    a.patch.set_alpha(0.4)
    a.set_xticklabels([])
    #
    ######### Dune orientations
    length = 400
    #
    # plot_arrow(ax, xy_arrows_BI[i], length, A_max, 'centered',  arrowprops = dict(arrowstyle="<|-|>", color = color_BI, shrinkA = 0, shrinkB = 0, lw = lw_arrow, mutation_scale = 10, alpha = 0.5))
    if i == 0:
        plot_arrow(
            ax,
            xy_arrows_BI[i],
            length,
            Alpha_BI,
            "centered",
            arrowprops=dict(arrowstyle="<|-|>", color=color_BI, shrinkA=0, shrinkB=0, lw=lw_arrow, mutation_scale=10),
        )
    if ~np.isnan(xy_arrows_F[i]).all():
        center = xy_arrows_F[i]
    else:
        center = xy_arrows_BI[i]
    if i == 1:
        plot_arrow(
            ax,
            center,
            length / 2,
            Alpha_E + 10,
            "not_centered",
            arrowprops=dict(arrowstyle="<|-", color=color_F, shrinkA=0, shrinkB=0, lw=lw_arrow, mutation_scale=10),
        )

figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), i)
fig.savefig(figname, dpi=1200)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
