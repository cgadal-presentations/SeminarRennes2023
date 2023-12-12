import glob
import os
import sys
from datetime import datetime

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from color_template_dune import color_max
from pydune.data_processing import velocity_to_shear
from pydune.math import (
    tand,
    vector_average,
)
from pydune.physics import quartic_transport_law
from pydune.physics.dune import bedinstability_1D
from scipy.optimize import curve_fit

matplotlib.rcParams["font.size"] = 8  # default 12
path = "/media/cyril/Backup_dune/DUNE/Clem_stuff/Clem_China/China_Clement"
# parameters
sec2yr = 365.25 * 24 * 3600
rho_s = 2.55e3  # grain density [kg/m3]
rho_f = 1.29  # fluid density   [kg/m3]
K = 0.4  # Von karman constant
z = 10  # heght of wind data [m]
z0 = 1e-3  # Aerodynamic roughness
g = 9.8  # gravity [m/s2]
grain_size = 190e-6  # grain size [m]
Lsat = 0.95  # Saturation length [m]
Lsatbis = 0.7  # Saturation length 2 [m]
Ames = 3  # Hydrodynamic coeff A0 measured
Bmes = 1.5  # Hydrodynamic coeff B0 measured
mu = tand(35)  # Avalanche slope
bed_porosity = 0.6

# %% wind_data
wind_data = np.loadtxt(os.path.join(path, "Data/all_201301_201710.txt"))
speed = wind_data[:, 4]  # wind velocity at 10 m
direction = wind_data[:, 5]  # wind direction at 10 m
#
sep = " "
dates_string = np.array(
    [
        f"{yr:04d}" + sep + f"{mth:02d}" + sep + f"{day:02d}" + sep + f"{hr:02d}"
        for [yr, mth, day, hr] in wind_data[:, :4].astype("int")
    ]
)
t_station = np.array([datetime.strptime(d, "%Y %m %d %H") for d in dates_string])

# filtering bad values
speed[direction > 360] = np.nan
direction[direction > 360] = np.nan
#


# %% computing regular quantities
shear_velocity = velocity_to_shear(speed, z)
#
Q = np.sqrt((rho_s - rho_f) * g * grain_size / rho_f) * grain_size  # characteristic flux [m2/s]
shield_th_quartic = 0.0035  # threshold shield numbers for the quartic

# shield number
shield = (rho_f / ((rho_s - rho_f) * g * grain_size)) * shear_velocity**2
# dimensional sand flux, [m2/day]
saturated_flux = (1 / bed_porosity) * Q * quartic_transport_law(shield, shield_th_quartic) * sec2yr

DP = np.nanmean(saturated_flux)  # Drift potential, [m2/day]
# Resultant drift direction [deg.] / Resultant drift potential, [m2/day]
RDD, RDP = vector_average(direction, saturated_flux)

# %% computing quantity for LSA

# threshold shear velocity [m/s]
shear_velocity_th = np.sqrt(shield_th_quartic / (rho_f / ((rho_s - rho_f) * g * grain_size)))
# average velocity ratio by angle bin
r_temporal = np.where(shear_velocity > shear_velocity_th, shear_velocity / shear_velocity_th, 1)
caracteristic_flux = saturated_flux / (1 - 1 / r_temporal**2)
DP_car = np.nanmean(caracteristic_flux)  # Drift potential, [m2/day]
RDD_car, RDP_car = vector_average(
    direction, caracteristic_flux
)  # Resultant drift direction [deg.] / Resultant drift potential, [m2/day]

# %%
# Quantities corresponding to Qcar = characteristic flux of instability - velocity period
tmin = datetime.strptime("2014 04 25", "%Y %m %d")
tmax = datetime.strptime("2014 05 13", "%Y %m %d")
mask = (t_station >= tmin) & (t_station <= tmax)
#
DP_car_short = np.nanmean(caracteristic_flux[mask])  # Drift potential, [m2/day]
RDD_car_short, RDP_car_short = vector_average(
    direction[mask], caracteristic_flux[mask]
)  # Resultant drift direction [deg.] / Resultant drift potential, [m2/day]

# %% Growth rate data
Datas = []
step_dirs = sorted(glob.glob(os.path.join(path, "Data/step*")))
for step in step_dirs:
    list_file = sorted(glob.glob(step + "/*"))
    for file in list_file:
        tp = np.loadtxt(file)
        Datas.append([float(step.split("_")[-1]), tp[0, 0], *tp[:, 1]])

Datas = np.array(Datas)  # Col1 = step, Col2 = lambda, Col3-end = growth rates adi (to be corrected by Q/lambda**2)
Datas[:, 2:] = Datas[:, 2:] * DP / 365.25 / 225 * np.log(10)  # Growth rates in [/day].


##########################################################
##################### Results ############################
##########################################################
Dim_fact_growth_rate = DP_car / Lsat**2 / 365.25
Dim_fact_growth_rate_bis = DP_car / Lsatbis**2 / 365.25
Dim_fact_celerity = DP_car_short / Lsat / 365.25
Dim_fact_celerity_bis = DP_car_short / Lsatbis / 365.25

#################################################### S = 0


def to_fit(k, A0, B0):
    r = 10000  # infinitely far from threshold, S = 0
    return bedinstability_1D.temporal_growth_rate(k, A0, B0, mu, r)


########### Fits
kc = 0.72  # cutoff wavenumber for fits [m]
inds = 2 * np.pi / Datas[:, 1] < 0.72  # mask for fitting only wavenumbers smaller than kc

# Non dimesnional fits:
p1, pcov = curve_fit(
    to_fit, 2 * np.pi / Datas[inds, 1] * Lsat, Datas[inds, 2:].mean(axis=1) / Dim_fact_growth_rate
)  # with Lsat = 0.95
p2, pcov = curve_fit(
    to_fit, 2 * np.pi / Datas[inds, 1] * Lsatbis, Datas[inds, 2:].mean(axis=1) / Dim_fact_growth_rate_bis
)  # with Lsatbis = 0.7


# #### figure
kplot = np.linspace(0, 1, 10000)  # wavenumbers for plot [sans dimensions]
sigma_th = to_fit(kplot, *p2) * Dim_fact_growth_rate_bis

xmax = kplot[sigma_th.argmax()] / Lsatbis
ymax = sigma_th.max()

xzero = kplot[kplot > 0.2][np.argmin(np.abs(sigma_th[kplot > 0.2]))] / Lsatbis

linewidth_small = 0.6

# %%
figsize = (0.7 * quarto.regular_fig_width, 0.6 * quarto.regular_fig_height)
fig, ax = plt.subplots(1, 1, layout="constrained", figsize=figsize)
ax.axhline(0, color="k", lw=linewidth_small)

# data
ax.plot(kplot / Lsatbis, to_fit(kplot, *p2) * Dim_fact_growth_rate_bis, label=r"theory")
lims = ax.get_ylim()

ax.errorbar(
    2 * np.pi / Datas[:, 1],
    Datas[:, 2:].mean(axis=1),
    yerr=Datas[:, 2:].std(axis=1) / np.sqrt(Datas[:, 2:].shape[1]),
    fmt=".",
    label="data",
)
# annotations
ax.plot(xmax, ymax, color=color_max, marker="+")
#
ax.vlines(xmax, ymin=lims[0], ymax=ymax, color=color_max, linewidth=linewidth_small, linestyle="--")
ax.text(1.03 * xmax, -5e-3, r"$\lambda_{\rm max}$", color=color_max, va="center")
#
ax.hlines(
    ymax,
    xmax=xmax,
    xmin=kplot.min() / Lsatbis,
    color=color_max,
    linewidth=linewidth_small,
    linestyle="--",
)
ax.text(0.02, 1.015 * ymax, r"$\sigma_{\rm max}$", color=color_max, ha="left", va="bottom")
#
ax.axvline(
    xzero,
    color="tab:red",
    # linewidth=linewidth_small,
    linestyle="--",
)
ax.text(1.02 * xzero, -10e-3, r"$\lambda_{\rm min}$", color="tab:red", va="center")
#
# ####
ax.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc="lower right", ncol=2, borderaxespad=0.0)
ax.set_xlabel(r"Wavenumber, $k~[\rm m^{-1}]$")
ax.set_ylabel(r"Growth rate, $\sigma~[\rm day ^{-1}]$")
ax.set_xlim([0, 0.9])
ax.set_ylim(bottom=-0.013, top=0.016)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

# %% saving figure
figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=1200)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
