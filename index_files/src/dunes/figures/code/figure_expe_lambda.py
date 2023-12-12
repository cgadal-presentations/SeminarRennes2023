import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import scipy.io
import svgutils.transform as sg
from pydune.math import tand
from pydune.physics.dune.bedinstability_1D import temporal_growth_rate
from scipy.optimize import curve_fit

Theta_exp = [0, 30, 50, 70, 90, 110, 130, 150, 180]
r_val = [1.005, 1.01, 1.2, 1.6, 2.5, 5, 7, 10, 15]
Dat = {}
pix2mil = 0.45
pathI = "/media/cyril/Backup_dune/DUNE/PhD_Parts/Part1_multidirectional/IPGP_intership/Experimental/"

delta = 0.0
beta = 0

A0 = 4.5
B0 = 3

# ############################################# importing data for courbe lambda/theta
# ######### Errorbars
# for th in Theta_exp :
#     path = pathI + '/Lit_plat/SerieV240/Theta' + str(th) + '/results_stat.txt'
#     a = np.loadtxt(path)
#     Dat[('Theta',th)] = a[0]
#     Dat[('Lambda',th)] = a[1]*pix2mil
#     Dat[('Std',th)] = max(2*a[2],1)*pix2mil
#     # Dat[('Time',th)] = a[3]
#     path = pathI + '/Lit_plat/SerieV240/Theta' + str(th) + '/Worskpace_finding_time.mat'
#     mat = scipy.io.loadmat(path)
#     Dat[('Lambda_corr',th)] = mat['Lambda'][0,-1]*pix2mil
#     Dat[('Time_corr',th)] = mat['kk'][0]+ (float(mat['Nstart'])-1)/3
#
#
# ######## correction pour Theta = 70
# Dat[('Std',70)] = 0.8*Dat[('Std',70)]


# ######### Theory
# path_dat = '/home/gadal/Documents/Work/Research/PhD_Parts/Part1_multidirectional/Paper/Resubmission/ReplytoRef/Analyse_spp/Other_transport_law/delta_0_beta_0/'
#
# n = 1
# Theta = np.linspace(0, 180, 181)
# Temp = np.load(path_dat + '/Data/DataN'+str(n)+'.npy', allow_pickle = True)
# Temp2 = Temp.item()
# Kmax = {}
# Sigmax = {}
# for th in Theta:
#     for r in r_val:
#         Sigmax[(n,th,r)] = Temp2['SigmaMax'][(n,th,r)]
#         #Amax[(n,th,r)] = Temp2['Amax'][(n,th,r)]
#         Kmax[(n,th,r)] = Temp2['Kmax'][(n,th,r)]

############################################## importing data for courbe lambda/V
Dattrans = {}

Vtrans = [120, 140, 180, 200, 240, 260]
################## Errorbars
for v in Vtrans:
    if v == 240 or v == 120:  #### 120 to fill with smth
        #        path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SerieV240/Theta0/results_stat.txt'
        path = pathI + "/Lit_plat/SerieV240/Theta0/Worskpace_finding_time.mat"
    elif v == 260:
        # path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SérieTheta0/V' + str(v) + 'bis/results_stat.txt'
        path = pathI + "/Lit_plat/SérieTheta0/V" + str(v) + "bis/Worskpace_finding_time.mat"
    else:
        #       path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SérieTheta0/V' + str(v) + '/results_stat.txt'
        path = pathI + "/Lit_plat/SérieTheta0/V" + str(v) + "/Worskpace_finding_time.mat"
    #    a = np.loadtxt(path)
    #    Dattrans[('Lambda',v)] = a[1]*pix2mil
    #    Dattrans[('Std',v)] = max(2*a[2],1)*pix2mil
    mat = scipy.io.loadmat(path)
    Dattrans[("Lambda_corr", v)] = mat["Lambda"][0, -1] * pix2mil
    Dattrans[("Time_corr", v)] = mat["kk"][0] + float(mat["Nstart"])
    if v == 240:  #### 120 to fill with smth
        path = pathI + "/Lit_plat/SerieV240/Theta0/results_stat.txt"
    #    path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SerieV240/Theta0/Worskpace_finding_time.mat'
    elif v == 260 or v == 120:
        path = pathI + "/Lit_plat/SérieTheta0/V" + str(v) + "bis/results_stat.txt"
    #    path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SérieTheta0/V' + str(v) + 'bis/Worskpace_finding_time.mat'
    else:
        path = pathI + "/Lit_plat/SérieTheta0/V" + str(v) + "/results_stat.txt"
    #    path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SérieTheta0/V' + str(v) + '/Worskpace_finding_time.mat'
    a = np.loadtxt(path)
    Dattrans[("Lambda", v)] = a[1] * pix2mil
    Dattrans[("Std", v)] = max(2 * a[2], 1) * pix2mil

mat = scipy.io.loadmat(pathI + "/Lit_plat/SérieTheta0/V120bis/Workspace_rot_img.mat")
Dattrans[("Lambda_corr", 120)] = mat["Lambda"][0, 32] * pix2mil
Dattrans[("Std", 120)] = 4 * pix2mil  #####""estimation via images

##################### Data final
##
# path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SérieTheta0/Lambda_Time_final.mat'
# mat = scipy.io.loadmat(path)
# Lambda_Trans_final = list(mat['Lambda_final'][0]*pix2mil)
# Lambda_Trans_final.insert(0,Dattrans[('Lambda',120)])
# Lambda_Trans_final.insert(4,Lambda_Theta[0][0])
#
# Time_Trans = mat['Time_final']

##################### Velocity scaled


Dat_scaling = np.loadtxt(pathI + "/Scaling_manip_velocity/scaling_velocity_manip.txt")
Vtrans_mes = np.empty(
    [
        len(Vtrans),
    ]
)
for j in range(len(Vtrans)):
    #    if np.argwhere(Vtrans[j] == Dat_scaling[0,:]).size == 0:
    #        #### lin interp
    #        Vtrans_mes[j] = Dat_scaling[1,np.argwhere(180 == Dat_scaling[0,:])] + (Vtrans[j] - 180)*(Dat_scaling[1,np.argwhere(240 == Dat_scaling[0,:])]-Dat_scaling[1,np.argwhere(180 == Dat_scaling[0,:])])/(240-180)
    #    else:
    #        Vtrans_mes[j] = Dat_scaling[1,np.argwhere(Vtrans[j] == Dat_scaling[0,:])]
    Vtrans_mes[j] = Dat_scaling[1, np.argwhere(Vtrans[j] == Dat_scaling[0, :])]

Datlong = {}
Vlong = [120, 140, 180, 240]

##################### Errorbars
for v in Vlong:
    if v == 240:
        #        path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SerieV240/Theta180/results_stat.txt'
        path = pathI + "/Lit_plat/SerieV240/Theta180/Worskpace_finding_time.mat"
    else:
        #        path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SérieTheta180/V' + str(v) + '/results_stat.txt'
        path = pathI + "/Lit_plat/SérieTheta180/V" + str(v) + "/Worskpace_finding_time.mat"
    #    a = np.loadtxt(path)
    #    Datlong[('Lambda',v)] = a[1]*pix2mil
    #    Datlong[('Std',v)] = max(2*a[2],1)*pix2mil
    mat = scipy.io.loadmat(path)
    Datlong[("Lambda_corr", v)] = mat["Lambda"][0, -1] * pix2mil
    Datlong[("Time_corr", v)] = mat["kk"][0] + float(mat["Nstart"])
    if v == 240:
        path = pathI + "/Lit_plat/SerieV240/Theta180/results_stat.txt"
    #        path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SerieV240/Theta180/Worskpace_finding_time.mat'
    else:
        path = pathI + "/Lit_plat/SérieTheta180/V" + str(v) + "/results_stat.txt"
    #        path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SérieTheta180/V' + str(v) + '/Worskpace_finding_time.mat'
    a = np.loadtxt(path)
    Datlong[("Lambda", v)] = a[1] * pix2mil
    Datlong[("Std", v)] = max(2 * a[2], 1) * pix2mil

##################### Data final
#
# path = '/home/cyril/Documents/Work/Research/IPGP_Intership/Experimental/Lit_plat/SérieTheta180/Lambda_Time_final.mat'
# mat = scipy.io.loadmat(path)
# Lambda_long_final = list(mat['Lambda_final'][0]*pix2mil)
# Lambda_long_final.insert(4,Lambda_Theta[0][-1])
#
Vlong_mes = np.empty(
    [
        len(Vlong),
    ]
)
for j in range(len(Vlong)):
    #    if np.argwhere(Vlong[j] == Dat_scaling[0,:]).size == 0:
    #        #### lin interp
    #        Vlong_mes[j] = Dat_scaling[1,np.argwhere(180 == Dat_scaling[0,:])] + (Vlong[j] - 180)*(Dat_scaling[1,np.argwhere(240 == Dat_scaling[0,:])]-Dat_scaling[1,np.argwhere(180 == Dat_scaling[0,:])])/(240-180)
    #    else:
    #        Vlong_mes[j] = Dat_scaling[1,np.argwhere(Vlong[j] == Dat_scaling[0,:])]
    Vlong_mes[j] = Dat_scaling[1, np.argwhere(Vlong[j] == Dat_scaling[0, :])]
#

###################

Vth = 14.1

# Lambda_trans = list({k:v for k,v in Dattrans.items() if k[0] == 'Lambda'}.values())
# Lambda_trans = Lambda_Trans_final
Lambda_trans = list({k: v for k, v in Dattrans.items() if k[0] == "Lambda_corr"}.values())
Std_trans = list({k: v for k, v in Dattrans.items() if k[0] == "Std"}.values())


def lambda_fit_test(v, G1, G2):
    mu = tand(35)
    r = G2 * (v - 1) + 1
    C4 = 1.0 / r
    C3 = 1.0 - C4**2
    C1 = 1.0 + delta * C3
    C2 = C4**2 * (1 - beta) + beta
    return G1 * C1 * A0 / (C1 * B0 - (C2 / mu))


p1, pconv1 = curve_fit(
    lambda_fit_test,
    Vtrans_mes / Vth,
    Lambda_trans,
    p0=[1, 10],
    check_finite=True,
    bounds=([0, 0], [20, 20]),
    sigma=Std_trans,
)


def lambda_fit(V, C1):
    C2 = p1[1]
    mu = tand(35)
    # R = C2 * (V -1) + 1
    k = np.linspace(0.0, 0.6, 1001)
    kM = []
    for v in V:
        Sigma = temporal_growth_rate(k, A0, B0, mu, C2 * (v - 1) + 1)
        kM.append(k[np.argwhere(Sigma == np.max(Sigma))[0][0]])
    return 2 * np.pi * C1 / np.array(kM)


p, pconv = curve_fit(lambda_fit, Vtrans_mes / Vth, Lambda_trans, p0=[1], check_finite=True, sigma=Std_trans)
np.savetxt("Lsat_value.txt", p)


def vtou(v, p):
    return p * (v - 1) + 1


Rplot = min(r_val, key=lambda x: abs(x - vtou(1.4, p1[1])))
# %%

#########################################

figsize = (0.6 * quarto.regular_fig_width, 0.5 * quarto.regular_fig_height)

fig, ax = plt.subplots(1, 1, layout="constrained", figsize=figsize)

a_mark = ax.errorbar(
    vtou(Vtrans_mes / Vth, p1[1]),
    Lambda_trans,
    yerr=Std_trans,
    fmt=".",
    # linewidth=linewidth_errorbar,
    label=r"$\lambda_{0}$",
)
(a_line,) = ax.plot(vtou(Vtrans_mes / Vth, p1[1]), Lambda_trans, color="#1f77b4", alpha=0.4, linestyle="--")

# Lambda_long = list({k:v for k,v in Datlong.items() if k[0] == 'Lambda'}.values())
Lambda_long = list({k: v for k, v in Datlong.items() if k[0] == "Lambda_corr"}.values())
# Lambda_long = Lambda_long_final
Std_long = list({k: v for k, v in Datlong.items() if k[0] == "Std"}.values())

# b_mark = ax.errorbar(
#     Vlong_mes / Vth,
#     Lambda_long,
#     yerr=Std_long,
#     fmt="s",
#     markerfacecolor="None",
#     # linewidth=linewidth_errorbar,
#     label=r"$\lambda_{180}$",
# )
# (b_line,) = ax.plot(Vlong_mes / Vth, Lambda_long, color="#ff7f0e", alpha=0.4, linestyle="--")

rplot = np.linspace(1, 1.5, 100)
(c_line,) = ax.plot(vtou(rplot, p1[1]), lambda_fit(rplot, *p), "k", label="fit")

ax.set_xlim([1, 3.5])
ax.set_yticks([8, 16, 24])
ax.set_ylim([7, 25])
ax.set_xlabel(r"$u_{*}/u_{\rm th}$")
ax.set_ylabel(r"$\lambda_{\rm max}$ [mm]")

# %% saving figure
figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=1200)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
