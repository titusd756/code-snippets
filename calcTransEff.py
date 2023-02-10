####
#### Script to calculate and visualise the transport efficiency of protons
#### through the focusing elements.
####

####
#### Import packages
####
import h5py
import numpy as np  #
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from matplotlib import rcParams
from matplotlib.ticker import LogFormatter
from copy import copy
from matplotlib.colors import LinearSegmentedColormap
import palettable
from matplotlib import gridspec
import pybdsim

# rcParams.update({"figure.autolayout": True})
plt.rcParams["axes.axisbelow"] = False

####
#### Set the context to 'talk' to make figure suitable for poster/presentation
####

sns.set_context("talk")
sns.set_style(
    {
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.bottom": True,
        "xtick.top": True,
        "ytick.left": True,
        "ytick.right": True,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "axes.spines.right": True,
        "axes.spines.top": True,
        "axes.grid": False,
        "grid.color": ".8",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)

####
#### Define physical constants
####
mpc2 = const.physical_constants["proton mass energy equivalent in MeV"][0]
gamma_15MeV_m2p = 1 + 15 * 0.98 / mpc2
gamma_15MeV_p2p = 1 + 15 * 1.02 / mpc2

####
#### Create a function to wrap the data extraction from file and the plotting
#### ---> returns histogram of all particles, transmission efficiency, sum of trans. eff.
#### ---> particle lost at the end, histogram of particles not lost
####
def analyse_one_file(fname):
    df = pd.read_csv(fname, header=0, dtype=np.float64, delimiter="\s+")
    print(df.head())

    # df = df.loc[
    #     ((df["OK_G_start"] >= gamma_15MeV_m2p) & (df["OK_G_start"] <= gamma_15MeV_p2p))
    #     | (
    #         (df["NOK_G_start"] >= gamma_15MeV_m2p)
    #         & (df["NOK_G_start"] <= gamma_15MeV_p2p)
    #     )
    # ]

    ok_vx = df["OK_Bx_start"].to_numpy()
    ok_vy = df["OK_By_start"].to_numpy()
    ok_vz = df["OK_Bz_start"].to_numpy()
    ok_G_start = df["OK_G_start"].to_numpy()
    nok_vx = df["NOK_Bx_start"].to_numpy()
    nok_vy = df["NOK_By_start"].to_numpy()
    nok_vz = df["NOK_Bz_start"].to_numpy()
    nok_G_start = df["NOK_G_start"].to_numpy()

    df_after_screen = df.loc[df["NOK_z_end"] > 5.753]
    nok_after_screen_G_start = df_after_screen["NOK_G_start"].to_numpy()
    nok_after_screen_vx = df_after_screen["NOK_Bx_start"].to_numpy()
    nok_after_screen_vy = df_after_screen["NOK_By_start"].to_numpy()
    nok_after_screen_vz = df_after_screen["NOK_Bz_start"].to_numpy()

    nok_z_end = df["NOK_z_end"].to_numpy()

    ok_v_abs = np.sqrt(np.power(ok_vx, 2) + np.power(ok_vy, 2) + np.power(ok_vz, 2))
    nok_v_abs = np.sqrt(np.power(nok_vx, 2) + np.power(nok_vy, 2) + np.power(nok_vz, 2))
    nok_after_screen_v_abs = np.sqrt(
        np.power(nok_after_screen_vx, 2)
        + np.power(nok_after_screen_vy, 2)
        + np.power(nok_after_screen_vz, 2)
    )

    ok_angle_start = (
        np.arccos(np.divide(ok_vz, ok_v_abs)) * 180 / np.pi
    )  # * np.sign(ok_vy)
    ok_ene_start = (ok_G_start - 1) * const.physical_constants[
        "proton mass energy equivalent in MeV"
    ][0]

    nok_after_screen_angle_start = (
        np.arccos(np.divide(nok_after_screen_vz, nok_after_screen_v_abs)) * 180 / np.pi
    )  # * np.sign(ok_vy)
    nok_after_screen_ene_start = (
        nok_after_screen_G_start - 1
    ) * const.physical_constants["proton mass energy equivalent in MeV"][0]

    nok_angle_start = (
        np.arccos(np.divide(nok_vz, nok_v_abs)) * 180 / np.pi  # * np.sign(nok_vy)
    )
    nok_ene_start = (nok_G_start - 1) * const.physical_constants[
        "proton mass energy equivalent in MeV"
    ][0]

    all_ene = np.concatenate((ok_ene_start, nok_ene_start))
    all_angle = np.concatenate((ok_angle_start, nok_angle_start))

    ok_ene_start = np.concatenate((ok_ene_start, nok_after_screen_ene_start))
    ok_angle_start = np.concatenate((ok_angle_start, nok_after_screen_angle_start))

    (h_all, xedges_all, yedgex_all, im_all) = plt.hist2d(
        all_angle, all_ene, range=[[0, 15], [0, 30]], bins=[60, 60]
    )

    plt.figure()
    (h_ok, xedges_ok, yedgex_ok, im_ok) = plt.hist2d(
        ok_angle_start, ok_ene_start, range=[[0, 15], [0, 30]], bins=[60, 60]
    )
    plt.colorbar()
    trans_eff = np.divide(h_ok, h_all, out=np.zeros_like(h_ok), where=h_all != 0) * 100

    plt.figure()
    plt.imshow(h_ok.T)

    h_ok_sum = np.sum(h_ok, axis=0)
    h_all_sum = np.sum(h_all, axis=0)
    trans_eff_sum = (
        np.divide(
            h_ok_sum, h_all_sum, out=np.zeros_like(h_ok_sum), where=h_all_sum != 0
        )
        * 100
    )
    return h_all, trans_eff, trans_eff_sum, nok_z_end, h_ok

####
#### Analyse the files from simulation based on their seed id
####
n_seeds = 1
for seed_idx in np.arange(1, 11, 1):
    if seed_idx == 1:
        h_all, trans_eff, trans_eff_sum, nok_z_end, h_ok = analyse_one_file(
            f"first3Solenoids_out_seedx{seed_idx}.txt"  # f"first3Solenoids_out_seedx{seed_idx}.txt"  # f"first3Lenses_out_500k_ptcls_seedx{seed_idx}.txt"  # f"first3Lenses_out_500k_ptcls_seedx{seed_idx}_screenAtGL1Exit.txt"  # f"test_out_seedx{seed_idx}.txt"
        )
    else:
        h_all_i, trans_eff_i, trans_eff_sum_i, nok_z_end_i, h_ok_i = analyse_one_file(
            f"first3Solenoids_out_seedx{seed_idx}.txt"  # f"first3Lenses_out_500k_ptcls_seedx{seed_idx}.txt"
        )
        h_all = h_all + h_all_i
        trans_eff = trans_eff + trans_eff_i
        trans_eff_sum = trans_eff_sum + trans_eff_sum_i
        nok_z_end = np.concatenate((nok_z_end, nok_z_end_i))
        h_ok = h_ok + h_ok_i
        n_seeds += 1
trans_eff /= n_seeds
trans_eff_sum /= n_seeds

# trans_eff = np.divide(h_ok, h_all, out=np.zeros_like(h_ok), where=h_all != 0) * 100

# h_all, trans_eff, trans_eff_sum, nok_z_end = analyse_one_file(
#     "first3Lenses_out_500k_ptlcs_seedx1_acc5_SCgrid50.txt"
# )

####
#### Visualise the results and save relevant figures
####

fig = plt.figure(figsize=(10, 7), dpi=120)

gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
gs.update(hspace=0.01)  # set the spacing between axes.
ax1 = plt.subplot(gs[0])
plt.gca().set_prop_cycle("color", palettable.cartocolors.qualitative.Bold_4.mpl_colors)
plt.plot(np.arange(0.25, 30, 0.5), trans_eff_sum)
plt.gca().set_yticks([25, 50, 75, 100])
plt.title(r"Transmission efficiency $\mathrm{\eta}$")
plt.ylim([0, 100])
plt.ylabel(r"$\mathrm{d\eta / dE_p} \,[\%]$")

cmap = mpl.cm.get_cmap("nipy_spectral_r")
# cmap.set_under(color="white")

colors = plt.cm.get_cmap("jet")(np.linspace(0, 0.9, 128))
colors = list(zip(np.linspace(0.0, 1.0, 128), colors))
# colors_ = [(0, "white")] + colors
cmap = mpl.colors.LinearSegmentedColormap.from_list("mycmap", colors)

wh = np.where(h_all > 0, np.ones_like(h_all), np.zeros_like(h_all))
# plt.imshow(
#     wh,
#     # cmap=,
#     vmin=0,
#     origin="lower",
#     extent=[0, 30, 0, 15],
#     aspect=1.2,
#     # alpha=wh,
#     # zorder=10,
# )
plt.subplot(gs[1], sharex=ax1)
plt.imshow(
    trans_eff,
    cmap=cmap,  # cmap="cividis",  # cmap=palettable.scientific.sequential.Bamako_10_r.mpl_colormap,
    vmin=0,
    origin="lower",
    extent=[0, 30, 0, 15],
    aspect=1.2,
    alpha=wh,
    zorder=2,
)
plt.imshow(
    trans_eff,
    cmap=cmap,  # cmap="cividis",  # cmap=palettable.scientific.sequential.Bamako_10_r.mpl_colormap,
    vmin=0,
    origin="lower",
    extent=[0, 30, 0, 15],
    aspect=1.2,
    alpha=0.0,
    zorder=1,
)

plt.subplots_adjust(left=0.1, right=0.8, hspace=0.001)
plt.ylabel(r"Divergence angle $\mathrm{\theta_z}$ $\left[\mathrm{^{\circ}}\right]$")
plt.xlabel(r"$\mathrm{E_p}$ [MeV]")
plt.ylim([0, 10])
# plt.title(r"Transmission efficiency")
ax = plt.gca()
cax = fig.add_axes(
    [ax.get_position().x1 + 0.02, ax.get_position().y0, 0.04, ax.get_position().height]
)
# plt.colorbar(cax=cax)  # Similar to fig.colorbar(im, cax = cax)
cbar = plt.colorbar(cax=cax, label=r"$\mathrm{d\eta / d\theta_z / dE_p} [\%]$")
cbar.set_alpha(1)
cbar.draw_all()

# plt.savefig("transf_eff_sol_avg_10samples_aligned.png", dpi=300)

#########
fig = plt.figure(figsize=(10, 7), dpi=120)

gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
gs.update(hspace=0.01)  # set the spacing between axes.
ax1 = plt.subplot(gs[0])
plt.gca().set_prop_cycle("color", palettable.cartocolors.qualitative.Bold_4.mpl_colors)
plt.plot(np.arange(0.25, 30, 0.5), trans_eff_sum)
plt.gca().set_yticks([25, 50, 75, 100])
plt.title(r"Transmission efficiency $\mathrm{\eta}$")
plt.ylim([0, 100])
plt.ylabel(r"$\mathrm{d\eta / dE_p} \,[\%]$")

cmap = mpl.cm.get_cmap("nipy_spectral_r")
# cmap.set_under(color="white")

colors = plt.cm.get_cmap("jet")(np.linspace(0, 0.9, 128))
colors = list(zip(np.linspace(0.0, 1.0, 128), colors))
# colors_ = [(0, "white")] + colors
cmap = mpl.colors.LinearSegmentedColormap.from_list("mycmap", colors)

wh = np.where(h_all > 0, np.ones_like(h_all), np.zeros_like(h_all))
# plt.imshow(
#     wh,
#     # cmap=,
#     vmin=0,
#     origin="lower",
#     extent=[0, 30, 0, 15],
#     aspect=1.2,
#     # alpha=wh,
#     # zorder=10,
# )
plt.subplot(gs[1], sharex=ax1)
plt.imshow(
    h_all,
    cmap=cmap,  # cmap="cividis",  # cmap=palettable.scientific.sequential.Bamako_10_r.mpl_colormap,
    vmin=0,
    origin="lower",
    extent=[0, 30, 0, 15],
    aspect=1.2,
    alpha=wh,
    zorder=2,
)
plt.imshow(
    h_all,
    cmap=cmap,  # cmap="cividis",  # cmap=palettable.scientific.sequential.Bamako_10_r.mpl_colormap,
    vmin=0,
    origin="lower",
    extent=[0, 30, 0, 15],
    aspect=1.2,
    alpha=0.0,
    zorder=1,
)

plt.subplots_adjust(left=0.1, right=0.8, hspace=0.001)
plt.ylabel(r"divergence $\theta$ $\mathrm{\Omega}$ [$\mathrm{^\circ}$]")
plt.xlabel(r"$\mathrm{E_p}$ [MeV]")
plt.ylim([0, 10])
# plt.title(r"Transmission efficiency")
ax = plt.gca()
cax = fig.add_axes(
    [ax.get_position().x1 + 0.02, ax.get_position().y0, 0.04, ax.get_position().height]
)
# plt.colorbar(cax=cax)  # Similar to fig.colorbar(im, cax = cax)
cbar = plt.colorbar(cax=cax, label=r"$\mathrm{d\eta / d\theta / dE_p} [\%]$")
cbar.set_alpha(1)
cbar.draw_all()
#########


fig = plt.figure(figsize=(7, 4), dpi=120)

np.save("nok_z_end_sol.npy", nok_z_end)

nok_z_end_sol = np.load("nok_z_end_sol.npy")
nok_z_end_lens = np.load("nok_z_end_lens.npy")

N1 = len(nok_z_end_sol)
N2 = len(nok_z_end_sol)

plt.gca().set_prop_cycle("color", palettable.cartocolors.qualitative.Bold_4.mpl_colors)
n1, bins1, pat1 = plt.hist(
    nok_z_end_sol,
    bins=100,
    weights=1 / 1e7 * np.ones_like(nok_z_end_sol) * 100,
    facecolor="None",
    lw=2,
    edgecolor=palettable.cartocolors.qualitative.Bold_4.mpl_colors[0],
    histtype="step",
    label="solenoids",
    cumulative=True,
)
# plt.errorbar(
#     bins1[:-1] + (bins1[1] - bins1[0]) / 2 * np.ones_like(bins1[:-1]),
#     n1,
#     yerr=np.sqrt(n1 - 1 / 1 * np.power(n1 * 1, 2) / 1e7) / 1e7,
#     marker="s",
#     color=palettable.cartocolors.qualitative.Bold_4.mpl_colors[0],
#     capsize=5,
#     linestyle="-",
# )
n2, bins2, pat2 = plt.hist(
    nok_z_end_lens,
    bins=100,
    weights=1 / (1e7) * np.ones_like(nok_z_end_lens) * 100,
    facecolor="None",
    lw=2,
    edgecolor=palettable.cartocolors.qualitative.Bold_4.mpl_colors[1],
    histtype="step",
    label="plasma lenses",
    cumulative=True,
)
# plt.errorbar(
#     bins2[:-1] + (bins2[1] - bins2[0]) / 2 * np.ones_like(bins2[:-1]),
#     n2,
#     yerr=np.sqrt(n2 * 1e7) / 1e7,
#     marker="+",
#     color=palettable.cartocolors.qualitative.Bold_4.mpl_colors[1],
#     capsize=0,
#     linestyle="",
# )
plt.xlabel(r"z / m")
plt.ylabel(r"Cumulative fractional" + "\n" + r"beam loss ($\%$)")
# plt.title(r"Loss map")
# plt.semilogy()
# plt.ylim([1e-6, -1])
plt.legend(loc="lower right")
# plt.gca().set_yticks([1e-5, 1e-4, 1e-3, 1e-2])
pybdsim.Plot.AddMachineLatticeFromSurveyToFigure(
    fig, "survey.dat", tightLayout=True, sOffset=0.0
)

plt.xlim([0, 5.753])
plt.savefig("loss_map_combined_cumulative_15MeVpm2p.png", dpi=300)


plt.show()
