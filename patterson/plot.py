"""
Plotting utilities for XRR and Patterson function data.
"""

import glob

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import cycler

from .core import _fresnel_reflectivity, _wrap_label
from .utils import natural_keys


def plotXrr(dataFile, rho_sub=0.71, rho_pre=0.0, rffPlot=True, scale=1, ax=None):
    """
    Plot XRR data, optionally normalized by Fresnel reflectivity.

    Parameters
    ----------
    dataFile : str
        Path to a two-column (q, intensity) data file.
    rho_sub, rho_pre : float
        Electron densities for Fresnel reflectivity calculation.
    rffPlot : bool
        If True plot R/R_F; otherwise plot raw R.
    scale : float
        Multiplicative scaling factor.
    ax : matplotlib.axes.Axes, optional
        Target axes. Uses current axes if None.
    """
    data = np.loadtxt(dataFile)
    qq = data[:, 0]
    ii = data[:, 1]

    target_ax = ax if ax is not None else plt.gca()

    label = _wrap_label(dataFile)
    if rffPlot:
        rf = _fresnel_reflectivity(qq, rho_sub, rho_pre)
        target_ax.semilogy(qq, ii / rf * scale, label=label)
        target_ax.set_ylabel(r"R/R$_\mathrm{F}$")
    else:
        target_ax.semilogy(qq, ii * scale, label=label)
        target_ax.set_ylabel("R")

    target_ax.set_xlabel(r"q$_z$ (Å$^{-1}$)")
    target_ax.legend()


def plot_all(
    rho_sub=0.71,
    rho_pre=0.0,
    xmin1=0,
    xmax1=0.7,
    xmin2=0,
    xmax2=350,
    ymin2=0,
    ymax2=4e-4,
    save=True,
    plotLegend=True,
    doScale=True,
    norm=True,
    outputFile="foo",
    outputDir=".",
    vline=None,
):
    """
    Create a side-by-side plot of all XRR and Patterson function files in the
    current working directory.

    Parameters
    ----------
    rho_sub, rho_pre : float
        Electron densities for Fresnel reflectivity calculation.
    xmin1, xmax1 : float
        q-axis limits for the XRR panel.
    xmin2, xmax2 : float
        Distance-axis limits for the Patterson panel.
    ymin2, ymax2 : float
        y-axis limits for the Patterson panel.
    save : bool
        Save the figure to *outputDir/outputFile.png*.
    plotLegend : bool
        Show legend on the Patterson panel.
    doScale : bool
        Stack Patterson curves by an auto-calculated offset.
    norm : bool
        Normalize each Patterson curve to its maximum.
    outputFile : str
        Base name for the saved figure.
    outputDir : str
        Directory for the saved figure.
    vline : float or None
        Draw a vertical reference line on the Patterson panel at this x value.
    """
    fileList_xr = sorted(glob.glob("*.xrr"), key=natural_keys)
    fileList_ft = sorted(glob.glob("*.ftMan"), key=natural_keys)

    color = mpl.cm.gnuplot(np.linspace(0, 1, max(len(fileList_xr), 1)))
    mpl.rcParams["axes.prop_cycle"] = cycler.cycler("color", color)

    plt.rc("font", size=14)
    plt.rc("legend", fontsize=14)
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"

    fig = plt.figure(figsize=(9.3, 6))
    gs = gridspec.GridSpec(1, 2)
    gs.update(wspace=0.03)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    ax1.set_xlabel(r"q$_\mathregular{z}$ (Å$\mathregular{^{-1}}$)")
    ax1.set_ylabel(r"R/R$_\mathregular{F}$")
    ax1.yaxis.set_ticks_position("both")
    ax1.xaxis.set_ticks_position("both")

    ax2.set_xlabel("d (Å)")
    ax2.set_ylabel(r"FT(R/R$_\mathregular{F}$)", rotation=-90, labelpad=24)
    ax2.yaxis.set_ticks_position("both")
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()

    ax1.set_xlim(xmin1, xmax1)
    ax2.set_xlim(xmin2, xmax2)
    ax2.set_ylim(ymin2, ymax2)

    # XRR panel
    nn = 1
    for fileName in fileList_xr:
        data = np.loadtxt(fileName)
        qq = data[:, 0]
        ii = data[:, 1]
        rf = _fresnel_reflectivity(qq, rho_sub, rho_pre)
        ax1.semilogy(qq, ii / rf * nn, label=_wrap_label(fileName), linewidth=1)
        nn *= 100

    # Patterson panel
    offset = 0
    scale_step = 0
    for i, fileName in enumerate(fileList_ft):
        data = np.loadtxt(fileName)
        rr = data[:, 0]
        pp = data[:, 1]

        y = pp / np.max(pp) if norm else pp

        ax2.plot(rr, y + (offset if doScale else 0), label=_wrap_label(fileName), linewidth=1)

        if i == 0:
            scale_step = np.max(y) / 10
        if doScale:
            offset += scale_step

    if plotLegend:
        ax2.legend(fontsize=8)

    if fileList_xr:
        ax1.set_title(_wrap_label(fileList_xr[0]), fontsize=12)
        ax2.set_title(_wrap_label(fileList_xr[-1]), fontsize=12)

    if vline is not None:
        ax2.axvline(x=vline, color="k", alpha=0.3, linewidth=0.7)

    if save:
        import os
        outpath = os.path.join(outputDir, f"{outputFile}.png")
        fig.savefig(outpath, bbox_inches="tight", dpi=600)
        print(f"Saved: {outpath}")

    return fig, (ax1, ax2)
