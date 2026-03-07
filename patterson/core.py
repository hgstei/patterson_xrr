"""
Core Patterson function computation for XRR data analysis.
"""

import numpy as np


def _fresnel_reflectivity(qq, rho_sub=0.71, rho_pre=0.0):
    """Compute Fresnel reflectivity for a substrate/air interface."""
    delta = rho_sub - rho_pre
    qc_sq = (0.0375 * np.sqrt(delta)) ** 2
    rf = np.abs(
        (qq - np.sqrt(qq**2 - qc_sq)) / (qq + np.sqrt(qq**2 - qc_sq))
    ) ** 2
    return rf


def _normalize_intensity(qq, ii, rrf_norm, rho_sub, rho_pre):
    """Return R/R_F (rrfNorm=True) or R*q^4 (rrfNorm=False), NaN-safe."""
    if rrf_norm:
        rf = _fresnel_reflectivity(qq, rho_sub, rho_pre)
        rrf = ii / rf
    else:
        rrf = ii * qq**4
    np.place(rrf, np.isnan(rrf), 1)
    return rrf


def pSpectrum_nonPadded(data):
    """Power spectrum via FFT (no zero-padding)."""
    n = len(data)
    ps = np.abs(np.fft.fft(data)) ** 2
    return ps[: int(np.round(n / 2))]


def pSpectrum_padded(data):
    """Power spectrum via FFT with 2x zero-padding."""
    n = len(data)
    padded = np.zeros(2 * n, dtype=np.float64)
    padded[:n] = data
    ps = np.abs(np.fft.fft(padded)) ** 2
    return ps[:n]


def patterson(
    dataFile,
    save=True,
    maxZ=2000,
    qMin=0.03,
    qMax=0.2,
    rrfNorm=True,
    zCutOff=70,
    rho_sub=0.71,
    rho_pre=0.0,
    meanSub=False,
    windowFunc=False,
    padding=False,
    norm=True,
    plot=False,
    logPlot=True,
    scale=1,
    method=0,
    modo=True,
    outputFile="foo",
    ax=None,
):
    """
    Compute and optionally plot the Patterson function (Fourier transform of
    the normalized XRR reflectivity).

    Parameters
    ----------
    dataFile : str
        Path to a two-column whitespace-delimited text file (q, intensity).
    save : bool
        Write the result to disk (.fft for method=0, .ftMan for method=1).
    maxZ : int
        Maximum real-space range in Angstrom.
    qMin, qMax : float
        q-range (Å⁻¹) used for the transform.
    rrfNorm : bool
        If True, normalize by Fresnel reflectivity (R/R_F); otherwise use R·q⁴.
    zCutOff : float
        Distances below this value (Å) are set to NaN (removes near-origin
        artefacts).
    rho_sub, rho_pre : float
        Electron densities of substrate and pre-layer (for Fresnel calculation).
    meanSub : bool
        Subtract mean of R/R_F before transform.
    windowFunc : bool
        Apply a Bartlett (triangular) window before transform.
    padding : bool
        Zero-pad the data to 2× length before FFT (method=0 only).
    norm : bool
        Normalize the output to its maximum value.
    plot : bool
        Draw a linear-scale plot on *ax* (or current axes).
    logPlot : bool
        Draw a log-scale plot on *ax* (or current axes).
    scale : float
        Multiplicative scaling factor for plotting.
    method : {0, 1}
        0 = FFT (fast), 1 = manual DFT (slow, educational).
    modo : bool
        If True, return |FT|² (method=1 only); otherwise return complex FT.
    outputFile : str
        Base name for the saved output file.
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on. Uses current axes if None.

    Returns
    -------
    distance : ndarray
        Real-space distances (Å).
    patterson_out : ndarray
        Patterson function values.
    """
    import matplotlib.pyplot as plt

    data = np.loadtxt(dataFile)
    qq = data[:, 0]
    ii = data[:, 1]

    # ------------------------------------------------------------------ FFT
    if method == 0:
        qq_interp = np.arange(np.min(qq), np.max(qq), np.pi / (maxZ * 10))
        ii_interp = np.interp(qq_interp, qq, ii)
        qq, ii = qq_interp, ii_interp

        mask = (qq > qMin) & (qq < qMax)
        qq, ii = qq[mask], ii[mask]

        rrf = _normalize_intensity(qq, ii, rrfNorm, rho_sub, rho_pre)
        if meanSub:
            rrf -= np.mean(rrf)

        w = np.bartlett(len(rrf)) if windowFunc else np.ones(len(rrf))
        rrf_w = rrf * w

        ps = pSpectrum_padded(rrf_w) if padding else pSpectrum_nonPadded(rrf_w)
        freq = np.linspace(0, np.pi / (qq[3] - qq[2]), len(ps))

        cutoff = int(len(ps) * 0.1)
        distance = freq[:cutoff]
        patterson_out = ps[:cutoff]

    # ----------------------------------------------------------- Manual DFT
    elif method == 1:
        mask = (qq > qMin) & (qq < qMax)
        qq, ii = qq[mask], ii[mask]
        dq = np.diff(qq)

        rrf = _normalize_intensity(qq, ii, rrfNorm, rho_sub, rho_pre)
        if meanSub:
            rrf -= np.mean(rrf)

        w = np.bartlett(len(rrf)) if windowFunc else np.ones(len(rrf))
        rrf_w = rrf * w

        zz = np.arange(0, maxZ, 1)
        patt = np.zeros(len(zz), dtype=np.complex64)

        # zip over N-1 points so dq indices are valid
        for iz, z in enumerate(zz):
            for q_j, rrf_j, dq_j in zip(qq[:-1], rrf_w[:-1], dq):
                patt[iz] += rrf_j * np.exp(1j * z * q_j) * dq_j

        patterson_out = np.abs(patt) ** 2 if modo else patt
        distance = zz

    else:
        raise ValueError(f"method must be 0 (FFT) or 1 (manual DFT), got {method}")

    # ----------------------------------------------- post-processing & plot
    if norm:
        patterson_out = patterson_out / np.amax(np.abs(patterson_out))

    patterson_out = np.real(patterson_out).astype(float)
    patterson_out[distance < zCutOff] = np.nan

    target_ax = ax if ax is not None else plt.gca()
    if plot:
        target_ax.plot(distance, patterson_out * scale, label=dataFile)
    if logPlot:
        target_ax.semilogy(distance, patterson_out * scale, label=dataFile)
    if plot or logPlot:
        target_ax.legend()

    if save:
        ext = "fft" if method == 0 else "ftMan"
        np.savetxt(f"{outputFile}.{ext}", np.column_stack((distance, patterson_out)))
        print(f"Saved: {outputFile}.{ext}")

    return distance, patterson_out
