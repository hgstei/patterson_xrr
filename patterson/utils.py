"""
Utility functions for file handling and batch processing.
"""

import glob
import re

from .core import patterson


def atof(text):
    try:
        return float(text)
    except ValueError:
        return text


def natural_keys(text):
    """
    Sort key for human-friendly ordering of filenames containing numbers.

    Usage: ``alist.sort(key=natural_keys)``
    """
    return [
        atof(c)
        for c in re.split(r"[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)", text)
    ]


def get_all_ft(qMax=0.6, qMin=0.15, maxZ=500, pattern="*.xrr"):
    """
    Run the manual DFT Patterson function on all matching files in the current
    working directory.

    Parameters
    ----------
    qMax, qMin : float
        q-range (Å⁻¹) passed to :func:`patterson`.
    maxZ : int
        Maximum real-space range in Å.
    pattern : str
        Glob pattern for input files (default ``*.xrr``).
    """
    fileList = sorted(glob.glob(pattern), key=natural_keys)
    for fileName in fileList:
        patterson(
            fileName,
            scale=1.0,
            qMax=qMax,
            qMin=qMin,
            maxZ=maxZ,
            zCutOff=0,
            method=1,
            norm=False,
            plot=False,
            logPlot=False,
            outputFile=f"{fileName}_{qMin}_{qMax}",
        )
