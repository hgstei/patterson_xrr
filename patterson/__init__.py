"""
patterson — XRR Patterson function analysis toolkit.
"""

import os

from .core import patterson, pSpectrum_nonPadded, pSpectrum_padded
from .plot import plot_all, plotXrr
from .utils import get_all_ft, natural_keys

# Path to the bundled example dataset (IZO thin film on Si substrate)
example_data = os.path.join(os.path.dirname(__file__), "example_data.xrr")

__all__ = [
    "patterson",
    "pSpectrum_nonPadded",
    "pSpectrum_padded",
    "plotXrr",
    "plot_all",
    "get_all_ft",
    "natural_keys",
    "example_data",
]
