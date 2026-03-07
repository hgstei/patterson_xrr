# patterson

XRR Patterson function (Fourier transform) analysis toolkit.

Computes the Patterson function — the power spectrum of the normalized X-ray
reflectivity — to extract real-space layer distances from XRR data.

## Installation

```bash
pip install git+https://github.com/hgstei/patterson_xrr.git
```

Or clone and install locally:

```bash
git clone https://github.com/hgstei/patterson_xrr.git
cd patterson_xrr
pip install -e .
```

## Quick start

```python
from patterson import patterson, get_all_ft, plotXrr, plot_all

# Single file — FFT method
distance, patt = patterson(
    "sample.xrr",
    qMin=0.03, qMax=0.2,
    method=0,
    outputFile="sample",
)

# Batch process all *.xrr files in the current directory
get_all_ft(qMin=0.15, qMax=0.6, maxZ=500)

# Side-by-side XRR + Patterson plot of all files
fig, axes = plot_all(outputFile="overview", outputDir="figures")
```

## API

### `patterson(dataFile, ...)`

| Parameter   | Default | Description |
|-------------|---------|-------------|
| `dataFile`  | —       | Path to two-column (q, I) text file |
| `method`    | `0`     | `0` = FFT (fast), `1` = manual DFT |
| `qMin/qMax` | 0.03/0.2| q-range (Å⁻¹) for the transform |
| `maxZ`      | `2000`  | Maximum real-space range (Å) |
| `rrfNorm`   | `True`  | Normalize by Fresnel reflectivity |
| `zCutOff`   | `70`    | Mask distances below this value (Å) |
| `rho_sub`   | `0.71`  | Substrate electron density |
| `rho_pre`   | `0.0`   | Pre-layer electron density |
| `norm`      | `True`  | Normalize output to maximum |
| `save`      | `True`  | Write result to disk |
| `outputFile`| `"foo"` | Base name for output file |
| `plot`      | `False` | Linear plot |
| `logPlot`   | `True`  | Log-scale plot |
| `ax`        | `None`  | Target matplotlib axes |

Returns `(distance, patterson_values)` as NumPy arrays.

### `get_all_ft(qMax, qMin, maxZ, pattern)`

Batch-processes all files matching `pattern` (default `*.xrr`) in the current
directory using the manual DFT method.

### `plotXrr(dataFile, ...)`

Plot raw or Fresnel-normalized XRR data.

### `plot_all(...)`

Side-by-side XRR + Patterson panel for all `*.xrr` / `*.ftMan` files in the
current directory. Returns `(fig, (ax1, ax2))`.

## Output files

| Extension | Method | Contents |
|-----------|--------|----------|
| `.fft`    | 0      | `distance  patterson` |
| `.ftMan`  | 1      | `distance  patterson` |
