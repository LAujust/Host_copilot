<div align="center">

<img src="HC_logo.png" alt="Host Copilot logo" width="560">

# Host Copilot

**Automatic host-galaxy identification pipeline for astronomical transients**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![astropy](https://img.shields.io/badge/powered%20by-AstroPy-orange.svg)](https://www.astropy.org/)
[![DOI](https://img.shields.io/badge/REGALADE-10.1051%2F0004--6361%2F202556896-blue)](https://doi.org/10.1051/0004-6361/202556896)

<p align="center">
  <img src="https://img.shields.io/badge/status-beta-yellow" alt="Status: Beta">
  <img src="https://img.shields.io/badge/designed%20for-Einstein%20Probe-important" alt="Einstein Probe">
  <img src="https://img.shields.io/badge/catalogs-REGALADE%20%7C%20Pan--STARRS%20%7C%20LS%20%7C%20NED-green" alt="Catalogs">
</p>

</div>

---

## 📖 Overview

**Host Copilot** is an automated pipeline that identifies likely host galaxies for astronomical transients. Given a transient position (RA, Dec) and a search radius, it queries multiple galaxy catalogs and source catalogs to find candidate hosts, applies redshift cuts, and produces an interactive visualization.

The pipeline is designed for — but not limited to — transients from the **Einstein Probe (EP)** mission, including X-ray flashes, GRB afterglows, and other high-energy transients.

### Key Features

- 🔍 **Multi-catalog cone search** — Queries REGALADE (primary galaxy catalog), Pan-STARRS DR2, Legacy Survey DR10, and NED
- 📐 **Redshift filtering** — Filters candidates by user-defined redshift threshold
- 🗺️ **Interactive visualization** — Renders search results on an Aladin Lite sky viewer with galaxy ellipses
- 📸 **Image cutout retrieval** — Fetches Pan-STARRS and Legacy Survey FITS cutouts for the field
- 🚀 **Quick & full modes** — Quick mode uses REGALADE only; full mode spans all catalogs
- 📦 **Batch processing** — Process entire transient catalogs with summary statistics
- 🌐 **Standalone HTML output** — Generates portable, self-contained Aladin HTML pages for sharing

---

## 📦 Installation

```bash
# Clone the repository
git clone https://github.com/LAujust/host_copilot.git
cd host_copilot

# Install the package and runtime dependencies
pip install -e .

# Development/test dependencies (optional)
pip install -e '.[test]'
```

### Dependencies

| Package | Purpose |
|---------|---------|
| `astropy` | Coordinate handling, FITS I/O, unit conversions |
| `astroquery` | VizieR/Vizier catalog queries |
| `ipyaladin` | Interactive sky visualization (Jupyter) |
| `regions` | Sky region definitions |
| `requests` | HTTP queries for Pan-STARRS, NED |
| `pyvo` | TAP & SIA queries for Legacy Survey (NOIRLab) |
| `pandas` | Catalog manipulation |
| `photutils` | Image-only source detection and segmentation |
| `scipy` | Statistical and image-analysis support |
| `Pillow` | Image handling |
| `matplotlib` | Plotting |

---

## 🚀 Quick Start

### Basic Usage

```python
from host_copilot import HostPipeline

# Position of a transient (RA, Dec in degrees, radius in arcsec)
ra, dec, r = 337.069, 50.764, 180

# Initialize the pipeline
hostpipe = HostPipeline(ra, dec, r, zcutout=0.1, save_path='./results')

# Run in quick mode (REGALADE catalog only)
aladin, cat_table = hostpipe.run()

# `aladin` is an interactive ipyaladin widget
# `cat_table` is an Astropy Table of candidate galaxies
```

### Structured quick and full searches

The compatibility ``run()`` method above remains available.  New analysis code
should use the structured API, which exposes catalog provenance, partial-service
failures, geometry tests, redshift provenance, and auditable ranking terms:

```python
from host_copilot import (
    CircleLocalization,
    HostPipeline,
    SearchConfig,
    TransientContext,
)

transient = TransientContext(
    name="EP260321a",
    localization=CircleLocalization(
        ra_deg=149.9287,
        dec_deg=0.4177,
        radius_arcsec=10.0,
        confidence=0.90,
    ),
    optical_ra_deg=149.928704,
    optical_dec_deg=0.418445,
    redshift=0.0343,
)

pipeline = HostPipeline(localization=transient.localization, save_path="results")

# REGALADE only, z < 0.1, with a 30-second deadline.
quick = pipeline.search(
    transient,
    SearchConfig(mode="quick"),
)

# Multi-catalog search, z < 0.5 by default, including unknown-z and
# image-only candidates.  Successful partial results survive provider outages.
full = pipeline.search(
    transient,
    SearchConfig(mode="full"),
)

full.to_csv("results/host_candidates.csv")
print(full.provider_status)
print(full.warnings)
```

``relative_probability`` is a normalized, transparent ranking weight, not a
calibrated host probability.  Until calibration gates are satisfied,
``posterior_probability`` is deliberately unset and the pipeline does not emit
``secure`` or ``probable`` probability labels.

### Output Example

```
==================================================
[QUICK MODE]
Searching for galaxies...
Querying REGALADE catalog...
--------------------------------------------------
Found 3 galaxies within z < 0.1 and r = 180 arcsec
--------------------------------------------------
Galaxy: WISEA J222804.94+504346., z=0.054, sep=166.84"
Galaxy: LEDA 97246, z=0.015, sep=43.63"
Galaxy: 2MASX J22281677+5047231, z=0.026, sep=91.51"
HostPipeline run completed.
==================================================
```

---

## 🧬 Architecture

```
                  ┌─────────────────────────────────┐
                  │        HostPipeline              │
                  │ deadlines + partial-state output │
                  └──────────┬──────────────────────┘
                             │
              ┌──────────────┼──────────────┐
              ▼              ▼              ▼
    ┌──────────────────┐ ┌──────┐  ┌──────────────────┐
    │ Catalog providers │ │ Cache│  │ Image recovery   │
    │ (7 remote feeds)  │ │      │  │ + Photutils     │
    └────────┬─────────┘ └──────┘  └──────────────────┘
             │
    ┌────────┼────────┬───────┬──────────┐
    ▼        ▼        ▼       ▼          ▼
  Normalize → deduplicate → geometry/redshift/Gaia veto → ranking
                  → pandas / Astropy / CSV / Aladin
```

### Modules

| Module | Class / Function | Role |
|--------|-----------------|------|
| `pipeline.py` | `HostPipeline` | Orchestrates quick/full workflows, partial state, and Aladin compatibility |
| `models.py` | Localization, context, candidate, result models | Stable typed API and export records |
| `providers.py` | `CatalogProvider` implementations | Bounded queries and schema normalization for seven services |
| `cache.py` | `QueryCache` | Release-aware positive, negative, and stale-response cache |
| `association.py` | Deduplication and ranking functions | Geometry, redshift compatibility, Gaia veto, and score audit trail |
| `extraction.py` | `ImageRecovery` | Photutils detection of catalog-missed image sources |
| `catalog.py`, `image.py` | `GalaxyFinder`, `Imager` | Legacy API and cutout compatibility |

---

## 📚 Catalogs

### Galaxy Catalogs

| Catalog | Query Method | Description |
|---------|-------------|-------------|
| **REGALADE** | VizieR (`J/A+A/706/A284/regalade`) | Comprehensive galaxy morphology catalog with ellipticities, position angles, and redshifts. **Primary catalog for host identification.** [Paper](https://doi.org/10.1051/0004-6361/202556896) |
| **DESI DR1** | NOIRLab TAP | High-quality spectroscopic galaxy redshifts |
| **SDSS DR17** | NOIRLab TAP | Spectroscopic galaxy redshifts |
| **NED** | IPAC ConeSearch API | Cross-identifications and literature redshifts |

### Source Catalogs

| Catalog | Query Method | Description |
|---------|-------------|-------------|
| **Pan-STARRS DR2** | MAST API | Optical detections, morphology proxy, and northern imaging |
| **Legacy Survey DR10** | NOIRLab TAP + viewer cutouts | Tractor morphology, photometry, photo-z, and imaging |
| **Gaia DR3** | NOIRLab TAP | Secure stellar-contaminant veto from astrometry/classification |

### Pipeline Modes

- **Quick mode** — Queries REGALADE through bounded HTTP, uses `z_max=0.1`,
  caches positive and negative responses, and has a 30-second default deadline.
- **Full mode** — Concurrently queries REGALADE, Legacy Survey DR10 plus its
  photo-z table, Pan-STARRS DR2, DESI DR1, SDSS DR17, NED, and Gaia DR3.  It
  defaults to `z_max=0.5`, retains plausible galaxies without redshifts, runs
  optional image-source recovery, and returns provenance-rich partial results
  if one or more services fail.

Both modes search beyond the formal localization by a configurable host margin,
then separately report whether each galaxy center lies inside the localization,
whether its footprint overlaps, and its directional-light-radius offset.  A
returned catalog row is therefore not automatically a host association.

---

## 🧪 Batch Processing

The repository ships two batch-processing scripts for running the pipeline on entire transient catalogs:

### `examples/process_ep_transients.py`

Processes all EP transients from `EP_data/EP_transients.csv` and produces per-transient output directories with:

- `aladin.html` — Standalone interactive Aladin viewer (no Python needed!)
- `aladin.ipynb` — Jupyter notebook replicating the view
- `cat_table.csv` — Filtered candidate galaxies
- `summary.csv` — Aggregate statistics

```bash
python examples/process_ep_transients.py
```

### `examples/process_ep_identified_sources.py`

Batch processes the EP identified-source catalog with optional parallel workers:

```bash
# Serial processing
python examples/process_ep_identified_sources.py

# Parallel processing (4 workers)
python examples/process_ep_identified_sources.py --workers 4

# Test run with 20 sources
python examples/process_ep_identified_sources.py --workers 2 --limit 20
```

---

## 📊 Example Results

From a batch run of 71 EP transients:

| Metric | Value |
|--------|-------|
| Total transients processed | 71 |
| With host candidates (z < 0.1) | 17 (24%) |
| Without host candidates | 54 (76%) |

![Host separation distribution](https://img.shields.io/badge/sep%20range-1.7%E2%80%93188%20arcsec-lightgrey)
![Redshift range](https://img.shields.io/badge/z%20range-0.015%E2%80%930.098-lightgrey)

### Sample Detection: EP260321a

```
Nearest host: SDSS J095942.88+002506.2
Separation:   1.69"
Redshift:     z = 0.035
```

### Sample Detection: EP260227a

```
Nearest host: DESI J145224.47-112324.3
Separation:   16.0"
Redshift:     z = 0.062
```

---

## 📂 Project Structure

```
host_copilot/
├── host_copilot/              # Core Python package
│   ├── __init__.py            # Package entry point (v0.0.1)
│   ├── pipeline.py            # HostPipeline — main orchestrator
│   ├── catalog.py             # GalaxyFinder — multi-catalog queries
│   ├── image.py               # Imager — FITS cutout retrieval
│   └── utils.py               # Shared utilities / imports
├── examples/                  # Usage examples & batch scripts
│   ├── Test.ipynb             # Interactive demo notebook
│   ├── process_ep_transients.py        # EP transient batch processor
│   ├── process_ep_identified_sources.py # Source catalog batch processor
│   └── EP_data/               # Per-transient output directories
│       ├── summary.csv        # Aggregate batch results
│       ├── EP251214b/         # Aladin HTML + notebook + catalog
│       ├── EP260321a/         # ...
│       └── ...
├── EP_data/                   # Input transient catalogs
│   ├── EP_transients.csv      # EP transient list
│   └── EP_identified_source_list.csv
├── Catalogs/                  # Large catalog files (gitignored)
├── LICENSE                    # MIT License
└── README.md                  # This file
```

---

## 🛠️ API Reference

### `HostPipeline(ra, dec, r_arcsec, zcutout=None, quick=True, save_path='./', ...)`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ra` | `float` | — | Right Ascension in degrees (ICRS) |
| `dec` | `float` | — | Declination in degrees (ICRS) |
| `r_arcsec` | `float` | — | Search radius in arcseconds |
| `zcutout` | `float` | mode default | Maximum redshift: 0.1 quick, 0.5 full |
| `quick` | `bool` | `True` | If True, uses REGALADE only; otherwise queries all catalogs |
| `save_path` | `str` | `'./'` | Directory for cached catalogs and output |

**Returns:** `(Aladin widget, Astropy Table)` tuple.

For the new API, ``HostPipeline.search(context, config)`` returns a
``HostSearchResult`` with candidates, provider statuses, warnings, metadata,
and pandas/Astropy/CSV export methods.  Circular and elliptical localizations
are supported.  Existing EP radii are interpreted as 90% containment radii.

### `GalaxyFinder(ra, dec, r_arcsec, redo=False, save_path='./')`

| Method | Description |
|--------|-------------|
| `find_regalade()` | Query REGALADE via VizieR |
| `find_ps()` | Query Pan-STARRS DR2 via MAST |
| `find_ls()` | Query Legacy Survey DR10 via NOIRLab TAP |
| `find_ned()` | Query NED via IPAC ConeSearch |

### `Imager(ra, dec, r_arcsec, band='r', save_path=None)`

| Method | Description |
|--------|-------------|
| `get_cutout()` | Try PS1 first, fall back to LS DR10 |
| `PS_cutout()` | Fetch Pan-STARRS DR2 FITS cutout |
| `LS_cutout()` | Fetch Legacy Survey DR10 FITS cutout |

---

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **REGALADE** catalog and its authors ([Tranin et al. 2026](https://doi.org/10.1051/0004-6361/202556896))
- **Pan-STARRS** project, MAST archive, and the STScI
- **Legacy Survey** (NOIRLab) and the **DESI** collaboration
- **NED** (IPAC/Caltech)
- **CDS** for the Aladin Lite and VizieR services
- The **Einstein Probe** mission team

---

<div align="center">
  <sub>Built with ❤️ for transient astronomy</sub>
</div>
