# HIPSCutoutViewer

A graphical viewer for astronomical image cutouts using the HiPS (Hierarchical Progressive Survey) system from CDS (Centre de Données astronomiques de Strasbourg).

## Features

- Query astronomical objects by name (using Simbad) or direct coordinates (RA/Dec)
- Support for multiple astronomical surveys including:
  - 2MASS (J, H, K bands)
  - PanSTARRS DR1 (g, r, i, z, y bands)
  - DSS2 (blue, red)
  - And many other CDS HiPS surveys
- Customizable field of view (cutout size)
- Interactive survey selection
- Image collage generation with:
  - North-East orientation markers
  - Scale bar (1 arcminute)
  - Survey labels
- Save capabilities:
  - Save image collages as JPEG
  - Download FITS format files with WCS for scientific analysis

## Requirements

- Python 3.x
- PyQt6
- astropy
- astroquery
- matplotlib
- requests

## Installation

1. Clone the repository
2. Install required packages:
