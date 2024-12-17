# HIPSCutoutViewer

A tool for viewing astronomical image cutouts using HiPS (Hierarchical Progressive Survey).

Project repository: https://github.com/yogeshw/HIPSCutoutViewer

## Features

- Query astronomical objects by name (using Simbad) or direct coordinates (RA/Dec)
- Default object: M 51
- Default surveys automatically loaded:
  - CDS/P/2MASS/color
  - CDS/P/HST/EPO
  - CDS/P/SDSS9/color
- Support for additional astronomical surveys including:
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
  - Save image collages as PNG
  - Download FITS format files with WCS for scientific analysis

## Requirements

- Python 3.x
- PyQt6
- astropy
- astroquery
- matplotlib
- requests

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yogeshw/HIPSCutoutViewer.git
cd HIPSCutoutViewer
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Run the application:
```bash
python hips_cutout_viewer.py
```

## Technical Support and Contributing
For technical issues, feature requests, or contributions, please visit the [project repository](https://github.com/yogeshw/HIPSCutoutViewer) and submit an issue or pull request.

## License
This project is licensed under the GPL 3 License
Copyright (c) 2024 Yogesh Wadadekar
