[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](https://www.python.org)
[![numpy](https://img.shields.io/badge/NumPy-4DABCF?logo=numpy&logoColor=fff)](http://www.numpy.org/)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Overview

`CaliFinder` is a Python3 module which allows an observer to rapidly search for and verify calibrator stars for interferometric observations carried out at the [Center for High Angular Resolution Astronomy (CHARA) Array](https://chara.gsu.edu).

Viable stars need to meet certain requirements to work accurately as calibrators:
- They need to be bright (i.e. V<sub>mag</sub> < 9.0, 3.9 < H<sub>mag</sub> < 6.4)
- Their angular diameters need to be relatively small (UDD<sub>H</sub> < 0.4)
- They <b>MUST NOT</b> be binaries
- They <b>MUST NOT</b> be rapid rotators
- They shouldn't have significant IR excesses
- They must be within the CHARA Array's declination limit

Potential calibrators are checked against a set of published catalogues and parameters therein to verify that they meet these requirements. These catalogues include:
- [Jean-Marie Mariotti Center Stellar Diameters Catalogue v2](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=II/346/jsdc_v2) to check brightness and angular diameters
- [Gaia DR3](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3) to check binarity and rapid rotation
- [Kervella et al. 2019](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/623/A72) to check binarity
- [Cruzalebes et al. 2019](vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/361) to check infrared excess

### Install with GitHub and UV

The GitHub version is the most up-to-date. To install from the GitHub repository, first install uv using:
```
pip install uv
```
You can then clone the repository and isntall CaliFinder using:
```
git clone https://github.com/akshat-chaturvedi/CaliFinder.git
uv sync
``` 

### Dependencies

The latest version of `CaliFinder` has been developed with:
- astropy>=5.2.2
- astroquery>=0.4.7
- numpy>=1.24.4
