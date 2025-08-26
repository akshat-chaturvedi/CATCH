![banner](Banner/CATCH_banner.png)

[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](https://www.python.org)
[![numpy](https://img.shields.io/badge/NumPy-4DABCF?logo=numpy&logoColor=fff)](http://www.numpy.org/)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Overview

`CATCH` is a Python3 module which allows an observer to rapidly search for and verify calibrator stars for interferometric observations carried out at the [Center for High Angular Resolution Astronomy (CHARA) Array](https://chara.gsu.edu).

Viable stars need to meet certain requirements to work accurately as calibrators:
- They need to be bright (i.e. V<sub>mag</sub> < 9.0, 3.9 < H<sub>mag</sub> < 6.4)
- Their angular diameters need to be relatively small (UDD<sub>H</sub> < 0.4)
- They <b>MUST NOT</b> be binaries
- They <b>MUST NOT</b> be rapid rotators
- They shouldn't have significant IR excesses
- They ideally shouldn't have close field companions
- They must be within the CHARA Array's declination limit

Potential calibrators are checked against a set of published catalogues and parameters therein to verify that they meet these requirements. These catalogues include:
- [Jean-Marie Mariotti Center Stellar Diameters Catalogue v2](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=II/346/jsdc_v2) to check brightness and angular diameters
- [Gaia DR3](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3) to check binarity and rapid rotation
- [Kervella et al. 2019](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/623/A72) to check binarity
- [Cruzalebes et al. 2019](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/361) to check infrared excess

CaliFinder has two main capabilities. It can generate a list of viable, cross-checked calibrator stars for a given science target, and it can verify a potential calibrator independently chosen by the user.

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

### Running the code
You can either run the code using `uv run califinder.py` or with `python3 califinder.py`.

### Example Usage — Generating List of Verified Calibrators
You can use the built-in `main` function (the default run case) as shown in the example below to generate a list of calibrators for a given science target. The calibrators are saved in an `ascii` file titled **star_name_Calibrators.txt**. 
```
$ uv run califinder.py

        ############################ This is CATCH ############################
                      [C]HARA [A]rray's [T]hrifty [C]alibrator [H]unter
                                 Version: 1.0 | 2025/08/26                           
                       https://github.com/akshat-chaturvedi/CaliFinder                 
        #######################################################################
        
Would you like to find calibrators for a science target (type A), or check a possible calibrators viability (type B)?:
A
Please enter the name of your target (preferably its HD number):
HD29441
Would you like to filter calibrators by whether they have close companions in Gaia DR3 Y/[N]?
Y
Please enter the desired cutoff radius for Gaia companions:
10
Beginning calibration search for target: HD29441
-->Querying JMMC Stellar Diameters Catalogue...
-->Query complete!
-->Querying Gaia DR3 Catalogue...
-->Query complete!
-->Checking for close Gaia companions within 10.0"
-->Querying Kervella et al. 2019 Catalogue...
-->Query complete!
-->Querying Cruzalebes et al. 2019 Catalogue...
-->Query complete!
Found 17 viable calibrators in 3.99 seconds!
```

### Example Usage — Verifying a Provided Potential Calibrator
You can use the built-in `main` function (the default run case) as shown in the example below to generate a list of calibrators for a given science target.
```
$ uv run califinder.py

        ############################ This is CATCH ############################
                      [C]HARA [A]rray's [T]hrifty [C]alibrator [H]unter
                                 Version: 1.0 | 2025/08/26                           
                       https://github.com/akshat-chaturvedi/CaliFinder                 
        #######################################################################
        
Would you like to find calibrators for a science target (type A), or check a possible calibrators viability (type B)?:
B
Please enter the name of your calibrator (preferably its HD number):
HD  30913                                 
Would you like to filter calibrators by whether it has a companion within 5" in Gaia DR3 Y/[N]?
Y
Checking calibrator viability of: HD  30913
-->Querying JMMC Stellar Diameters Catalogue...
-->Query complete!
-->Querying Gaia DR3 Catalogue...
-->Query complete!
-->Querying Kervella et al. 2019 Catalogue...
-->Query complete!
-->Querying Cruzalebes et al. 2019 Catalogue...
-->Query complete!
-->HD  30913 passed 5/5 checks
Confirmed HD  30913 is a viable calibrator in 0.84 seconds!
```

### Dependencies

The latest version of `CaliFinder` has been developed for:
- astropy>=5.2.2
- astroquery>=0.4.7
- numpy>=1.24.4

## Acknowledgments
Banner background image: [The exotic stellar population of Westerlund 1 (ESA/Webb, NASA & CSA, M. Zamani (ESA/Webb), M. G. Guarcello (INAF-OAPA) and the EWOCS team)](https://esawebb.org/images/potm2409a/)
