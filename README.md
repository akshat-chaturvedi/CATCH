![banner](Banner/CATCH_banner.png)

[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](https://www.python.org)
[![numpy](https://img.shields.io/badge/NumPy-4DABCF?logo=numpy&logoColor=fff)](http://www.numpy.org/)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Overview

`CATCH` is a Python3 program which allows an observer to rapidly search for and verify calibrator stars for interferometric observations carried out at the [Center for High Angular Resolution Astronomy (CHARA) Array](https://chara.gsu.edu) with the **MIRC-X** (_H_-band) and **MYSTIC** (_K_-band) instruments. Future support for other instruments at the Array is under development. 

Viable stars need to meet certain requirements to work accurately as calibrators. These requirements are chosen to work for the majority of **MIRC-X/MYSTIC** use cases, but are easily modifiable within the source code to fit user-specific criteria. The default for the calibrators are:
- They need to be bright (i.e. V<sub>mag</sub> < 9.0, H<sub>mag</sub> < 6.4)
- Their angular diameters need to be relatively small (UDD<sub>H</sub> < 0.4)
- They <b>MUST NOT</b> be binaries
- They <b>MUST NOT</b> be rapid rotators (V•sin(i) < 100 km/s)
- They shouldn't have significant IR excesses
- They ideally shouldn't have close field companions
- They must be within the CHARA Array's declination limit

Potential calibrators are checked against a set of published catalogues and parameters therein to verify that they meet these requirements. These catalogues include:
- [Jean-Marie Mariotti Center Stellar Diameters Catalogue v2 (JSDC)](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=II/346/jsdc_v2) to check brightness and angular diameters
- [Gaia DR3](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3) to check multiplicity and rapid rotation
- [Kervella et al. 2019](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/623/A72) to check multiplicity
- [Cruzalebes et al. 2019 (MDFC) ](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/361) to check infrared excess

`CATCH` has two main capabilities. It can generate a list of viable, cross-checked calibrator stars for a given science target, and it can verify a potential calibrator — or list thereof — independently chosen by the user.

### Install with GitHub and UV

The GitHub version is the most up-to-date. To install from the GitHub repository, first install `uv` using:
```
$ pip install uv
```
You can then clone the repository and install `CATCH` using:
```
$ git clone https://github.com/akshat-chaturvedi/CATCH.git
``` 
Then `cd` into the `CATCH` directory and use `uv sync` to ensure all dependencies match 
``` 
$ cd CATCH && uv sync
``` 

**NOTE:** `CATCH` requires Python version >=3.10. If you do not have this installed on your device, you can use `uv` to update the virtual environment to a compatible python version with the following:
```
$ uv venv --python 3.10
```

### Running the code
You can either run the code using `uv run catch.py` or with `python3 catch.py`.

### Example Usage — Generating List of Verified Calibrators
You can use the built-in `main` function (the default run case) as shown in the example below to generate a list of calibrators for a given science target. The calibrators are saved in an `ascii` file titled **star_name_Calibrators.txt**. 
```
$ uv run catch.py

        ############################ This is CATCH ############################
                   [C]HARA [A]rray's [T]hrifty [C]alibrator [H]unter
                                 Version: 1.1 | 2025/08/27                           
                       https://github.com/akshat-chaturvedi/CATCH                 
        #######################################################################
        
Would you like to find calibrators for a science target (type A), or check a possible calibrator's viability (type B)?:
A
Please enter the name of your target (please ensure the name is resolvable in SIMBAD):
HD29441
Would you like to filter calibrators by whether they have close companions in Gaia DR3 Y/[N]?
Y
Please enter the desired cutoff radius (in arcseconds) for Gaia companions:
10
Beginning calibration search for target: HD29441
-->Querying JMMC Stellar Diameters Catalogue (JSDC)...
-->Query complete!
-->Querying Gaia DR3 Catalogue...
-->Query complete!
-->Checking for close Gaia companions within 10.0"
-->Querying Kervella et al. 2019 Catalogue...
-->Query complete!
-->Querying Cruzalebes et al. 2019 Catalogue (MDFC)...
-->Query complete!
Found 17 viable calibrators in 3.41 seconds!
```

### Example Usage — Verifying a Provided Potential Calibrator
You can use the built-in `main` function (the default run case) as shown in the example below to generate a list of calibrators for a given science target.
```
$ uv run catch.py

        ############################ This is CATCH ############################
                   [C]HARA [A]rray's [T]hrifty [C]alibrator [H]unter
                                 Version: 1.1 | 2025/08/27                           
                       https://github.com/akshat-chaturvedi/CATCH                 
        #######################################################################
        
Would you like to find calibrators for a science target (type A), or check a possible calibrator's viability (type B)?:
B
Would you like to check a single calibrator, or multiple calibrators? [S]/M
S
Please enter the name of your calibrator (please ensure the name is resolvable in SIMBAD):
HD 30913
Would you like to filter calibrators by whether it has a companion within 10" in Gaia DR3 Y/[N]?
Y
Checking calibrator viability of: HD 30913
-->Querying JMMC Stellar Diameters Catalogue (JSDC)...
-->Query complete!
---->Query results: V mag = 6.82, H mag = 5.87, UDDH = 0.270
-->HD 30913 passes JMMC Stellar Diameters Catalogue (JSDC) checks!
-->Querying Gaia DR3 Catalogue...
-->Query complete!
---->Query results: IPDfmp = 0.0, RUWE = 1.32, Vbroad = 97.82 km/s
-->HD 30913 passes Gaia DR3 Catalogue checks!
-->Querying Kervella et al. 2019 Catalogue...
-->Query complete!
---->Query results: DMS = 0, W = 0, BinH = 0, BinG2 = 0
-->HD 30913 passes Kervella et al. 2019 Catalogue checks!
-->Querying Cruzalebes et al. 2019 Catalogue (MDFC)...
-->Query complete!
---->Query results: CalFlag = 0, IRflag = 0
-->HD 30913 passes Cruzalebes et al. 2019 Catalogue (MDFC) checks!
-->HD 30913 passed 5/5 checks
Confirmed HD 30913 is likely an ideal calibrator in 9.4 seconds!
```
### Known Issues
Sometimes (more so recently) the default Vizier site (https://vizier.cds.unistra.fr/) is down. In such a case, the Vizier query will time out. If this happens you can change which Vizier server is being used by uncommenting the following line from the `main` function:
```commandline
# If you get a timeout error, un-comment the following line to use the Japanese mirror site for Vizier
# conf.server = "vizier.nao.ac.jp"
print("Current Vizier server:", conf.server)
```

### Dependencies

The latest version of `CATCH` has been developed for:
- python>=3.10
- astropy>=5.2.2
- astroquery>=0.4.7
- numpy>=1.24.4

## Acknowledgments
- We would like to thank Becky Flores and Dr. Gail Schaefer at the CHARA Array for their expert guidance in defining the necessary criteria for calibrator determination (and all the observing help!)
- Banner background image: [The exotic stellar population of Westerlund 1 (ESA/Webb, NASA & CSA, M. Zamani (ESA/Webb), M. G. Guarcello (INAF-OAPA) and the EWOCS team)](https://esawebb.org/images/potm2409a/)
