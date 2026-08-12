![banner](Banner/CATCH_Banner_New.png)

[![Version](https://img.shields.io/badge/V-2.0-blue)]()
[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](https://www.python.org)
[![numpy](https://img.shields.io/badge/NumPy-4DABCF?logo=numpy&logoColor=fff)](http://www.numpy.org/)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Overview

`CATCH` is a Python3 program which allows an observer to rapidly search for and verify calibrator stars for interferometric observations carried out at the [Center for High Angular Resolution Astronomy (CHARA) Array](https://chara.gsu.edu). These include observations conducted with the **MIRC-X** (_H_), **MYSTIC** (_K_), **Silmaril** (_H_ & _K_), and **SPICA** (_R_) beam combiners. Future support for other instruments at the Array is under development. 

Viable stars need to meet certain requirements to work accurately as calibrators. These requirements are chosen to work for the majority of **MIRC-X/MYSTIC** use cases, but are easily modifiable within the source code to fit user-specific criteria. The default for the calibrators are:
- They need to be bright, i.e. V<sub>mag</sub> < 9.0 for all beam combiners
- H<sub>mag</sub> < 6.4 for MIRC-X/MYSTIC; H<sub>mag</sub> < 10.5 for Silmaril; RP<sub>mag</sub> < 5.4 for SPICA
- Their angular diameters need to be relatively small (UDD<sub>H</sub> < 0.4)
- They <b>MUST NOT</b> be binaries
- They <b>MUST NOT</b> be rapid rotators (V•sin(i) < 100 km/s)
- They shouldn't have significant IR excesses
- They ideally shouldn't have close field companions
- They must be within the CHARA Array's declination limit

Potential calibrators are checked against a set of published catalogues and parameters therein to verify that they meet these requirements. These catalogues include:
- [Jean-Marie Mariotti Center Stellar Diameters Catalogue v2 (JSDC)](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=II/346/jsdc_v2) to check brightness and angular diameters
- [Gaia DR3](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3) to check multiplicity and rapid rotation
- [Kervella et al. 2022](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A%2bA/657/A7) to check multiplicity
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
You can either run the code using `uv run catch.py` or with `python3 catch.py` for an interactive run. 

If you are a more experienced user, and would prefer a faster use-case, you can use `CATCH` as a CLI. You can run it as such using `uv run catch.py -[FLAGS]`. For a list of usable flags, use `uv run catch.py -h`.

### Example Usage — Generating List of Verified Calibrators
You can use the built-in `main` function (the default run case) as shown in the example below to generate a list of viable calibrators for a given science target. The calibrators are saved in an `ascii` file titled **star_name_Calibrators.txt**. 
```
$ uv run catch.py

        ############################ This is CATCH ############################
                   [C]HARA [A]rray's [T]hrifty [C]alibrator [H]unter
                                 Version: 2.0 | 2026/08/12                           
                       https://github.com/akshat-chaturvedi/CATCH                 
        #######################################################################
        
What instrument is this for? 
Type HK for MIRC-X and MYSTIC; type S for Silmaril; type R for Spica
HK
Vizier server: vizier.cds.unistra.fr
vizier.cds.unistra.fr server up (HTML Response Code: 200)
Would you like to find calibrators for a science target (type A), or check a possible calibrator's viability (type B)?:
A
Please enter the name of your target (please ensure the name is resolvable in SIMBAD):
HD191610
Would you like to filter calibrators by whether they have close companions in Gaia DR3 Y/[N]?
Y
Please enter the desired cutoff radius (in arcseconds) for Gaia companions:
10
Beginning calibration search for target: HD191610
-->Querying JMMC Stellar Diameters Catalogue (JSDC)...
-->Query complete!
-->Querying Gaia DR3 Catalogue...
-->Query complete!
-->Checking for close Gaia companions within 10.0"
-->Querying Kervella et al. 2022 Catalogue...
-->Query complete!
-->Querying Cruzalebes et al. 2019 Catalogue (MDFC)...
-->Query complete!
Found 11 viable calibrators in 91.37 seconds!
```

### Example Usage — Verifying a Provided Potential Calibrator
You can use the built-in `main` function (the default run case) as shown in the example below to generate a list of calibrators for a given science target.
```
$ uv run catch.py

        ############################ This is CATCH ############################
                   [C]HARA [A]rray's [T]hrifty [C]alibrator [H]unter
                                 Version: 2.0 | 2026/08/12                           
                       https://github.com/akshat-chaturvedi/CATCH                 
        #######################################################################
        
What instrument is this for? 
Type HK for MIRC-X and MYSTIC; type S for Silmaril; type R for Spica
HK
Vizier server: vizier.cds.unistra.fr
vizier.cds.unistra.fr server up (HTML Response Code: 200)
Would you like to find calibrators for a science target (type A), or check a possible calibrator's viability (type B)?:
B
Would you like to check a single calibrator, or multiple calibrators? [S]/M
S
Please enter the name of your calibrator (please ensure the name is resolvable in SIMBAD):
HD190009
Would you like to filter calibrators by whether it has a companion within 10" in Gaia DR3 Y/[N]?
Y
Checking calibrator viability of: HD190009
-->Querying JMMC Stellar Diameters Catalogue (JSDC)...
-->Query complete!
---->Query results: V mag = 6.44, H mag = 5.32, UDDH = 0.365
-->HD190009 passes JMMC Stellar Diameters Catalogue (JSDC) checks!
-->Querying Gaia DR3 Catalogue...
-->Query complete!
---->Query results: IPDfmp = 0.0, RUWE = 1.08, Vbroad = 7.09 km/s
-->HD190009 passes Gaia DR3 Catalogue checks!
-->Querying Kervella et al. 2022 Catalogue...
-->Query complete!
---->Query results: DMS = 0, W = 0, BinHG1 = 0, BinH2G2 = 0, BinH2EG3b = 0, snrPMaHG1 = 1.94, snrPMaH2G2 = 1.92, snrPMaH2EG3b = 0.99
-->HD190009 passes Kervella et al. 2022 Catalogue checks!
-->Querying Cruzalebes et al. 2019 Catalogue (MDFC)...
-->Query complete!
---->Query results: CalFlag = 0, IRflag = 0
-->HD190009 passes Cruzalebes et al. 2019 Catalogue (MDFC) checks!
-->HD190009 passed 5/5 checks
Confirmed HD190009 is likely an ideal calibrator in 4.87 seconds!
```

### Example Usage — Generating List of Verified Calibrators Using the CLI Mode
You can use `CATCH` in the CLI mode to find a list of viable calibrators for a given science target as shown in the example below
```
$ uv run catch.py -f -g -t "HD89484"

        ############################ This is CATCH ############################
                   [C]HARA [A]rray's [T]hrifty [C]alibrator [H]unter
                                 Version: 1.4 | 2026/03/20                           
                       https://github.com/akshat-chaturvedi/CATCH                 
        #######################################################################
        
Vizier server: vizier.cds.unistra.fr
vizier.cds.unistra.fr server up (HTML Response Code: 200)
Beginning calibration search for target: HD89484
-->Querying JMMC Stellar Diameters Catalogue (JSDC)...
-->Query complete!
-->Querying Gaia DR3 Catalogue...
-->Query complete!
-->Checking for close Gaia companions within 5"
-->Querying Kervella et al. 2019 Catalogue...
-->Query complete!
-->Querying Cruzalebes et al. 2019 Catalogue (MDFC)...
-->Query complete!
Found 33 viable calibrators in 1.37 seconds!
```

### Known Issues
#### Vizier server timeouts
Sometimes (more so recently) the default Vizier site (https://vizier.cds.unistra.fr/) is down. In such a case, the Vizier 
query will either time out or simply not go through. At present, there doesn't seem to be a way to resolve this, so we 
recommend finding calibrators well in advance of observing, and not relying on `CATCH` for "on-the-fly" calibrator 
searches during observing runs.

#### Reliability of Gaia flags
The two parameters that `CATCH` checks for in the Gaia DR3 database is the **Renormalised Unit Weight Error (RUWE)** 
and the **Image Parameter Determination fraction of multiple peaks (IPDfmp)**. **RUWE** is a measure of the goodness of fit
of a single star astrometric model to each Gaia source. While **RUWE** is a very good indicator of the possible multiplicity 
of a source, it usually peaks at separations between about 0.04 to 1 arcsecond, so anything closer or wider won't be 
picked up by it.
**IPDfmp** measures the fraction of Gaia images that have more than one peak in the PSF. This parameter is pretty 
reliable, and usually IPDfmp > 2 is indicative of a companion. It is however less useful for very close companions though,
as it peaks between ~ 0.1 to 1 arcsecond separations. As such, a possible calibrator that passes the Gaia checks may still
be a very close separation binary star.

#### Silmaril Calibrators
The calibrator constraints on Silmaril can make it difficult to find a large sample of calibrators. If this issue arises for 
your science target, increase the vizier row limit on **line 88** in `catch_silmaril.py`. If the issue still persists, it is likely
you will need to change the constraints.

### Dependencies

The latest version of `CATCH` has been developed for:
- python>=3.10
- astropy>=5.2.2
- astroquery>=0.4.7
- numpy>=1.24.4

## Acknowledgments
- If you used `CATCH` to help with your observations, please consider starring the GitHub repository and including the following acknowledgement in your publications: 
  - _This work made use of CATCH, CHARA Array's Thrifty Calibrator Hunter, available at https://github.com/akshat-chaturvedi/CATCH_
- We would like to thank Becky Flores, Dr. Gail Schaefer, and Dr. Cyprien Lanthermann at the CHARA Array for their guidance in defining the necessary criteria for good calibrators for MIRC-X/MYSTIC and Silmaril (and all the observing help!)
- We would like to thank Dr. Denis Mourard at Observatoire de la Côte d’Azur for his guidance in defining the necessary criteria for good calibrators for SPICA
- Banner background image: [The exotic stellar population of Westerlund 1 (ESA/Webb, NASA & CSA, M. Zamani (ESA/Webb), M. G. Guarcello (INAF-OAPA) and the EWOCS team)](https://esawebb.org/images/potm2409a/)
