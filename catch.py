from astroquery.vizier import Vizier, conf
from astroquery.simbad import Simbad
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import collections
import numpy as np
import warnings
from astroquery.exceptions import NoResultsWarning

warnings.simplefilter("ignore", NoResultsWarning)

RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
MAGENTA = '\033[95m'
ORANGE = '\033[38;2;255;128;0m'
RESET = '\033[0m'
ITALIC = '\033[3m'
BLINK = '\033[5m'

# __version__ = '1.0 | 2025/08/26' # First release version :)
__version__ = '1.1 | 2025/08/27' # Changes to file output and printout formats


def cal_finder(star_name: str, gaia_comp_check: int | float | None = None) -> None:
    """
    Finds viable calibrator stars within 10 degrees for CHARA Array interferometric targets. Successful calibrators pass
    magnitude and diameter checks from the JMMC Stellar Diameters Catalogue
    (https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=II/346/jsdc_v2) and binarity checks from the Gaia DR3
    (https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3), Kervella et al. 2019
    (https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/623/A72), and Cruzalebes et al. 2019
    (https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/361) catalogues.

    :param star_name: The name of the target star
    :param gaia_comp_check: Flag that checks for nearby Gaia DR3 companions. Default = None. Enter a radius if you want
    to check calibrator viability based on nearby Gaia DR3 companions.

    :return: None. Saves a file containing calibrator information.
    """
    t1 = time.perf_counter()

    Simbad.add_votable_fields("V")

    star_details_table = Simbad.query_object(f"{star_name}")
    if len(star_details_table) == 0:
        exit(f"Error: {YELLOW}{star_name}{RESET} not found in SIMBAD. Please check that you typed it in correctly!")
    else:
        star_ra = star_details_table['ra'].value[0]
        star_dec = star_details_table['dec'].value[0]
        star_v_mag = star_details_table['V'].value[0]

        if star_dec <= -25:
            exit("This star is outside the declination limits (DEC > -25) for the CHARA Array!")

        star_coords = SkyCoord(ra=star_ra * u.deg, dec=star_dec * u.deg, frame='icrs')
        star_ra = star_coords.ra.to_string(unit=u.hour, sep=' ', pad=True, precision=2)
        star_dec = star_coords.dec.to_string(unit=u.deg, sep=' ', alwayssign=True, pad=True, precision=2)

    print(f"Beginning calibration search for target: {YELLOW}{star_name}{RESET}")
    # Check with JMMC Stellar Diameters Catalogue (Vmag < 9.0, Hmag < 6.4, UDDH < 0.4)
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "Name", "SpType", "Vmag", "Rmag","Hmag", "Kmag", "UDDH", "UDDK",
                             "+_r"], catalog="II/346/jsdc_v2")

    # The default vizier query row limit is set here to 100. If you would like to search for more, increase this number
    # NOTE: Increasing it will increase run-time
    vizier.ROW_LIMIT = 100
    print(f"-->Querying {BLUE}JMMC Stellar Diameters Catalogue (JSDC){RESET}...")

    # By default, CATCH queries the JMMC catalog for stars within 10 degrees, but it can be increased as required. The
    # default constraints are described in the README file, but can be edited by changing the column_filters parameter
    # in the query below. Guidance on syntax can be found at https://vizier.cds.unistra.fr/vizier/vizHelp/cst.htx
    jmmc_result = vizier.query_region(f"{star_name}", radius="10d", column_filters={"Vmag":"<9.0", "Hmag":"<=6.4",
                                                                           "UDDH": "<0.4", "_DEJ2000": ">-25"})
    print(f"-->{GREEN}Query complete!{RESET}")
    if len(jmmc_result) > 0:
        jmmc_result = jmmc_result[0]
    else:
        exit("ERROR: No calibrators found within 10 degrees of your target in JSDC. Consider modifying your "
             "constraints!")
    # Cross-check with Gaia DR3 catalogue for IPDfmp (<2), RUWE (<1.4), RVamp, and Vbroad<100 binarity and rapid
    # rotation checks
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "IPDfmp", "RUWE", "RVamp", "Vbroad", "+_r"],
                    catalog="I/355/gaiadr3")
    print(f"-->Querying {BLUE}Gaia DR3 Catalogue{RESET}...")
    # The default vizier request times out at 60 seconds. You usually will not hit this limit at a row limit of 100, but
    # if you increased the row limit above, uncomment the following line and increase the timeout to a larger number.
    # WARNING: Doing this, depending on how much you increased the row limit by, can crash the program!!
    # vizier.TIMEOUT = 120
    gaia_result = vizier.query_region(jmmc_result, radius="10s", column_filters={"IPDfmp": "<2", "RUWE": "<1.4",
                                                                           "Vbroad": "<100"})
    print(f"-->{GREEN}Query complete!{RESET}")
    if len(gaia_result) > 0:
        gaia_result = gaia_result[0]
    else:
        exit("ERROR: No calibrators found within 10 degrees of your target in Gaia DR3. Consider modifying your "
             "constraints!")
    if gaia_comp_check:
        vizier_neighbors = Vizier(columns=["_RAJ2000", "_DEJ2000", "IPDfmp", "RUWE", "RVamp", "Vbroad", "+_r"],
                                  catalog="I/355/gaiadr3")

        # Now can print out each entry and catch Gaia DR3 companions
        vizier_neighbors.ROW_LIMIT = -1

        print(f"-->Checking for close Gaia companions within {gaia_comp_check}\"")
        neighbors = vizier_neighbors.query_region(gaia_result, radius=f"{gaia_comp_check}s")[0]

        removal_list = ([item for item, count in collections.Counter(neighbors['_q']).items() if count > 1])

        if len(removal_list) > 0:
            gaia_result = gaia_result[~np.isin(gaia_result['_q'], removal_list)]

    ind = gaia_result['_q'] - 1
    jmmc_cols = Table([jmmc_result['Name'][ind], jmmc_result['_r'][ind], jmmc_result['_RAJ2000'][ind],
                       jmmc_result['_DEJ2000'][ind], jmmc_result['SpType'][ind], jmmc_result['Vmag'][ind],
                       jmmc_result['Rmag'][ind], jmmc_result['Hmag'][ind], jmmc_result['Kmag'][ind],
                       jmmc_result['UDDH'][ind], jmmc_result['UDDK'][ind]])

    gaia_cols = Table([gaia_result['IPDfmp'], gaia_result['RUWE'], gaia_result['RVamp'], gaia_result['Vbroad']])

    first_cross_check_table = hstack([jmmc_cols, gaia_cols])

    # Cross-check with Kervella catalogue for binarity (should all be 0)
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "Name", "DMS", "W", "BinH", "BinG2"], catalog="J/A+A/623/A72")
    print(f"-->Querying {BLUE}Kervella et al. 2019 Catalogue{RESET}...")
    kervella_result = vizier.query_region(gaia_result, radius="10s", column_filters={"DMS": "=0", "W": "=0",
                                                                           "BinH": "=0", "BinG2": "=0"})
    print(f"-->{GREEN}Query complete!{RESET}")
    if len(kervella_result) > 0:
        kervella_result = kervella_result[0]
    else:
        exit("ERROR: No calibrators found within 10 degrees of your target in the Kervella et al. 2019 Catalogue. "
             "Consider modifying your constraints!")

    kervella_cols = Table([kervella_result['DMS'], kervella_result['W'], kervella_result['BinH'],
                           kervella_result['BinG2']])

    ind = kervella_result['_q'] - 1
    second_cross_check_table = hstack([first_cross_check_table[ind], kervella_cols])

    # Cross-check with Cruzalebes catalogue for possible use as calibrators (CalFlag and IRflag should be 0)
    vizier = Vizier(columns=["Diam-GAIA", "CalFlag", "IRflag"], catalog="II/361/mdfc-v10")
    print(f"-->Querying {BLUE}Cruzalebes et al. 2019 Catalogue (MDFC){RESET}...")
    cruzalebes_result = vizier.query_region(kervella_result, radius="10s", column_filters={"CalFlag": "=0",
                                                                                           "IRflag": "=0"})
    print(f"-->{GREEN}Query complete!{RESET}")
    if len(cruzalebes_result) > 0:
        cruzalebes_result = cruzalebes_result[0]
    else:
        exit("ERROR: No calibrators found within 10 degrees of your target in MDFC. Consider modifying your "
             "constraints!")
    cruzalebes_cols = Table([cruzalebes_result['Diam-GAIA'], cruzalebes_result['CalFlag'], cruzalebes_result['IRflag']])

    ind = cruzalebes_result['_q'] - 1

    # fcct = final cross-check table
    fcct = hstack([second_cross_check_table[ind], cruzalebes_cols])

    # Check that Gaia estimated angular diameter is within 0.075 mas of JMMC's reported UDDH and UDDK
    fcct = fcct[(abs(fcct['Diam-GAIA']-fcct['UDDH']) < 0.075) & (abs(fcct['Diam-GAIA']-fcct['UDDK']) < 0.075)]

    # By default, CATCH will save the calibrator coordinates in sexagesimal. If you prefer them in degrees, comment out
    # the next 3 lines
    coords = SkyCoord(ra=fcct['_RAJ2000'].data * u.deg, dec=fcct['_DEJ2000'].data * u.deg, frame='icrs')
    fcct['_RAJ2000'] = coords.ra.to_string(unit=u.hour, sep=' ', pad=True, precision=2)
    fcct['_DEJ2000'] = coords.dec.to_string(unit=u.deg, sep=' ', alwayssign=True, pad=True, precision=2)

    # Add one or more comment lines
    fcct.meta['comments'] = [f'Calibrators for {star_name} (RA: {star_ra}, DEC: {star_dec}, '
                             f'V Mag: {star_v_mag:.2f})']
    fcct.write(f'{star_name}_Calibrators.txt', format='ascii.fixed_width', delimiter="", overwrite=True)

    t2 = time.perf_counter()
    if len(fcct['Name']) > 0:
        print(f"Found {YELLOW}{len(fcct['Name'])}{RESET} viable calibrators in {round(t2 - t1, 2)} seconds!")
    else:
        print(f"{RED}Found no viable calibrators!{RESET}")

    return


def cal_checker(calibrator_name: str, gaia_comp_check: bool = False) -> None:
    """
    Checks chosen calibrator stars using magnitude and diameter checks from the JMMC Stellar Diameters Catalogue
    (https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=II/346/jsdc_v2) and binarity checks from the Gaia DR3
    (https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3), Kervella et al. 2019
    (https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/623/A72), and Cruzalebes et al. 2019
    (https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/361) catalogues.

    :param calibrator_name: The name of the chosen calibrator star
    :param gaia_comp_check: Flag that checks for nearby Gaia DR3 companions. Default = False. Change to True if you want
    to check calibrator viability based on Gaia DR3 companions within 5".
    check.

    :return: None. Prints whether the chosen calibrator is viable or not.
    """
    t1 = time.perf_counter()

    star_details_table = Simbad.query_object(f"{calibrator_name}")

    if len(star_details_table) == 0:
        exit(f"Error: {YELLOW}{calibrator_name}{RESET} not found in SIMBAD. Please check that you typed it in correctly!")
    else:
        star_dec = star_details_table['dec'].value[0]

        if star_dec <= -25:
            exit("This star is outside the declination limits (DEC > -25) for the CHARA Array!")
        else:
            pass

    if gaia_comp_check:
        init_check_pass_count = 5
    else:
        init_check_pass_count = 4

    check_pass_count = init_check_pass_count

    print(f"Checking calibrator viability of: {YELLOW}{calibrator_name}{RESET}")
    # Check with JMMC Stellar Diameters Catalogue (Vmag < 9.0, Hmag < 6.4, UDDH < 0.4, SpType = GKM)
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "Name", "SpType", "Vmag", "Rmag", "Hmag", "Kmag", "UDDH", "UDDK",
                             "+_r"], catalog="II/346/jsdc_v2")

    vizier.ROW_LIMIT = 100
    print(f"-->Querying {BLUE}JMMC Stellar Diameters Catalogue (JSDC){RESET}...")
    jmmc_result = vizier.query_region(f"{calibrator_name}", radius="5s")

    if len(jmmc_result) != 0:
        jmmc_result = jmmc_result[0]
        print(f"-->{GREEN}Query complete!{RESET}")
        print(f"---->Query results: V mag = {jmmc_result['Vmag'][0]:.2f}, H mag = {jmmc_result['Hmag'][0]:.2f}, "
              f"UDDH = {jmmc_result['UDDH'][0]:.3f}")
        if (jmmc_result['Vmag'] > 9) |  (jmmc_result['Hmag'] > 6.4) | (jmmc_result['UDDH'] > 0.5):
            print(f"-->{RED}{calibrator_name} fails JMMC Stellar Diameters Catalogue (JSDC) checks!{RESET}")
            check_pass_count -= 1
            if jmmc_result['Vmag'] > 9:
                print(f"---->{RED}Calibrator Vmag > 9!{RESET}")
            if jmmc_result['Hmag'] > 6.4:
                print(f"---->{RED}Calibrator Hmag > 6.4!{RESET}")
            if jmmc_result['UDDH'] > 0.5:
                print(f"---->{RED}Calibrator UDDH > 0.5!{RESET}")
        else:
            print(f"-->{GREEN}{calibrator_name} passes JMMC Stellar Diameters Catalogue (JSDC) checks!{RESET}")
            pass

    else:
        print(f"-->{ORANGE}{calibrator_name} not found in JMMC Stellar Diameters Catalogue (JSDC) {RESET} — "
              f"Check against other catalogues!{RESET}")
        check_pass_count -= 1

    # Cross-check with Gaia DR3 catalogue for IPDfmp (<2), RUWE (<1.4), RVamp, and Vbroad<100 binarity and rapid
    # rotation checks
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "IPDfmp", "RUWE", "RVamp", "Vbroad", "+_r"],
                    catalog="I/355/gaiadr3")
    print(f"-->Querying {BLUE}Gaia DR3 Catalogue{RESET}...")
    gaia_result = vizier.query_region(f"{calibrator_name}", radius="10s")

    if len(gaia_result) == 0:
        exit(f"ERROR: {calibrator_name} not found in Gaia DR3!")
    else:
        gaia_result = gaia_result[0]
        print(f"-->{GREEN}Query complete!{RESET}")
        if len(gaia_result) > 1:
            print(f"-->{RED}Warning: Potential calibrator has Gaia DR3 companions within 10\"{RESET}")
            if gaia_comp_check:
                check_pass_count -= 1
                print(f"---->Companions shown below. _r corresponds to distance from calibrator in arcseconds")
                print(gaia_result[1:])

        if gaia_result[0]['Vbroad'] is np.ma.masked:
            print(f"-->{RED}Warning: Calibrator does not have a listed Vbroad!{RESET}")
            check_pass_count -= 5
            print(f"---->Query results: IPDfmp = {gaia_result[0]['IPDfmp']:.1f}, RUWE = {gaia_result[0]['RUWE']:.2f}")

        else:
            print(f"---->Query results: IPDfmp = {gaia_result[0]['IPDfmp']:.1f}, "
                  f"RUWE = {gaia_result[0]['RUWE']:.2f}, Vbroad = {gaia_result[0]['Vbroad']:.2f} km/s")

            if (((gaia_result[0]['IPDfmp'] is not np.ma.masked) and (gaia_result[0]['IPDfmp'] > 2)) |
                    ((gaia_result[0]['RUWE'] is not np.ma.masked) and (gaia_result[0]['RUWE'] > 1.4)) |
                    ((gaia_result[0]['Vbroad'] is not np.ma.masked) and (gaia_result[0]['Vbroad'] > 100))):

                print(f"-->{RED}{calibrator_name} fails Gaia DR3 Catalogue checks!{RESET}")
                check_pass_count -= 1
                if gaia_result[0]['IPDfmp'] > 2:
                    print(f"---->{RED}Calibrator IPDfmp > 2!{RESET}")
                if gaia_result[0]['RUWE'] > 1.4:
                    print(f"---->{RED}Calibrator RUWE > 1.4!{RESET}")
                if gaia_result[0]['Vbroad'] > 100:
                    check_pass_count -= 5
                    print(f"---->{RED}Calibrator Vbroad > 100!{RESET}")

            else:
                print(f"-->{GREEN}{calibrator_name} passes Gaia DR3 Catalogue checks!{RESET}")
                pass

    # Cross-check with Kervella catalogue for binarity (should all be 0)
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "Name", "DMS", "W", "BinH", "BinG2"], catalog="J/A+A/623/A72")
    print(f"-->Querying {BLUE}Kervella et al. 2019 Catalogue{RESET}...")
    kervella_result = vizier.query_region(f"{calibrator_name}", radius="5s")
    if len(kervella_result) > 0:
        print(f"-->{GREEN}Query complete!{RESET}")
        print(f"---->Query results: DMS = {kervella_result[0]['DMS'][0]}, W = {kervella_result[0]['W'][0]}, "
              f"BinH = {kervella_result[0]['BinH'][0]}, BinG2 = {kervella_result[0]['BinG2'][0]}")
        if ((kervella_result[0]['DMS'] != 0) | (kervella_result[0]['W'] != 0) | (kervella_result[0]['BinH'] != 0)
                | kervella_result[0]['BinG2'] != 0):
            print(f"-->{RED}{calibrator_name} fails Kervella et al. 2019 Catalogue checks!{RESET}")
            check_pass_count -= 5
            if kervella_result[0]['DMS'] != 0:
                print(f"---->{RED}Calibrator DMS ≠ 0!{RESET}")
            if kervella_result[0]['W'] != 0:
                print(f"---->{RED}Calibrator W ≠ 0!{RESET}")
            if kervella_result[0]['BinH'] != 0:
                print(f"---->{RED}Calibrator BinH ≠ 0!{RESET}")
            if kervella_result[0]['BinG2'] != 0:
                print(f"---->{RED}Calibrator BinG2 ≠ 0!{RESET}")
        else:
            print(f"-->{GREEN}{calibrator_name} passes Kervella et al. 2019 Catalogue checks!{RESET}")
            pass

    else:
        print(f"-->{ORANGE}Warning: {calibrator_name} not found in Kervella et al. 2019 Catalogue "
              f"— check against other catalogues!{RESET}")
        check_pass_count -= 1


    # Cross-check with Cruzalebes catalogue for possible use as calibrators (CalFlag and IRflag should be 0)
    vizier = Vizier(columns=["Diam-GAIA", "CalFlag", "IRflag"], catalog="II/361/mdfc-v10")
    print(f"-->Querying {BLUE}Cruzalebes et al. 2019 Catalogue (MDFC){RESET}...")
    cruzalebes_result = vizier.query_region(f"{calibrator_name}", radius="5s")

    if len(cruzalebes_result) != 0:
        cruzalebes_result = cruzalebes_result[0]
        print(f"-->{GREEN}Query complete!{RESET}")
        print(f"---->Query results: CalFlag = {cruzalebes_result['CalFlag'][0]}, IRflag = {cruzalebes_result['IRflag'][0]}")
        if (cruzalebes_result['CalFlag'] != 0) | (cruzalebes_result['IRflag'] == 7):
            print(f"-->{RED}{calibrator_name} fails Cruzalebes et al. 2019 Catalogue (MDFC) checks!{RESET}")
            check_pass_count -= 1
            if cruzalebes_result['CalFlag'] != 0:
                print(f"---->{RED}Calibrator CalFlag ≠ 0!{RESET}")
            if cruzalebes_result['IRflag'] == 7:
                print(f"---->{RED}Calibrator IRFlag ≠ 7!{RESET}")
        else:
            print(f"-->{GREEN}{calibrator_name} passes Cruzalebes et al. 2019 Catalogue (MDFC) checks!{RESET}")
            pass

    else:
        print(f"-->{ORANGE}{calibrator_name} not found in Cruzalebes et al. 2019 Catalogue (MDFC) — "
              f"check against other catalogues!{RESET}")
        check_pass_count -= 1

    t2 = time.perf_counter()
    if check_pass_count / init_check_pass_count == 1:
        print(f"-->{YELLOW}{calibrator_name}{RESET} passed {GREEN}{check_pass_count}/{init_check_pass_count}{RESET} checks")
        print(f"Confirmed {YELLOW}{calibrator_name}{RESET} is likely an {GREEN}ideal{RESET} calibrator in {round(t2 - t1, 2)} seconds!")
    elif (check_pass_count / init_check_pass_count < 1) & (check_pass_count / init_check_pass_count >= 0.7):
        print(f"-->{YELLOW}{calibrator_name}{RESET} passed {ORANGE}{check_pass_count}/{init_check_pass_count}{RESET} checks")
        print(f"Confirmed {YELLOW}{calibrator_name}{RESET} is likely a {ORANGE}usable{RESET} calibrator in {round(t2 - t1, 2)} seconds!")
    else:
        print(f"{YELLOW}{calibrator_name}{RESET} {RED}is unlikely to be a viable calibrator!{RESET}")
        print("We recommend submitting this star to the JMMC Bad Calibrators Database: https://www.jmmc.fr/badcal/\n")

    return


def main():
    # If you get a timeout error, un-comment the following line to use the Japanese mirror site for Vizier
    # conf.server = "vizier.nao.ac.jp"
    print("Current Vizier server:", conf.server)
    main_question = (
        input(f"Would you like to find calibrators for a science target {BLUE}(type A){RESET}, "
              f"or check a possible calibrator's viability {BLUE}(type B){RESET}?:\n"))
    if main_question in ["A", "a"]:
        target_star_name = input("Please enter the name of your target (please ensure the name is resolvable in SIMBAD):\n")
        gaia_question = input("Would you like to filter calibrators by whether they have close companions in Gaia DR3 "
                              "Y/[N]?\n").strip()
        if gaia_question in ["Y", "y"]:
            while True:
                gaia_radius = input("Please enter the desired cutoff radius (in arcseconds) for Gaia companions:\n").strip()
                if not gaia_radius:
                    print(f"{RED}Invalid cutoff radius, please enter a number!{RESET}")
                    continue
                try:
                    gaia_radius = float(gaia_radius)
                    cal_finder(target_star_name, gaia_radius)
                    break
                except ValueError:
                    print(f"{RED}Invalid cutoff radius, please enter a number!{RESET}")

        else:
            cal_finder(target_star_name)


    elif main_question in ["B", "b"]:
        multiple_cal_check = input(f"Would you like to check a single calibrator, or multiple calibrators? [S]/M\n")
        if multiple_cal_check in ["S", "s", ""]:
            target_star_name = input("Please enter the name of your calibrator (please ensure the name is resolvable in SIMBAD):\n")
            gaia_question = input("Would you like to filter calibrators by whether it has a companion within 10\" in Gaia DR3 "
                                  "Y/[N]?\n").strip()
            if gaia_question in ["Y", "y"]:
                cal_checker(target_star_name, gaia_comp_check=True)

            else:
                cal_checker(target_star_name)

        elif multiple_cal_check in ["M", "m"]:
            gaia_question = input(
                "Would you like to filter calibrators by whether it has a companion within 10\" in Gaia DR3 "
                "Y/[N]?\n").strip()
            target_star_name_list = input(
                "Please enter the names of your calibrators as a comma-separated list (please ensure the names are "
                "resolvable in SIMBAD):\n")

            target_star_name_list = [target_star_name.strip() for target_star_name in target_star_name_list.split(',')]

            for target_star_name in target_star_name_list:
                if gaia_question in ["Y", "y"]:
                    cal_checker(target_star_name, gaia_comp_check=True)

                else:
                    cal_checker(target_star_name)

    else:
        exit(f"{YELLOW}No option selected. Have a good day!{RESET}")


if __name__ == '__main__':
    print(
        f"""
        ############################ This is {RED}C{RESET}{GREEN}A{RESET}{YELLOW}T{RESET}{BLUE}C{RESET}{MAGENTA}H{RESET} ############################
                   [{RED}C{RESET}]HARA [{GREEN}A{RESET}]rray's [{YELLOW}T{RESET}]hrifty [{BLUE}C{RESET}]alibrator [{MAGENTA}H{RESET}]unter
                                 Version: {__version__}                           
                       https://github.com/akshat-chaturvedi/CATCH                 
        #######################################################################
        """
        )
    main()