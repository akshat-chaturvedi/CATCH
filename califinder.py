from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astropy.table import Table, hstack
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
RESET = '\033[0m'
ITALIC = '\033[3m'
BLINK = '\033[5m'

__version__ = '0.1 | 2025/08/25' # Created main function to run everything


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

    print(f"Beginning calibration search for target: {YELLOW}{star_name}{RESET}")
    # Check with JMMC Stellar Diameters Catalogue (Vmag < 9.0, 3.9 < Hmag < 6.4, UDDH < 0.4, SpType = GKM)
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "Name", "SpType", "Vmag", "Rmag","Hmag", "Kmag", "UDDH", "UDDK",
                             "+_r"], catalog="II/346/jsdc_v2")

    vizier.ROW_LIMIT = 100
    print(f"-->Querying {BLUE}JMMC Stellar Diameters Catalogue{RESET}...")
    jmmc_result = vizier.query_region(f"{star_name}", radius="10d", column_filters={"Vmag":"<9.0", "Hmag":"3.9 .. 6.4",
                                                                           "UDDH": "<0.4", "_DEJ2000": ">-25"})[0]
    print(f"-->{GREEN}Query complete!{RESET}")
    # Cross-check with Gaia DR3 catalogue for IPDfmp (<2), RUWE (<1.4), RVamp, and Vbroad<100 binarity and rapid
    # rotation checks
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "IPDfmp", "RUWE", "RVamp", "Vbroad", "+_r"],
                    catalog="I/355/gaiadr3")
    print(f"-->Querying {BLUE}Gaia DR3 Catalogue{RESET}...")
    gaia_result = vizier.query_region(jmmc_result, radius="10s", column_filters={"IPDfmp": "<2", "RUWE": "<1.4",
                                                                           "Vbroad": "<100"})[0]
    print(f"-->{GREEN}Query complete!{RESET}")
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
                                                                           "BinH": "=0", "BinG2": "=0"})[0]
    print(f"-->{GREEN}Query complete!{RESET}")
    kervella_cols = Table([kervella_result['DMS'], kervella_result['W'], kervella_result['BinH'],
                           kervella_result['BinG2']])

    ind = kervella_result['_q'] - 1
    second_cross_check_table = hstack([first_cross_check_table[ind], kervella_cols])

    # Cross-check with Cruzalebes catalogue for possible use as calibrators (CalFlag and IRflag should be 0)
    vizier = Vizier(columns=["Diam-GAIA", "CalFlag", "IRflag"], catalog="II/361/mdfc-v10")
    print(f"-->Querying {BLUE}Cruzalebes et al. 2019 Catalogue{RESET}...")
    cruzalebes_result = vizier.query_region(kervella_result, radius="10s", column_filters={"CalFlag": "=0",
                                                                                           "IRflag": "=0"})[0]
    print(f"-->{GREEN}Query complete!{RESET}")
    cruzalebes_cols = Table([cruzalebes_result['Diam-GAIA'], cruzalebes_result['CalFlag'], cruzalebes_result['IRflag']])

    ind = cruzalebes_result['_q'] - 1

    # fcct = final cross-check table
    fcct = hstack([second_cross_check_table[ind], cruzalebes_cols])

    # Check that Gaia estimated angular diameter is within 0.075 mas of JMMC's reported UDDH and UDDK
    fcct = fcct[(abs(fcct['Diam-GAIA']-fcct['UDDH']) < 0.075) & (abs(fcct['Diam-GAIA']-fcct['UDDK']) < 0.075)]

    # Add one or more comment lines
    fcct.meta['comments'] = [f'Calibrators for {star_name} (RA: {star_ra:.5f}, DEC: {star_dec:.5f}, '
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
        pass

    if gaia_comp_check:
        init_check_pass_count = 5
    else:
        init_check_pass_count = 4

    check_pass_count = init_check_pass_count

    print(f"Checking calibrator viability of: {YELLOW}{calibrator_name}{RESET}")
    # Check with JMMC Stellar Diameters Catalogue (Vmag < 9.0, 3.9 < Hmag < 6.4, UDDH < 0.4, SpType = GKM)
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "Name", "SpType", "Vmag", "Rmag", "Hmag", "Kmag", "UDDH", "UDDK",
                             "+_r"], catalog="II/346/jsdc_v2")

    vizier.ROW_LIMIT = 100
    print(f"-->Querying {BLUE}JMMC Stellar Diameters Catalogue{RESET}...")
    jmmc_result = vizier.query_region(f"{calibrator_name}", radius="5s")[0]
    print(f"-->{GREEN}Query complete!{RESET}")

    if jmmc_result:
        if ((jmmc_result['Vmag'] > 9) | (jmmc_result['Hmag'] < 3.9) | (jmmc_result['Hmag'] > 6.4) |
                (jmmc_result['UDDH'] > 0.5)):
            print(f"-->{RED}{calibrator_name} fails JMMC Stellar Diameters Catalogue checks!{RESET}")
            check_pass_count -= 1
        else:
            pass

    else:
        print(f"-->{YELLOW}{calibrator_name} not found in JMMC Stellar Diameters Catalogue{RESET}")
        print(f"-->{YELLOW}Check against other catalogues!{RESET}")

    # Cross-check with Gaia DR3 catalogue for IPDfmp (<2), RUWE (<1.4), RVamp, and Vbroad<100 binarity and rapid
    # rotation checks
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "IPDfmp", "RUWE", "RVamp", "Vbroad", "+_r"],
                    catalog="I/355/gaiadr3")
    print(f"-->Querying {BLUE}Gaia DR3 Catalogue{RESET}...")
    gaia_result = vizier.query_region(f"{calibrator_name}", radius="5s")[0]
    print(f"-->{GREEN}Query complete!{RESET}")

    if len(gaia_result) > 1:
        print(f"-->{RED}Warning: Potential calibrator has Gaia DR3 companions within 5\"{RESET}")
        if gaia_comp_check:
            check_pass_count -= 1

    if ((gaia_result[0]['IPDfmp'] > 2) | (gaia_result[0]['RUWE'] > 1.4) |
            ((gaia_result[0]['Vbroad'] is not np.ma.masked) and (gaia_result[0]['Vbroad'] > 100))):
        print(f"-->{RED}{calibrator_name} fails Gaia DR3 Catalogue checks!{RESET}")
        check_pass_count -= 1
    else:
        pass

    # Cross-check with Kervella catalogue for binarity (should all be 0)
    vizier = Vizier(columns=["_RAJ2000", "_DEJ2000", "Name", "DMS", "W", "BinH", "BinG2"], catalog="J/A+A/623/A72")
    print(f"-->Querying {BLUE}Kervella et al. 2019 Catalogue{RESET}...")
    kervella_result = vizier.query_region(f"{calibrator_name}", radius="5s")
    if len(kervella_result) > 0:
        print(f"-->{GREEN}Query complete!{RESET}")
        if ((kervella_result[0]['DMS'] != 0) | (kervella_result[0]['W'] != 0) | (kervella_result[0]['BinH'] != 0)
                | kervella_result[0]['BinG2'] != 0):
            print(f"-->{RED}{calibrator_name} fails Kervella et al. 2019 Catalogue checks!{RESET}")
            check_pass_count -= 1
        else:
            pass

    else:
        print(f"-->{RED}Warning: {calibrator_name} not found in Kervella et al. 2019 Catalogue "
              f"— Check against other catalogues!{RESET}")
        check_pass_count -= 1


    # Cross-check with Cruzalebes catalogue for possible use as calibrators (CalFlag and IRflag should be 0)
    vizier = Vizier(columns=["Diam-GAIA", "CalFlag", "IRflag"], catalog="II/361/mdfc-v10")
    print(f"-->Querying {BLUE}Cruzalebes et al. 2019 Catalogue{RESET}...")
    cruzalebes_result = vizier.query_region(f"{calibrator_name}", radius="5s")[0]
    print(f"-->{GREEN}Query complete!{RESET}")

    if cruzalebes_result:
        if (cruzalebes_result['CalFlag'] != 0) | (cruzalebes_result['IRflag'] == 7):
            print(f"-->{RED}{calibrator_name} fails Cruzalebes et al. 2019 Catalogue checks!{RESET}")
            check_pass_count -= 1
        else:
            pass

    else:
        print(f"-->{YELLOW}{calibrator_name} not found in Cruzalebes et al. 2019 Catalogue — "
              f"Check against other catalogues!{RESET}")

    t2 = time.perf_counter()
    if check_pass_count == init_check_pass_count:
        print(f"-->{YELLOW}{calibrator_name}{RESET} passed {GREEN}{check_pass_count}/{init_check_pass_count}{RESET} checks")
        print(f"Confirmed {YELLOW}{calibrator_name}{RESET} is a viable calibrator in {round(t2 - t1, 2)} seconds!")
    else:
        print(f"-->{YELLOW}{calibrator_name}{RESET} passed {RED}{check_pass_count}/{init_check_pass_count}{RESET} checks")
        print(f"{YELLOW}{calibrator_name}{RESET} {RED}is not a viable calibrator!{RESET}\n")
        print("We recommend submitting this star to the JMMC Bad Calibrators Database: https://www.jmmc.fr/badcal/")

    return


def main():
    main_question = (
        input(f"Would you like to find calibrators for a science target {BLUE}(type A){RESET}, "
              f"or check a possible calibrators viability {BLUE}(type B){RESET}?:\n"))
    if main_question in ["A", "a"]:
        target_star_name = input("Please enter the name of your target (preferably its HD number):\n")
        gaia_question = input("Would you like to filter calibrators by whether they have close companions in Gaia DR3 "
                              "Y/[N]?\n").strip()
        if gaia_question in ["Y", "y"]:
            while True:
                gaia_radius = input("Please enter the desired cutoff radius for Gaia companions:\n").strip()
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
        target_star_name = input("Please enter the name of your calibrator (preferably its HD number):\n")
        gaia_question = input("Would you like to filter calibrators by whether it has a companion within 5\" in Gaia DR3 "
                              "Y/[N]?\n").strip()
        if gaia_question in ["Y", "y"]:
            cal_checker(target_star_name, gaia_comp_check=True)

        else:
            cal_checker(target_star_name)

    else:
        exit(f"{YELLOW}No option selected. Have a good day!{RESET}")


if __name__ == '__main__':
    print(
        f"""
        ############################ This is {RED}Ca{RESET}{GREEN}li{RESET}{YELLOW}Fi{RESET}{BLUE}nd{RESET}{MAGENTA}er{RESET} ############################
                                 Version: {__version__}                           
                       https://github.com/akshat-chaturvedi/CaliFinder                 
        ############################################################################
        """
        )
    main()