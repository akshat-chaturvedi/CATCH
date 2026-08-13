#!/usr/bin/env python

"""catch.py: A program to cross-match across multiple Vizier databases to return suitable calibrator stars for
observations with the CHARA Array"""

__author__ = "Akshat S. Chaturvedi and Mahir M. Patel"
__credits__ = ["Akshat S. Chaturvedi", "Mahir M. Patel", "Colin Kane", "Becky Flores", "Jeremy Jones"]
__license__ = "MIT"
__version__ = "2.0 | 2026/08/12" #Added Instruments for MRIC-X, MYSTIC, SILMARIL
__maintainer__ = "Akshat S. Chaturvedi"
__email__ = "achaturvedi3@gsu.edu"
__status__ = "Production"

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
from astropy.utils.metadata import MergeConflictWarning
import requests
from argparse import ArgumentParser
import sys
from catch_mircx_mystic import *
from catch_silmaril import *
from catch_spica import *

Vizier.clear_cache()

warnings.simplefilter("ignore", NoResultsWarning)
warnings.simplefilter("ignore", MergeConflictWarning)

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
# __version__ = '1.1 | 2025/08/27' # Changes to file output and printout formats
# __version__ = '1.2 | 2025/10/31' # Added e_LDD to print output, vizier server error messages
# __version__ = '1.3 | 2026/02/11' #  Added server switch capability if normal Vizier server is down, updated README

def questions(instrument: str):
    if len(sys.argv) <= 1:
        print(f"Vizier server: {GREEN}{conf.server}{RESET}")
        vizier_web_status = requests.get(
            "https://" + f"{conf.server}"
        ).status_code

        if vizier_web_status != 200:
            print(
                f"Vizier is currently down "
                f"(HTML Response Code: {RED}{vizier_web_status}{RESET}), "
                f"switching to {GREEN}vizier.nao.ac.jp{RESET}"
            )
            Vizier.VIZIER_SERVER = "vizier.nao.ac.jp"

        else:
            print(
                f"{GREEN}{conf.server}{RESET} server up "
                f"(HTML Response Code: {GREEN}{vizier_web_status}{RESET})"
            )

        main_question = input(
            f"Would you like to find calibrators for a science target "
            f"{BLUE}(type A){RESET}, or check a possible calibrator's "
            f"viability {BLUE}(type B){RESET}?:\n"
        )

        if main_question in ["A", "a"]:
            target_star_name = input(
                "Please enter the name of your target "
                "(please ensure the name is resolvable in SIMBAD):\n"
            )

            gaia_question = input(
                "Would you like to filter calibrators by whether they "
                "have close companions in Gaia DR3 Y/[N]?\n"
            ).strip()

            if gaia_question in ["Y", "y"]:
                while True:
                    gaia_radius = input(
                        "Please enter the desired cutoff radius "
                        "(in arcseconds) for Gaia companions:\n"
                    ).strip()

                    if not gaia_radius:
                        print(
                            f"{RED}Invalid cutoff radius, "
                            f"please enter a number!{RESET}"
                        )
                        continue

                    try:
                        gaia_radius = float(gaia_radius)
                        if instrument == "HK" or instrument == "" :
                            hk_cal_finder(target_star_name, gaia_radius)
                        elif instrument == "S":
                            s_cal_finder(target_star_name, gaia_radius)
                        elif instrument == "R":
                            r_cal_finder(target_star_name, gaia_radius)

                        break

                    except ValueError:
                        print("b")
                        print(f"{RED}Invalid cutoff radius, please enter a number!{RESET}")
                        break

            else:
                if instrument == "HK" or instrument == "" :
                    hk_cal_finder(target_star_name)
                elif instrument == "S":
                    s_cal_finder(target_star_name)
                elif instrument == "R":
                    r_cal_finder(target_star_name)

        elif main_question in ["B", "b"]:
            multiple_cal_check = input(
                "Would you like to check a single calibrator, "
                "or multiple calibrators? [S]/M\n"
            )

            if multiple_cal_check in ["S", "s", ""]:
                target_star_name = input(
                    "Please enter the name of your calibrator "
                    "(please ensure the name is resolvable in SIMBAD):\n"
                )

                gaia_question = input(
                    "Would you like to filter calibrators by whether "
                    "it has a companion within 10\" in Gaia DR3 Y/[N]?\n"
                ).strip()

                if gaia_question in ["Y", "y"]:
                    if instrument == "HK" or instrument == "" :
                        hk_cal_checker(target_star_name, gaia_comp_check=True)
                    elif instrument == "S":
                        s_cal_checker(target_star_name, gaia_comp_check=True)
                    elif instrument == "R":
                        r_cal_checker(target_star_name, gaia_comp_check=True)

                else:
                    if instrument == "HK" or instrument == "" :
                        hk_cal_checker(target_star_name)
                    elif instrument == "S":
                        s_cal_checker(target_star_name)
                    elif instrument == "R":
                        r_cal_checker(target_star_name)

            elif multiple_cal_check in ["M", "m"]:
                gaia_question = input(
                    "Would you like to filter calibrators by whether "
                    "it has a companion within 10\" in Gaia DR3 Y/[N]?\n"
                ).strip()

                target_star_name_list = input(
                    "Please enter the names of your calibrators as a "
                    "comma-separated list (please ensure the names are "
                    "resolvable in SIMBAD):\n"
                )

                target_star_name_list = [
                    target_star_name.strip()
                    for target_star_name
                    in target_star_name_list.split(",")
                ]

                for target_star_name in target_star_name_list:
                    if gaia_question in ["Y", "y"]:
                        if instrument == "HK" or instrument == "" :
                            hk_cal_checker(target_star_name, gaia_comp_check=True)
                        elif instrument == "S":
                            s_cal_checker(target_star_name, gaia_comp_check=True)
                        elif instrument == "R":
                            r_cal_checker(target_star_name, gaia_comp_check=True)

                    else:
                        if instrument == "HK" or instrument == "" :
                            hk_cal_checker(target_star_name)
                        elif instrument == "S":
                            s_cal_checker(target_star_name)
                        elif instrument == "R":
                            r_cal_checker(target_star_name)

        else:
            exit(
                f"{YELLOW}No option selected. "
                f"Have a good day!{RESET}"
            )

    else:
        parser = ArgumentParser(epilog="Happy fringing!")

        main_choice_group = parser.add_mutually_exclusive_group()

        main_choice_group.add_argument(
            "-i",
            "--instrument",
            help=(
                "Choose an instrument "
                f"{BLUE}HK{RESET} for MIRC-X and MYSTIC; type {YELLOW}S{RESET} for Silmaril; type {RED}R{RESET} for Spica"
            ),
            action="store_true"
        )

        main_choice_group.add_argument(
            "-f",
            "--find",
            help=(
                "Choose this option if you would like to find "
                "calibrators for a science target"
            ),
            action="store_true"
        )

        parser.add_argument(
            "-t",
            "--target",
            help="Add the science target to find calibrators for",
            dest="science_target"
        )

        parser.add_argument(
            "-g",
            "--gaia",
            help=(
                "Filter calibrators by whether they have companions "
                "within 5 arcseconds in Gaia DR3"
            ),
            action="store_true"
        )

        main_choice_group.add_argument(
            "-v",
            "--verify",
            help=(
                "Choose this option if you would like to check a "
                "possible calibrator's viability"
            ),
            action="store_true"
        )

        parser.add_argument(
            "-s",
            "--single",
            help="Verify a single calibrator",
            action="store_true"
        )

        parser.add_argument(
            "-m",
            "--multiple",
            help="Verify multiple calibrators",
            action="store_true"
        )

        parser.add_argument(
            "-c",
            "--calibrator",
            help=(
                "The name of the possible calibrator you want "
                "to verify, in quotes"
            ),
            dest="cand_calibrator"
        )

        parser.add_argument(
            "-cl",
            "--cal_list",
            help=(
                "Add comma separated list of possible calibrators "
                "in quotes"
            ),
            dest="list_of_cals"
        )

        args = parser.parse_args()

        if args.find:
            print(f"Vizier server: {GREEN}{conf.server}{RESET}")

            vizier_web_status = requests.get(
                "https://" + f"{conf.server}"
            ).status_code

            if vizier_web_status != 200:
                print(
                    f"Vizier is currently down "
                    f"(HTML Response Code: "
                    f"{RED}{vizier_web_status}{RESET}), "
                    f"switching to {GREEN}vizier.nao.ac.jp{RESET}"
                )

                Vizier.VIZIER_SERVER = "vizier.nao.ac.jp"

            else:
                print(
                    f"{GREEN}{conf.server}{RESET} server up "
                    f"(HTML Response Code: "
                    f"{GREEN}{vizier_web_status}{RESET})"
                )

            if args.gaia:
                if args.instrument == "HK":
                    hk_cal_checker(args.science_target, gaia_comp_check=5)
                elif args.instrument == "S":
                    s_cal_checker(args.science_target, gaia_comp_check=5)
                elif args.instrument == "R":
                    r_cal_checker(args.science_target, gaia_comp_check=5)

            else:
                if args.instrument == "HK":
                    hk_cal_checker(args.science_target)
                elif args.instrument == "S":
                    s_cal_checker(args.science_target)
                elif args.instrument == "R":
                    r_cal_checker(args.science_target)

        elif args.verify:
            print(f"Vizier server: {GREEN}{conf.server}{RESET}")

            vizier_web_status = requests.get(
                "https://" + f"{conf.server}"
            ).status_code

            if vizier_web_status != 200:
                print(
                    f"Vizier is currently down "
                    f"(HTML Response Code: "
                    f"{RED}{vizier_web_status}{RESET}), "
                    f"switching to {GREEN}vizier.nao.ac.jp{RESET}"
                )

                Vizier.VIZIER_SERVER = "vizier.nao.ac.jp"

            else:
                print(
                    f"{GREEN}{conf.server}{RESET} server up "
                    f"(HTML Response Code: "
                    f"{GREEN}{vizier_web_status}{RESET})"
                )

            if args.single:
                if args.gaia:
                    if args.instrument == "HK":
                        hk_cal_checker(args.cand_calibrator, gaia_comp_check=True)
                    elif args.instrument == "S":
                        s_cal_checker(args.cand_calibrator, gaia_comp_check=True)
                    elif args.instrument == "R":
                        r_cal_checker(args.cand_calibrator, gaia_comp_check=True)

                else:
                    if args.instrument == "HK":
                        hk_cal_checker(args.cand_calibrator)
                    elif args.instrument == "S":
                        s_cal_checker(args.cand_calibrator)
                    elif args.instrument == "R":
                        r_cal_checker(args.cand_calibrator)

            elif args.multiple:
                target_star_name_list = [
                    target_star_name.strip()
                    for target_star_name
                    in args.list_of_cals.split(",")
                ]

                for target_star_name in target_star_name_list:
                    if args.gaia:
                        if args.instrument == "HK":
                            hk_cal_checker(target_star_name, gaia_comp_check=True)
                        elif args.instrument == "S":
                            s_cal_checker(target_star_name, gaia_comp_check=True)
                        elif args.instrument == "R":
                            r_cal_checker(target_star_name, gaia_comp_check=True)

                    else:
                        if args.instrument == "HK":
                            hk_cal_checker(target_star_name)
                        elif args.instrument == "S":
                            s_cal_checker(target_star_name)
                        elif args.instrument == "R":
                            r_cal_checker(target_star_name)



def main():
    instrument_question = input(
        f"What instrument is this for? \n"
        f"Type [{BLUE}HK{RESET}] for MIRC-X and MYSTIC; type {YELLOW}S{RESET} for Silmaril; type {RED}R{RESET} for Spica\n"
    ).strip().upper()
    questions(instrument_question)

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

