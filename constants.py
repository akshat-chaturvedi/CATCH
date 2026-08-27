#!/usr/bin/env python

"""constants.py: Constants and setup parameters for CATCH"""

# Colors and formatting strings for stdout
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
MAGENTA = '\033[95m'
ORANGE = '\033[38;2;255;128;0m'
RESET = '\033[0m'
ITALIC = '\033[3m'
BLINK = '\033[5m'

# Output file format (change to ascii.tab for TSV or ascii.csv for CSV; will need to change delimiter below accordingly)
outfile_format = "ascii.fixed_width"

# Output file delimiter (change to "\t" for TSV or "," for CSV)
outfile_delimiter = ""

