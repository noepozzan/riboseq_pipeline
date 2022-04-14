#!/usr/bin/env python

_version__ = "0.1"
__author__ = "Michi Jackson"
__contact__ = "noe.pozzan@stud.unibas.ch"
__doc__ = "Parse a parameter file and change some parameters."

# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import os
from argparse import ArgumentParser, RawTextHelpFormatter


def main():
    """ Main function """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--database",
        dest="database",
        help="database as generated by the philosopher database command",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--param",
        dest="param",
        help="params file as generated by msfragger --config command",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="output file name",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        required=False,
        help="Verbose"
    )

    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if options.verbose:
        sys.stdout.write(
            "Parsing params file: {} {}".format(
                options.param,
                os.linesep
            )
        )

    lines = []
    with open(options.param, "r") as f:
        for position, line in enumerate(f):
            line = line.strip()
            # search for specific lines and change them
            if line.startswith('database_name'):
                line = "database_name = {} \n".format(options.database)
            elif line.startswith('calibrate_mass'):
                line = "calibrate_mass = 0 \n"
            else:
                line = line + "\n"

            lines.append(line)

    with open(options.out, "w") as out_file:
        out_file.writelines(lines)


# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)

