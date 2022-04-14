#!/usr/bin/env python

_version__ = "0.1"
__author__ = "Noe Pozzan"
__contact__ = "noe.pozzan@stud.unibas.ch"
__doc__ = """
            Based on a nextflow config file containing addresses to dockerhub
            registries, these container images are pulled to the user's machine
            and a second config file with the container's path is written.
            ATTENTION: This process needs >=python3.7 and singularity installed
          """

# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import os
from argparse import ArgumentParser, RawTextHelpFormatter
import subprocess
from subprocess import PIPE

def main():
    """ Main function """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--config_file",
        dest="config_file",
        help="""nextflow configuration file containing the addresses
                of the docker images necessary for the pipeline""",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="file to be written with the 'docker pull' commands",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--dest",
        dest="dest",
        help="system path where the singularity images should be placed",
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
                options.config_file,
                os.linesep
            )
        )

    # define a list where all lines will be appended after editing
    lines = []
    # open the config file
    with open(options.config_file, 'r') as f:
        # for each line in the config, check whether it starts with "container"
        for pos_in, line in enumerate(f):
            # line = line.strip()
            if line.strip().startswith('container'):
                # if line starts with "container" extract only interesting part
                line = line.strip()
                new_line = line.split("= ", 1)[1]
                # do some beautyfying
                line_wa = new_line.replace("'", "")
                if options.verbose:
                    print(line_wa)
                # pull the images that we found in the file above
                print(line_wa)
                # p0 = subprocess.run("export SINGULARITY_DOCKER_USERNAME=noepozzan \
                #                     && export SINGULARITY_DOCKER_PASSWORD=testpwd47",
                #                     shell=True,
                #                     capture_output=True,
                #                     text=True,
                #                     check=True)
                #os.environ['SINGULARITY_DOCKER_USERNAME'] = 'noepozzan'
                #os.environ['SINGULARITY_DOCKER_PASSWORD'] = 'testpwd47'
                os.putenv("SINGULARITY_DOCKER_USERNAME", "noepozzan")
                os.putenv("SINGULARITY_DOCKER_PASSWORD", "testpwd47")
                p1 = subprocess.run("singularity pull {}".format(line_wa),
                                    shell=True,
                                    capture_output=True,  # stdout=PIPE, stderr=PIPE,
                                    text=True,
                                    check=True)
                if options.verbose:
                    print("stdout:", p1.stdout)
                    print("stderr:", p1.stderr)

                print("stdout:", p1.stdout)
                print("stderr:", p1.stderr)

                # search for whatever has been pulled
                p2 = subprocess.run("ls -l *.sif",
                                    shell=True,
                                    capture_output=True,  # stdout=PIPE, stderr=PIPE,
                                    text=True,
                                    check=True)
                if options.verbose:
                    print("stdout:", p2.stdout)
                    print("stderr:", p2.stderr)
                # extract the image name out of the "ls -l" string
                if p2.stdout:
                    stripped_stdout = p2.stdout.rstrip()
                    split_stdout = stripped_stdout.split(" ")
                    if options.verbose:
                        print("object to be moved:", split_stdout[-1])
                else:
                    raise ValueError("""The singularity image that was supposed
                                    to be downloaded cannot be found in the
                                    current directory""")

                # check whether the directory to be moved to exists
                # if it does not, create it
                if not os.path.exists(options.dest):
                    subprocess.run("mkdir -p {}"
                                   .format(options.dest),
                                   shell=True,
                                   capture_output=True,  # stdout=PIPE, stderr=PIPE,
                                   text=True,
                                   check=True)

                p3 = subprocess.run("mv {} {}"
                                    .format(split_stdout[-1], options.dest),
                                    shell=True,
                                    capture_output=True,  # stdout=PIPE, stderr=PIPE,
                                    text=True,
                                    check=True)
                # do a bit of string manipulation
                # to have the new string look like the original
                new_line = "        container = "
                tmp_line = os.path.join("file://",
                                        str(options.dest),
                                        str(split_stdout[-1]))
                tmp_line = '"{}"'.format(tmp_line)
                new_line += tmp_line
                new_line += "\n"

                if p3.returncode == 0:
                    lines.append(new_line)

            # uncomment cpus nextflow option for slurm not to fail
            elif line.strip().startswith('cpus'):
                line = line.strip()
                lines.append("        // " + line + "\n")
            # append lines that did not need editing
            else:
                lines.append(line)

    # now write the image's name into the offline slurm config,
    # together with all the unedited lines
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
