#!/bin/bash

# first of all, set some env variables
# the environment variables will be needed to access the private docker registry
cat ./data/bash_scripts/echo_env.sh
source ./data/bash_scripts/echo_env.sh

# convert proteomics raw to mzML files
nextflow run \
	main.nf \
	-profile slurm \
	-entry RAW_TO_MZML

# this process pulls singularity packages
# to the place in your system where they belong
# nextflow implementation currently not working
# since containerization with singularity running
# on singularity is difficult
: '
nextflow run \
    main.nf \
    -entry PULLING \
	-profile slurm_offline
'
# for now, instead use this direct python call by running it in the ${projectDir} directory
# ATTENTION: needs >=python3.7 to run
python ~/riboseq_pipeline/data/python_scripts/pull_containers.py \
	--config ~/riboseq_pipeline/conf/slurm.config \
	--out ~/riboseq_pipeline/conf/slurm_offline.config \
	--dest ~/.singularity/cache/library/

# we might also want to pull the newest fasta genome and annotation files:
# This subworkflow takes a special profile, since it has to run on the login node
# which is the only one connected to the internet
nextflow run \
	main.nf \
	-entry PULL_FILES \
	-profile pull

# the following nf process generates test files
nextflow run \
	main.nf \
	-entry TEST_FILES

# if your proteomics filenames do not clearly indicate the experimental run
# and are not in the right format, this process will fix it.
# In order for this to work, the file experimental_conditions.csv in the data folder
# has to be quickly filled out.
nextflow run \
	main.nf \
	-entry MAP_NAMES \
	-profile slurm_offline

# finally call sbatch to make the main workflow run
sbatch slurm_script
# or just: nextflow run main.nf -profile slurm_offline
