#!/bin/bash

nextflow run \
    main.nf \
    -entry PULLING

sbatch slurm_script
