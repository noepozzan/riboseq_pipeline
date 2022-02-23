nextflow.enable.dsl=2

process PULL_SINGULARITY {

	label "htseq"

    input:
    path slurm_online
    path python_script

    output:
    path 'slurm_offline.config'

    script:
    """
	# extract docker addresses out of conf file
	# and write to output file, also pull the singularity images
	python ${python_script} \
		--config_file ${slurm_online} \
		--out slurm_offline.config \
		--dest ${params.singularity_store}

	cp slurm_offline.config ${projectDir}/conf
    
	"""

}


workflow PULL_CONTAINERS {

    take:
    online_slurm_config
    pull_file_ch

    main:
    PULL_SINGULARITY(
		online_slurm_config,
		pull_file_ch
	)

    emit:
    PULL_SINGULARITY.out
}

