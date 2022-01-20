nextflow.enable.dsl=2

process PULL_SINGULARITY {

    input:
    path config_file
    path script_py

    output:
    path 'finished.txt'

    script:
    """
    workd=\$(pwd)
    python ${script_py} ${config_file} > containers.txt

    input="containers.txt"
    while IFS= read -r line
    do
      #echo Pulling singularity image: \$line
      singularity pull \$line
    done < "\$input"
    
    mv *.sif ${params.singularity_store}

    # pseudo process to make other processes wait	
    echo finished > finished.txt
    """

}


workflow PULL_CONTAINERS {

    take:
      config_file_ch
      pull_file_ch

    main:
      PULL_SINGULARITY(  config_file_ch,
			 pull_file_ch  )

    emit:
      PULL_SINGULARITY.out
}
