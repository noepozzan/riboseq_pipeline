nextflow.enable.dsl=2

process CHECK_FILE_PRESENCE {

    publishDir "./results/check_files", mode: 'copy'

    input:
    path riboseq_reads
    path proteomics_reads
    path genome
    path genome_fai
    path gtf
    path other_RNAs_sequence
    path script_py

    output:
    path 'check_files.log'

    script:
    """
    python ${script_py} 6 \
	${riboseq_reads} \
	${proteomics_reads} \
	${genome} \
	${genome_fai} \
	${gtf} \
	${other_RNAs_sequence} \
	&> check_files.log

    """

}


workflow CHECK_FILES_PIPE {

    take:
      riboseq_reads_ch
      proteomics_reads_ch
      genome_ch
      genome_fai_ch
      gtf_ch
      other_RNAs_sequence_ch
      check_files_script_ch

    main:
      CHECK_FILE_PRESENCE(  riboseq_reads_ch,
			    proteomics_reads_ch,
			    genome_ch,
			    genome_fai_ch,
			    gtf_ch,
			    other_RNAs_sequence_ch,
			    check_files_script_ch  )

}
