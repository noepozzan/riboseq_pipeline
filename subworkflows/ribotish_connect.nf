nextflow.enable.dsl=2

// https://github.com/zhpn1024/ribotish

process RIBOTISH_QUALITY {

    label "ribotish"

    publishDir "${params.ribotish_connect_outDir}/ribotish_quality", mode: 'copy'

    input:
    each(path(bam_file))
    path gtf_file

    output:
    path '*'

    script:
    """
    ribotish quality \
	-b ${bam_file} \
	-g ${gtf_file} \
	--th ${params.ribotish_quality_th}
    """

}

process RIBOTISH_PREDICT {

    echo true

    label "ribotish"
    label "heavy_computation"

    publishDir "${params.ribotish_connect_outDir}/ribotish_predict", mode: 'copy'

    input:
    each(path(bam_sort_index_folder))
    path gtf_file
    path genome

    output:
    //path '*.ribotish.predict'
    path '*'

    script:
    """
    input=\$(basename ${bam_sort_index_folder})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    cp ${bam_sort_index_folder}/* .

    ribotish predict \
	-b \${prefix}.${bam_sort_index_folder} \
	-g ${gtf_file} \
	-f ${genome} \
	-o \${prefix}.ribotish_pred.txt

    : '
    ribotish predict \
	-b \${prefix}.transcripts.mapped.unique.a_site_profile.sorted.bam \
	-g ${gtf_file} \
	-f ${genome} \
	${params.ribotish_predict_mode} \
	-o \${prefix}.ribotish.predict \
	&> \${prefix}_ribotish_predict.log
    '
    """


}


workflow RIBOTISH_PIPE {

    take:
      gtf_ch
      bam_sort_index_folder_ch
      genome_ch

    main:

      RIBOTISH_QUALITY(  bam_sort_index_folder_ch,
			 gtf_ch  )

      ribotish_predict = RIBOTISH_PREDICT(  bam_sort_index_folder_ch,
					    gtf_ch,
					    genome_ch  )

    emit:
      ribotish_predict

}











