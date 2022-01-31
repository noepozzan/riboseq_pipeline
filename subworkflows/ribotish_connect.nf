nextflow.enable.dsl=2

// https://github.com/zhpn1024/ribotish

process RIBOTISH_QUALITY {

    label "ribotish"

    publishDir "${params.ribotish_connect_outDir}/ribotish_quality", mode: 'copy'

    input:
    each(path(bam_sort_index_folder))
    path gtf_file

    output:
	path '*.para.py', emit: offset
    path '*.pdf', emit: ribo_pdf
	path '*.txt', emit: ribo_txt

    script:
    """
    input=\$(basename ${bam_sort_index_folder})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    cp ${bam_sort_index_folder}/* .

    ribotish quality \
		-b \${prefix}.*.bam \
		-g ${gtf_file} \
		--th ${params.ribotish_quality_th}
    """

}

process RIBOTISH_PREDICT {

    label "ribotish"
    //label "heavy_computation"

    publishDir "${params.ribotish_connect_outDir}/ribotish_predict", mode: 'copy'

    input:
    each(path(bam_sort_index_folder))
    path gtf_file
    path genome

    output:
    path '*.ribotish.pred_all.txt', emit: ribo_pred

    script:
    """
    input=\$(basename ${bam_sort_index_folder})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    cp ${bam_sort_index_folder}/* .

	ribotish quality \
        -b \${prefix}.*.bam \
        -g ${gtf_file} \
        --th ${params.ribotish_quality_th}

    ribotish predict \
		-b \${prefix}.*.bam \
		-g ${gtf_file} \
		-f ${genome} \
		${params.ribotish_predict_mode} \
		-o \${prefix}.ribotish.pred.txt \
		&> \${prefix}_ribotish_predict.log

    """

}

process COMBINE {

    echo true

    //label "pandas"

    publishDir "${params.ribotish_connect_outDir}/ribotish_combine", mode: 'copy'

    input:
    path ribo_pred

    output:
	path 'combined_shORFs_position_preds.txt', emit: prediction

    script:
    """

	WRITE_HEADER="true"
    for VAR in ${ribo_pred}
    do

		if [ "\$WRITE_HEADER" == "true" ]; then

	    	head -1 \$VAR > combined_shORFs_position_preds.txt
	    	WRITE_HEADER="false"

		fi

		tail -n +2 -q \$VAR >> combined_shORFs_position_preds.txt

    done
	
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

      RIBOTISH_PREDICT(  bam_sort_index_folder_ch,
						 gtf_ch,
						 genome_ch  )
      
	  COMBINE(RIBOTISH_PREDICT.out.ribo_pred.collect())
	  ribo_pred = COMBINE.out.prediction

    emit:
      ribo_pred

}











