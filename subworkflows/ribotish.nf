nextflow.enable.dsl=2

// https://github.com/zhpn1024/ribotish

process RIBOTISH_QUALITY {

    label "ribotish"

    publishDir "${params.ribotish_dir}/ribotish_quality", mode: 'copy'

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
    label "heavy_computation"

    publishDir "${params.ribotish_dir}/ribotish_predict", mode: 'copy'

    input:
    each(path(bam_sort_index_folder))
    path gtf_file
    path genome
	each(path(offsets))

    output:
    path '*.ribotish_pred.txt', emit: ribo_pred, optional: true

    script:
    """
	input1=\$(basename ${bam_sort_index_folder})
    prefix1=\$(echo \$input1 | cut -d '.' -f 1)
    input2=\$(basename ${offsets})
    prefix2=\$(echo \$input2 | cut -d '.' -f 1)

    if [ "\$prefix1" == "\$prefix2" ]; then

		cp ${bam_sort_index_folder}/* .

		ribotish predict \
			-b \${prefix1}.*.bam \
			-g ${gtf_file} \
			-f ${genome} \
			${params.ribotish_predict_mode} \
			-o \${prefix1}.ribotish_pred.txt \
			&> \${prefix1}_ribotish_predict.log

	fi

	"""

}

process GFFREAD {

	label "gffread"

	publishDir "${params.ribotish_dir}/gffread", mode: 'copy'

	input:
	path genome
	path genome_fai
	path gtf

	output:
	path 'transcripts.fa', emit: fasta	

	script:
	"""

	gffread \
        -w transcripts.fa \
        -g ${genome} \
        ${gtf}

	"""

}

process SORF_TO_PEPTIDE {

	//label "python..."

    publishDir "${params.ribotish_dir}/ribotish_combine", mode: 'copy'

    input:
    each(path(ribo_pred))
	path fasta
	path python_script

    output:
	path '*.speptide', emit: prediction

    script:
    """
	input=\$(basename ${ribo_pred})
    prefix=\$(echo \$input | cut -d '.' -f 1)

	python ${python_script} \
		--ribo_pred ${ribo_pred} \
		--fasta ${fasta} \
		--out \${prefix}.speptide

	
	"""

}


workflow RIBOTISH_PIPE {

    take:
    gtf_ch
    bam_sort_index_folder_ch
    genome_ch
	genome_fai_ch
	

    main:
    RIBOTISH_QUALITY(
		bam_sort_index_folder_ch,
		gtf_ch
	)

    RIBOTISH_PREDICT(
		bam_sort_index_folder_ch,
		gtf_ch,
		genome_ch,
		RIBOTISH_QUALITY.out.offset
	)
    
	GFFREAD(
		genome_ch,
		genome_fai_ch,
		gtf_ch
	)
	transcripts_fa = GFFREAD.out.fasta

	SORF_TO_PEPTIDE(
		RIBOTISH_PREDICT.out.ribo_pred,
		transcripts_fa,
		params.sorf_peptide_script
	)
	speptide = SORF_TO_PEPTIDE.out.prediction

    emit:
    speptide

}











