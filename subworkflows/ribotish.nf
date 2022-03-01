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

	echo true

    label "ribotish"
    label "heavy_computation"

    publishDir "${params.ribotish_dir}/ribotish_predict", mode: 'copy'

    input:
	each(path(bam_sort_index_folder))
    path gtf
    path genome
	each(path(offsets))
 
    output:
    path '*.ribotish_pred_all.txt', emit: ribo_pred, optional: true

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
			-g ${gtf} \
			-f ${genome} \
			${params.ribotish_predict_mode} \
			-o \${prefix1}.ribotish_pred.txt \
			&> \${prefix1}_ribotish_predict.log

	fi
	
	"""

}

/*
process RIBOTISH_PREDICT {

	echo true

    label "ribotish"
    label "heavy_computation"

    publishDir "${params.ribotish_dir}/ribotish_predict", mode: 'copy'

    input:
	tuple val(sample_id), path(reads)
	each(path(gtf))
	each(path(genome))
 
    output:
    path '*.ribotish_pred_all.txt', emit: ribo_pred, optional: true

    script:
    """

	cp ${reads[1]}/* .

	ribotish predict \
		-b *.bam \
		-g ${gtf} \
		-f ${genome} \
		--ribopara ${reads[0]} \
		${params.ribotish_predict_mode} \
		-o ${sample_id}.ribotish_pred.txt \
		&> ${sample_id}_ribotish_predict.log

	"""

}
*/

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

	label "sorf_to_speptide"

    publishDir "${params.ribotish_dir}/sorf_to_peptide", mode: 'copy'

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

	python3 ${python_script} \
		--ribo_pred ${ribo_pred} \
		--fasta ${fasta} \
		--out \${prefix}.speptide

	
	"""

}

process COMBINE {
    
	label "sorf_to_speptide"

    publishDir "${params.ribotish_dir}/combine", mode: 'copy'

    input:
    path ribo_pred

    output:
	path 'combined_speptide.txt', emit: combined_prediction

    script:
    """
	WRITE_HEADER="true"
    for VAR in ${ribo_pred}
    do
		if [ "\$WRITE_HEADER" == "true" ]; then
	    	head -1 \$VAR > combined_speptide.txt
	    	WRITE_HEADER="false"
		fi
		tail -n +2 -q \$VAR >> combined_speptide.txt
    done
	
	"""

}


process MOVE_RENAME {

	echo true

	publishDir "${params.ribotish_dir}/move_rename", mode: 'copy'

	input:
	path files1
	path files2
	path python_script

	output:
	path "Sample*"

	shell:
	params.run_predict = true
	'''
	workd=$(pwd)
	
	python !{python_script} \
		--dir $workd

	echo $workd

	: '
	# sort file lists
	files1_sorted=$(echo !{files1} | xargs -n1 | sort | xargs)
	files2_sorted=$(echo !{files2} | xargs -n1 | sort | xargs)

	# transform to array
	read -a arr1 <<< "$files1_sorted"
	read -a arr2 <<< "$files2_sorted"

	# go through arrays and join
	for i in "${!arr1[@]}"
	do
		echo ["${arr1[$i]}" "${arr2[$i]}"]
	done
	'
	'''
	
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
	offsets = RIBOTISH_QUALITY.out.offset

	/*
	MOVE_RENAME(
		RIBOTISH_QUALITY.out.offset.collect(),
		bam_sort_index_folder_ch.collect(),
		params.rename_script
	)
	
	params.run_predict = false
	if ( params.run_predict == true ){
	reads_ch = Channel.fromFilePairs("${projectDir}/results/ribotish/move_rename/Sample_*_{1,2}", type: 'any', checkIfExists: true)
	
	RIBOTISH_PREDICT(
        reads_ch,
		gtf_ch,
        genome_ch
    )
	}
	*/

	RIBOTISH_PREDICT(
		bam_sort_index_folder_ch,
		gtf_ch,
		genome_ch,
		offsets
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

	COMBINE(
		speptide.collect()
	)
	speptide_combined = COMBINE.out.combined_prediction

	emit:
	speptide_combined

}



