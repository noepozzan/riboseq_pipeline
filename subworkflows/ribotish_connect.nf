nextflow.enable.dsl=2

// https://github.com/zhpn1024/ribotish

process RIBOTISH_CONNECT {

    echo true

    label "ribotish"
    label "heavy_computation"

    publishDir "${params.ribotish_connect_outDir}/ribotish_connect", mode: 'copy'

    input:
    //each(path(bam_file))
    each(path(bam_sort_index_folder))
    path gtf_file
    path genome
    path genome_fai

    output:
    //path '*.ribotish.predict'
    path '*'

    script:
    """
    input=\$(basename ${bam_sort_index_folder})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    cp ${bam_sort_index_folder}/* .

    echo \$(pwd) > \${prefix}.esel

    ribotish predict \
	-b \${prefix}.transcripts.mapped.unique.a_site_profile.sorted.bam \
	-g ${gtf_file} \
	-f ${genome} \
	--longest \
	-o \${prefix}.ribotish.predict \
	&> \${prefix}_ribotish_predict.log

    """


}


workflow RIBOTISH_PIPE {

    take:
      //longest_pc_transcript_per_gene_gtf_ch
      gtf_ch
      //transcripts_mapped_unique_a_site_profile_sorted_bam_ch
      bam_sort_index_folder_ch
      genome_ch
      genome_fai_ch

    main:
      ribotish_predict = RIBOTISH_CONNECT(  //transcripts_mapped_unique_a_site_profile_sorted_bam_ch,
					    bam_sort_index_folder_ch,
					    //longest_pc_transcript_per_gene_gtf_ch,
					    gtf_ch,
					    genome_ch,
					    genome_fai_ch  )

    emit:
      ribotish_predict

}











