#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process STAR_INDEX_GENOME {

    label 'star'

    publishDir "${params.map_dir}/star_index_transcriptome", mode: 'copy'

    input:
    path sequence

    output:
    path 'starIndex'

    script:
    """
    mkdir starIndex

    STAR --runThreadN ${params.index_threads} \
    	--runMode genomeGenerate \
        --genomeDir starIndex \
        --genomeFastaFiles ${sequence}

    """

}

process MAP_GENOME_STAR {

    label "star"
    label "heavy_computation"

    publishDir "${params.map_dir}/map_genome_star", mode: 'copy'

    input:
    each(path(reads))
    path index
    path gtf

    output:
    path '*.Aligned.out.sam', emit: aligned
    path '*.Unmapped*', emit: unmapped

    script:
    """
    for VAR in ${reads}
    do

        input=\$(basename \$VAR)
        prefix=\$(echo \$input | cut -d '.' -f 1)

		STAR --runThreadN ${params.star_map_threads} \
			--genomeDir ${index} \
			--sjdbGTFfile ${gtf} \
			--outSAMattributes All \
			--quantMode GeneCounts \
			--readFilesIn \$VAR \
			--outReadsUnmapped Fastx \
			--outFileNamePrefix \${prefix}.

	: '
		--outFilterMismatchNmax 2 \
		--alignEndsType EndToEnd \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		--alignIntronMax 20000 \
		--outMultimapperOrder Random \
		--outSAMmultNmax 1

	mv *.Unmapped.out.* \${prefix}.unmapped
	'

    done

    """

}


process SAM_TO_BAM_SORT_AND_INDEX_STAR {

    label "samtools"

    publishDir "${params.map_dir}/sam_to_bam_sort_and_index_star", mode: 'copy'

    input:
    path sam

    output:
    path '*.sorted_indexed.bam', emit: bam
    path '*.sorted_indexed.bam.bai', emit: bai
    path '*.folder_sorted_indexed_bam', emit: folder

    script:
    """
    input=\$(basename ${sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    samtools view -bS ${sam} \
        | samtools sort - > \${prefix}.sorted_indexed.bam; \
        samtools index \${prefix}.sorted_indexed.bam; \
        2> \${prefix}_bam_sorted_and_indexed.log

    mkdir \${prefix}.folder_sorted_indexed_bam
    cp \${prefix}.sorted* \${prefix}.folder_sorted_indexed_bam

    """

}

workflow GENOME_PIPE {

    take:
	genome_ch
    gtf_ch
	other_genes_unmapped_fasta

    main:   
	// generate index for later mapping
	genome_index = STAR_INDEX_GENOME(genome_ch)

	// map reads to genome with star and use samtools
	MAP_GENOME_STAR(
		other_genes_unmapped_fasta,
		genome_index,
		gtf_ch
	)
	star_mapped_sam = MAP_GENOME_STAR.out.aligned
	star_unmapped_fasta = MAP_GENOME_STAR.out.unmapped

	SAM_TO_BAM_SORT_AND_INDEX_STAR(
        star_mapped_sam
    )
    star_mapped_unique_sorted_bam = SAM_TO_BAM_SORT_AND_INDEX_STAR.out.bam
    star_mapped_unique_sorted_bam_bai = SAM_TO_BAM_SORT_AND_INDEX_STAR.out.bai
    bam_sort_index_folder = SAM_TO_BAM_SORT_AND_INDEX_STAR.out.folder

    emit:
	bam_sort_index_folder

}

