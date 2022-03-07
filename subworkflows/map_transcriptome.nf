#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process SEGEMEHL_INDEX_TRANSCRIPTOME {

    label "segemehl"

    publishDir "${params.map_dir}/segemehl_index_transcriptome", mode: 'copy'

    input:
    path sequence

    output:
    path 'segemehl_index_transcriptome'

    script:
    """
    segemehl.x \
    	-x segemehl_index_transcriptome \
    	-d ${sequence} \
    	2> segemehl_index_transcriptome.log

    """

}

process MAP_TRANSCRIPTOME_SEGEMEHL {
    
    label "segemehl"

    publishDir "${params.map_dir}/map_transcriptome_segemehl", mode: 'copy'

    input:
    each(path(reads))
    path index
    path sequence

    output:
    path '*.transcripts_mapped_sam', emit: mapped
    path '*.transcripts_unmapped_fa', emit: unmapped

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    segemehl.x \
		-i ${index} \
		-d ${sequence} \
		-q ${reads} \
		${params.segemehl_args} \
		-o \${prefix}.transcripts_mapped_sam \
		-u \${prefix}.transcripts_unmapped_fa \
		2> \${prefix}_map_to_transcripts.log
	
	"""

}

process REMOVE_MULTIMAPPERS {
    
    // basic container like alpine
    label "pysam"

    publishDir "${params.map_dir}/remove_multimappers", mode: 'copy'

    input:
    path transcripts_mapped_sam

    output:
    path '*.transcripts_mapped_unique_sam'

    script:
    """
    input=\$(basename ${transcripts_mapped_sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    grep -P \"^@|\tNH:i:1\t\" \
		${transcripts_mapped_sam} \
		> \${prefix}.transcripts_mapped_unique_sam \
		2> \${prefix}_remove_multimappers.log

    """

}

process SAM_TO_BAM_SORT_AND_INDEX {
    
    label "samtools"

    publishDir "${params.map_dir}/sam_to_bam_sort_and_index", mode: 'copy'

    input:
    path sam

    output:
    path '*.sorted_indexed.bam', emit: bam
    path '*.sorted_indexed.bam.bai', emit: bai
    path '*.folder_sorted_indexed_bam',emit: folder

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


workflow TRANSCRIPTOME_PIPE {

    take:
    longest_pc_transcript_per_gene_fa  
    other_genes_unmapped_fasta

    main:   
	// generate index for later mapping
	transcriptome_index = SEGEMEHL_INDEX_TRANSCRIPTOME(longest_pc_transcript_per_gene_fa)

	// map the filtered reads to the transcriptome to then do QC
	MAP_TRANSCRIPTOME_SEGEMEHL(
		other_genes_unmapped_fasta,
		transcriptome_index,
		longest_pc_transcript_per_gene_fa
	)
	transcripts_mapped_sam = MAP_TRANSCRIPTOME_SEGEMEHL.out.mapped
	transcripts_unmapped_fasta = MAP_TRANSCRIPTOME_SEGEMEHL.out.unmapped

	REMOVE_MULTIMAPPERS(
		transcripts_mapped_sam
	)
	transcripts_mapped_unique_sam = REMOVE_MULTIMAPPERS.out
	
	SAM_TO_BAM_SORT_AND_INDEX(
		transcripts_mapped_unique_sam
	)
	transcripts_mapped_unique_sorted_bam = SAM_TO_BAM_SORT_AND_INDEX.out.bam
	transcripts_mapped_unique_sorted_bam_bai = SAM_TO_BAM_SORT_AND_INDEX.out.bai
	bam_bai_folder = SAM_TO_BAM_SORT_AND_INDEX.out.folder

	emit:
	transcripts_mapped_unique_sam
	bam_bai_folder
}


