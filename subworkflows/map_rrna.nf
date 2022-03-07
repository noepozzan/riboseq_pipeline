#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SEGEMEHL_INDEX_rRNA {

    label "segemehl"

    publishDir "${params.map_dir}/segemehl_index_rrna", mode: 'copy'

    input:
    path sequence

    output:
    path 'segemehl_index_other_rnas.idx'

    script:
    """
    segemehl.x \
    	-x segemehl_index_other_rnas.idx \
    	-d ${sequence} \
    	2> segemehl_index_rrna.log

    """

}

process MAP_rRNA_SEGEMEHL {
    
    label "segemehl"

    publishDir "${params.map_dir}/map_rrna_segemehl", mode: 'copy'

    input:
    each(path(reads))
    path index
    path sequence

    output:
    path '*.other_genes_mapped_sam', emit: mapped
    path '*.other_genes_unmapped', emit: unmapped

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    segemehl.x \
		-i ${index} \
		-d ${sequence} \
		-q ${reads} \
		${params.segemehl_args} \
		-o \${prefix}.other_genes_mapped_sam \
		-u \${prefix}.other_genes_unmapped \
		&> \${prefix}_map_rRNA.log
	
    """

}


workflow rRNA_PIPE {

    take:
	genome_ch
    other_RNAs_sequence_ch
    fasta_reads_ch

    main:   
	// generate index for later mapping
	other_RNAs_index = SEGEMEHL_INDEX_rRNA(other_RNAs_sequence_ch)

	// map reads to rRNA, continue with unmapped reads
	MAP_rRNA_SEGEMEHL(
		fasta_reads_ch,
		other_RNAs_index,
		other_RNAs_sequence_ch
	)
	other_genes_mapped_sam = MAP_rRNA_SEGEMEHL.out.mapped
	other_genes_unmapped_fasta = MAP_rRNA_SEGEMEHL.out.unmapped

    emit:
	other_genes_unmapped_fasta
	other_genes_mapped_sam
}


