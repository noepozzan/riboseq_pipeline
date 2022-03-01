#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SEGEMEHL_INDEX_rRNA {

    label "segemehl"

    publishDir "${params.riboseq_dir}/segemehl_index_rrna", mode: 'copy'

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

process SEGEMEHL_INDEX_TRANSCRIPTOME {

    label "segemehl"

    publishDir "${params.riboseq_dir}/segemehl_index_transcriptome", mode: 'copy'

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

process STAR_INDEX_GENOME {

    label 'star'

    publishDir "${params.riboseq_dir}/star_index_transcriptome", mode: 'copy'

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

process TRIM_FIRST_BASES {

    label "cutadapt"

    publishDir "${params.riboseq_dir}/trim_first_bases", mode: 'copy'

    input:
    path reads

    output:
    path '*.trimmed_first_bases'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    (cutadapt \
		--cut ${params.cut} \
		--minimum-length ${params.minimum_length} \
		${reads} | gzip > \
		\${prefix}.trimmed_first_bases) \
		&> \${prefix}_trim_first_bases.log

    """

}

process CLIP_READS {
    
    label "fastx"

    publishDir "${params.riboseq_dir}/clip_reads", mode: 'copy'

    input:
    path reads

    output:
    path '*.pro_clipped'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastx_clipper \
		${params.clip_reads_args} \
		-i <(zcat ${reads}) \
		-o \${prefix}.pro_clipped \
		&> \${prefix}_clip_reads.log
    
    """

}


process TRIM_READS {
    
	//errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	//maxRetries 5

    label "fastx"

    publishDir "${params.riboseq_dir}/trim_reads", mode: 'copy'

    input:
    path reads

    output:
    path '*.pro_trimmed'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastq_quality_trimmer \
		${params.trim_reads_args} \
		-i <(zcat ${reads}) \
		-o \${prefix}.pro_trimmed \
		&> \${prefix}_trim_reads.log

    """

}

process FILTER_READS {
    
    label "fastx"

    publishDir "${params.riboseq_dir}/filter_reads", mode: 'copy'

    input:
    path reads

    output:
    path '*.pro_filtered'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastq_quality_filter \
		${params.filter_reads_args} \
		-i <(zcat ${reads}) \
		-o \${prefix}.pro_filtered \
		&> \${prefix}_filter_reads.log

    """

}

process FASTQ_TO_FASTA {
    
    label "fastx"

    publishDir "${params.riboseq_dir}/fastq_to_fasta", mode: 'copy'

    input:
    path reads

    output:
    path '*.pro_filtered_fasta'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastq_to_fasta \
		${params.fastq_to_fasta_args} \
		-i <(zcat ${reads}) \
		-o \${prefix}.pro_filtered_fasta \
		&> \${prefix}_fastq_to_fasta.log

    """

}

process MAP_rRNA_SEGEMEHL {
    
    label "segemehl"
    label "heavy_computation"

    publishDir "${params.riboseq_dir}/map_rrna_segemehl", mode: 'copy'

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
		${params.segemehl_silent} \
		-i ${index} \
		-d ${sequence} \
		-q ${reads} \
		--accuracy ${params.segemehl_accuracy} \
		--threads ${params.segemehl_threads} \
		-o \${prefix}.other_genes_mapped_sam \
		-u \${prefix}.other_genes_unmapped \
		&> \${prefix}_map_rRNA.log
	
    """

}

process MAP_GENOME_STAR {

    label "star"
    label "heavy_computation"

    publishDir "${params.riboseq_dir}/map_genome_star", mode: 'copy'

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

    publishDir "${params.riboseq_dir}/sam_to_bam_sort_and_index_star", mode: 'copy'

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

workflow RIBOSEQ_PIPE {

    take:
	genome_ch
    riboseq_reads_ch
    other_RNAs_sequence_ch
    gtf_ch

    main:   
	// generate indexes for later mapping
	other_RNAs_index = SEGEMEHL_INDEX_rRNA(other_RNAs_sequence_ch)
	genome_index = STAR_INDEX_GENOME(genome_ch)

	// prepare riboseq reads
	trimmed_first_bases_fastq = TRIM_FIRST_BASES(riboseq_reads_ch)
	pro_clipped_fastq = CLIP_READS(trimmed_first_bases_fastq)
	pro_trimmed_fastq = TRIM_READS(pro_clipped_fastq)
	pro_filtered_fastq = FILTER_READS(pro_trimmed_fastq)
	pro_filtered_fasta = FASTQ_TO_FASTA(pro_filtered_fastq)

	// map reads to rRNA, continue with unmapped reads
	MAP_rRNA_SEGEMEHL(
		pro_filtered_fasta,
		other_RNAs_index,
		other_RNAs_sequence_ch
	)
	other_genes_mapped_sam = MAP_rRNA_SEGEMEHL.out.mapped
	other_genes_unmapped_fasta = MAP_rRNA_SEGEMEHL.out.unmapped

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


