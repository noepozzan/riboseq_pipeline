#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TRIM_FIRST_BASES {

    label "cutadapt"

    publishDir "${params.reads_dir}/trim_first_bases", mode: 'copy'

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

    publishDir "${params.reads_dir}/clip_reads", mode: 'copy'

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

    publishDir "${params.reads_dir}/trim_reads", mode: 'copy'

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

    publishDir "${params.reads_dir}/filter_reads", mode: 'copy'

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

    publishDir "${params.reads_dir}/fastq_to_fasta", mode: 'copy'

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


workflow READS_PIPE {

    take:
    riboseq_reads_ch

    main:   
	// prepare riboseq reads
	TRIM_FIRST_BASES(riboseq_reads_ch)
	
	CLIP_READS(TRIM_FIRST_BASES.out)
	
	TRIM_READS(CLIP_READS.out)
	
	FILTER_READS(TRIM_READS.out)
	
	FASTQ_TO_FASTA(FILTER_READS.out)

    emit:
	FASTQ_TO_FASTA.out

}


