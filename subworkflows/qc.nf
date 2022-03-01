#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process COUNT_OLIGOS {
    
    label "htseq_biopython"

    publishDir "${params.riboseq_dir}/count_oligos", mode: 'copy'

    input:
    each(path(reads))
    path oligos
    path py_script

    output:
    path '*_oligos_counts'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    python ${py_script} --fastq <(gunzip -c ${reads}) \
        --oligos ${oligos} \
        --out \${prefix}_oligos_counts 2> \${prefix}_count_oligos.log

    """

}

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
		2> \${prefix}_trim_reads.log

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
		-o \${prefix}.pro_filtered 2> \${prefix}_filter_reads.log

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
		-o \${prefix}.pro_filtered_fasta 2> \${prefix}_fastq_to_fasta.log

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


process COUNT_OVERREPRESENTED_SEQUENCES_OTHER {
    
    label "pysam"

    publishDir "${params.riboseq_dir}/count_overrepresented_sequences_other", mode: 'copy'

    input:
    each(path(sam))
    path script_py

    output:
    path '*.overrepresented_sequences_counts', emit: sequences

    script:
    """
    input=\$(basename ${sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    python ${script_py} \
		--sam ${sam} \
		--out \${prefix}.overrepresented_sequences_counts 2> \
		\${prefix}_count_overrepresented_sequences_other.log
    
    : '
    grep -P -v \"^@\" ${sam} \
	| cut -f 10 | sort | uniq -c \
	| sort -n -r > \
	\${prefix}.overrepresented_sequences_counts \
	2> \${prefix}_count_overrepresented_sequences_other.log
    '

    """

}

process MAP_TRANSCRIPTOME_SEGEMEHL {
    
    label "segemehl"
    label "heavy_computation"

    publishDir "${params.riboseq_dir}/map_transcriptome_segemehl", mode: 'copy'

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
		${params.segemehl_silent} \
		-i ${index} \
		-d ${sequence} \
		-q ${reads} \
		--accuracy ${params.segemehl_accuracy} \
		--threads ${params.segemehl_threads} \
		-o \${prefix}.transcripts_mapped_sam \
		-u \${prefix}.transcripts_unmapped_fa \
		2> \${prefix}_map_to_transcripts.log
	
	"""

}

process REMOVE_MULTIMAPPERS {
    
    // basic container like alpine
    label "pysam"

    publishDir "${params.riboseq_dir}/remove_multimappers", mode: 'copy'

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

    publishDir "${params.riboseq_dir}/sam_to_bam_sort_and_index", mode: 'copy'

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

process READ_LENGTH_HISTOGRAM {
    
    label "rcrunch_python"

    publishDir "${params.riboseq_dir}/read_length_histogram", mode: 'copy'

    input:
    each(path(sam))
    path script_py

    output:
    path '*'

    script:
    """
    workd=\$(pwd)
    input=\$(basename ${sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    python ${script_py} \
	--sam ${sam} \
	--outdir \${prefix} \
	2> \${prefix}_read_length_histogram.log

    """

}

process DETERMINE_P_SITE_OFFSET {
    
    label "pysam"

    publishDir "${params.riboseq_dir}/determine_p_site_offset", mode: 'copy'

    input:
    each(path(bam_folder))
    path transcript_id_gene_id_CDS
    path script_py

    output:
    path '*.alignment_offset.json', emit: offsets

    script:
    """
    workd=\$(pwd)
    input=\$(basename ${bam_folder})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    echo \$workd
    echo \${prefix}

    cp ${bam_folder}/* .

    python ${script_py} \
		--bam \${prefix}.*.bam \
		--cds_coordinates ${transcript_id_gene_id_CDS} \
		--outdir \${prefix}_p_site_offset \
		2> \${prefix}_p_site_offset.log

    cp \${prefix}_p_site_offset/* .
    mv alignment_offset.json \${prefix}.alignment_offset.json

    """

}

process COUNT_READS {

    label "pysam"

    publishDir "${params.riboseq_dir}/count_reads", mode: 'copy'

    input:
    each(path(bam_folder))
    path transcript_id_gene_id_CDS
    each(path(p_site_offsets))
    path script_py

    output:
    path '*.counts.tsv' optional true

    script:
    """
    workd=\$(pwd)
    input1=\$(basename ${bam_folder})
    prefix1=\$(echo \$input1 | cut -d '.' -f 1)
    input2=\$(basename ${p_site_offsets})
    prefix2=\$(echo \$input2 | cut -d '.' -f 1)

    if [ "\$prefix1" == "\$prefix2" ]; then

	cp ${bam_folder}/* .

    	mkdir \${prefix1}.count_reads
	
	python ${script_py} \
    	--bam \${prefix1}.*.bam \
        --tsv ${transcript_id_gene_id_CDS} \
        --json ${p_site_offsets} \
        --outdir \${prefix1}.count_reads

    	cp \${prefix1}.count_reads/* .
    	mv counts.tsv \${prefix1}.counts.tsv

    fi

    """

}

process CHECK_PERIODICITY {

    label "rcrunch_python"

    publishDir "${params.riboseq_dir}/check_periodicity", mode: 'copy'

    input:
    each(path(bam_folder))
    path transcript_id_gene_id_CDS
    each(path(p_site_offsets))
    path script_py

    output:
    path '*.periodicity_start.pdf' optional true
    path '*.periodicity_stop.pdf' optional true
    path '*.Periodicity_Analysis_Start_Ribo_Seq.txt' optional true

    script:
    """
    workd=\$(pwd)
    input1=\$(basename ${bam_folder})
    prefix1=\$(echo \$input1 | cut -d '.' -f 1)
    input2=\$(basename ${p_site_offsets})
    prefix2=\$(echo \$input2 | cut -d '.' -f 1)

    if [ "\$prefix1" == "\$prefix2" ]; then

        cp ${bam_folder}/* .

        mkdir \${prefix1}.check_periodicity

	python ${script_py} \
    	--bam \${prefix1}.*.bam \
      	--tsv ${transcript_id_gene_id_CDS} \
       	--json ${p_site_offsets} \
       	--outdir \${prefix1}.check_periodicity \
       	--codnum ${params.check_peridocitiy_codnum} \
        2> \${prefix1}_check_periodicity.log

    	cp \${prefix1}.check_periodicity/* .
    	mv periodicity_start.pdf \${prefix1}.periodicity_start.pdf
		mv periodicity_stop.pdf \${prefix1}.periodicity_stop.pdf
		mv Periodicity_Analysis_Start_Ribo_Seq.txt \${prefix1}.Periodicity_Analysis_Start_Ribo_Seq.txt

    fi

    """

}

process FILTER_LENGTHS_OFFSETS {
    
    label "pysam"

    publishDir "${params.riboseq_dir}/filter_lengths_offsets", mode: 'copy'

    input:
    each(path(bam_folder))
    each(path(p_site_offsets))
    path script_py

    output:
    path '*.unique_a_site.bam' optional true

    script:
    """
    workd=\$(pwd)
    input1=\$(basename ${bam_folder})
    prefix1=\$(echo \$input1 | cut -d '.' -f 1)
    input2=\$(basename ${p_site_offsets})
    prefix2=\$(echo \$input2 | cut -d '.' -f 1)

    if [ "\$prefix1" == "\$prefix2" ]; then

        cp ${bam_folder}/* .

        python ${script_py} \
        	--bam \${prefix1}.*.bam \
            --p_site_offsets ${p_site_offsets} \
			--bam_out \${prefix1}.unique_a_site.bam \
            2> \${prefix1}_filter_lengths_offsets.log
	
    fi

    """

}

process BAM_SORT_AND_INDEX {
    
    label "samtools"

    publishDir "${params.riboseq_dir}/bam_sort_and_index", mode: 'copy'

    input:
    path bam

    output:
    path '*.bam'
    path '*.bam.bai'
    path '*.bam_sort_index'

    script:
    """
    input=\$(basename ${bam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    samtools sort ${bam} \
		> \${prefix}.unique_a_site_sorted.bam; \
		samtools index \${prefix}.unique_a_site_sorted.bam; \
		2> \${prefix}.bam_sort_and_index.log

    mkdir \${prefix}.bam_sort_index
    cp \${prefix}.unique* \${prefix}.bam_sort_index
    """

}




workflow QC_PIPE {

    take:
	genome_ch
    riboseq_reads_ch
    oligos_ch
    other_RNAs_sequence_ch
    longest_pc_transcript_per_gene_fa  
	transcript_id_gene_id_CDS
    gtf_ch

    main:   
	if ( params.run_count_oligos == "true" ) {
	COUNT_OLIGOS(
		reads_ch,
		oligos_ch,
		params.count_oligos_script
	)
    }

	// generate indexes for later mapping
	other_RNAs_index = SEGEMEHL_INDEX_rRNA(other_RNAs_sequence_ch)

	transcriptome_index = SEGEMEHL_INDEX_TRANSCRIPTOME(longest_pc_transcript_per_gene_fa)

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

	COUNT_OVERREPRESENTED_SEQUENCES_OTHER(
		other_genes_mapped_sam,
		params.find_overrepresented_sequences_script
	)
	overrepresented_sequences_other = COUNT_OVERREPRESENTED_SEQUENCES_OTHER.out.sequences
	
	// map the filtered reads to the transcriptome to then do QC
	MAP_TRANSCRIPTOME_SEGEMEHL(
		other_genes_unmapped_fasta,
		transcriptome_index,
		longest_pc_transcript_per_gene_fa
	)
	transcripts_mapped_sam = MAP_TRANSCRIPTOME_SEGEMEHL.out.mapped
	transcripts_unmapped_fasta = MAP_TRANSCRIPTOME_SEGEMEHL.out.unmapped

	transcripts_mapped_unique_sam = REMOVE_MULTIMAPPERS(transcripts_mapped_sam)
	
	SAM_TO_BAM_SORT_AND_INDEX(
		transcripts_mapped_unique_sam
	)
	transcripts_mapped_unique_sorted_bam = SAM_TO_BAM_SORT_AND_INDEX.out.bam
	transcripts_mapped_unique_sorted_bam_bai = SAM_TO_BAM_SORT_AND_INDEX.out.bai
	bam_bai_folder = SAM_TO_BAM_SORT_AND_INDEX.out.folder

	READ_LENGTH_HISTOGRAM(
		transcripts_mapped_unique_sam,
		params.plot_read_lengths_script
	)
      
	DETERMINE_P_SITE_OFFSET(
		bam_bai_folder,
		transcript_id_gene_id_CDS,
		params.determine_p_site_offsets_script
	)
    alignment_offset_json = DETERMINE_P_SITE_OFFSET.out.offsets
  
	COUNT_READS(
		bam_bai_folder,
		transcript_id_gene_id_CDS,
		alignment_offset_json,
		params.count_reads_script
	)
      
	CHECK_PERIODICITY(
		bam_bai_folder,
		transcript_id_gene_id_CDS,
		alignment_offset_json,
		params.check_periodicity_script
	)
      
	FILTER_LENGTHS_OFFSETS(
		bam_bai_folder,
		alignment_offset_json,
		params.filter_lengths_offsets_script
	)
    transcripts_mapped_unique_a_site_profile_bam = FILTER_LENGTHS_OFFSETS.out  
    
	BAM_SORT_AND_INDEX(
		transcripts_mapped_unique_a_site_profile_bam
	)	


}


