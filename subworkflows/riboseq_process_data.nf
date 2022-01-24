#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process COUNT_OLIGOS {
    
    label "htseq_biopython"

    publishDir "${params.riboseq_process_data_outDir}/count_oligos", mode: 'copy'

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

process TRIM_FIRST_BASES {

    label "cutadapt"

    publishDir "${params.riboseq_process_data_outDir}/trim_first_bases", mode: 'copy'

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

    publishDir "${params.riboseq_process_data_outDir}/clip_reads", mode: 'copy'

    input:
    path reads

    output:
    path '*.pro_clipped'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastx_clipper \
	${params.clip_reads_v} \
	${params.clip_reads_n} \
	-l ${params.clip_reads_l} \
	${params.clip_reads_c} \
	${params.clip_reads_z} \
	-a ${params.clip_reads_adapter} \
	-i <(zcat ${reads}) \
	-o \${prefix}.pro_clipped \
	&> \${prefix}_clip_reads.log
    
    """

}


process TRIM_READS {
    
    label "fastx"

    publishDir "${params.riboseq_process_data_outDir}/trim_reads", mode: 'copy'

    input:
    path reads

    output:
    path '*.pro_trimmed'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastq_quality_trimmer \
	${params.trim_reads_v} \
	-l ${params.trim_reads_l} \
	-t ${params.trim_reads_t} \
	-Q ${params.trim_reads_Q} \
	${params.trim_reads_z} \
	-i <(zcat ${reads}) \
	-o \${prefix}.pro_trimmed 2> \${prefix}_trim_reads.log

    """

}

process FILTER_READS {
    
    label "fastx"

    publishDir "${params.riboseq_process_data_outDir}/filter_reads", mode: 'copy'

    input:
    path reads

    output:
    path '*.pro_filtered'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastq_quality_filter \
	${params.filter_reads_v} \
	-q ${params.filter_reads_q} \
	-p ${params.filter_reads_p} \
	-Q ${params.filter_reads_Q} \
	${params.filter_reads_z} \
	-i <(zcat ${reads}) \
	-o \${prefix}.pro_filtered 2> \${prefix}_filter_reads.log

    """

}

process FASTQ_TO_FASTA {
    
    label "fastx"

    publishDir "${params.riboseq_process_data_outDir}/fastq_to_fasta", mode: 'copy'

    input:
    path reads

    output:
    path '*.pro_filtered_fasta'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    fastq_to_fasta \
	${params.fq_to_fa_v} \
	${params.fq_to_fa_n} \
	${params.fq_to_fa_r} \
	-i <(zcat ${reads}) \
	-o \${prefix}.pro_filtered_fasta 2> \${prefix}_fastq_to_fasta.log

    """

}

process MAP_TO_OTHER_GENES_SEGEMEHL {
    
    label "segemehl"
    label "heavy_computation"

    publishDir "${params.riboseq_process_data_outDir}/map_to_other_genes_segemehl", mode: 'copy'

    input:
    each(path(reads))
    path index
    path sequence

    output:
    path '*.other_genes_mapped_sam'
    path '*.other_genes_unmapped'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    segemehl.x ${params.map_to_other_genes_silent} \
	-i ${index} \
	-d ${sequence} \
	-q ${reads} \
	--accuracy ${params.map_to_other_genes_accuracy} \
	--threads ${params.map_to_other_genes_threads} \
	-o \${prefix}.other_genes_mapped_sam \
	-u \${prefix}.other_genes_unmapped \
	&> \${prefix}_map_to_other_genes.log
	
    """

}

process MAP_TO_OTHER_GENES_STAR {

    label "star"
    label "heavy_computation"

    publishDir "${params.riboseq_process_data_outDir}/map_to_other_genes_star", mode: 'copy'

    input:
    each(path(reads))
    path index

    output:
    path '*.Aligned.out.sam'
    path '*.Unmapped.out*'

    script:
    """
    for VAR in ${reads}
    do

        input=\$(basename \$VAR)
        prefix=\$(echo \$input | cut -d '.' -f 1)

        STAR --runThreadN ${params.map_to_other_genes_threads} \
                --genomeDir ${index} \
                --outFileNamePrefix \${prefix}. \
                --readFilesIn \$VAR \
                --outSAMattributes All \
		--outReadsUnmapped Fastx

    done

    """


}


process COUNT_OVERREPRESENTED_SEQUENCES_OTHER {
    
    label "pysam"

    publishDir "${params.riboseq_process_data_outDir}/count_overrepresented_sequences_other", mode: 'copy'

    input:
    each(path(sam))
    path script_py

    output:
    path '*.overrepresented_sequences_counts'

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

process MAP_TO_TRANSCRIPTS_SEGEMEHL {
    
    label "segemehl"
    label "heavy_computation"

    publishDir "${params.riboseq_process_data_outDir}/map_to_transcripts_segemehl", mode: 'copy'

    input:
    each(path(reads))
    path index
    path sequence

    output:
    path '*.transcripts_mapped_sam'
    path '*.transcripts_unmapped_fa'

    script:
    """
    input=\$(basename ${reads})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    segemehl.x \
	${params.map_to_transcripts_silent} \
	-i ${index} \
	-d ${sequence} \
	-q ${reads} \
	--accuracy ${params.map_to_transcripts_accuracy} \
	--threads ${params.map_to_transcripts_threads} \
	-o \${prefix}.transcripts_mapped_sam \
	-u \${prefix}.transcripts_unmapped_fa 2> \${prefix}_map_to_transcripts.log
    """

}

process MAP_TO_TRANSCRIPTS_STAR {

    label "star"
    label "heavy_computation"

    publishDir "${params.riboseq_process_data_outDir}/map_to_transcripts_star", mode: 'copy'

    input:
    each(path(reads))
    path index
    path gtf

    output:
    path '*.Aligned.out.sam'
    path '*.unmapped'

    script:
    """
    for VAR in ${reads}
    do

        input=\$(basename \$VAR)
        prefix=\$(echo \$input | cut -d '.' -f 1)

	STAR --runThreadN 12 \
		--genomeDir ${index} \
		--sjdbGTFfile ${gtf} \
		--outSAMattributes All \
		--quantMode GeneCounts \
		--outSAMtype BAM SortedByCoordinate \
		--readFilesIn \$VAR \
		--outFileNamePrefix \${prefix}.

	: '
	STAR --runThreadN ${params.map_to_transcripts_threads} \
                --genomeDir ${index} \
                --outFileNamePrefix \${prefix}. \
                --readFilesIn \$VAR \
                --outSAMattributes All \
		--outReadsUnmapped Fastx \
		--outFilterMismatchNmax 2 \
		--alignEndsType EndToEnd \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		--alignIntronMax 20000 \
		--outMultimapperOrder Random \
		--outSAMmultNmax 1
	'

	mv *.Unmapped.out.* \${prefix}.unmapped

    done

    """

}

process REMOVE_MULTIMAPPERS {
    
    // basic container like alpine
    label "pysam"

    publishDir "${params.riboseq_process_data_outDir}/remove_multimappers", mode: 'copy'

    input:
    path transcripts_mapped_sam

    output:
    path '*.transcripts_mapped_unique_sam'

    script:
    """
    input=\$(basename ${transcripts_mapped_sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    grep -P \"^@|\tNH:i:1\t\" ${transcripts_mapped_sam} \
	> \${prefix}.transcripts_mapped_unique_sam \
	2> \${prefix}_remove_multimappers.log

    """

}

process SAM2BAM_SORT_AND_INDEX {
    
    label "samtools"

    publishDir "${params.riboseq_process_data_outDir}/sam2bam_sort_and_index", mode: 'copy'

    input:
    //path sam
    path bam

    output:
    path '*.transcripts_mapped_unique_sorted_bam'
    path '*.transcripts_mapped_unique_sorted_bam.bai'
    path '*.sam2bam'

    script:
    """
    input=\$(basename ${sam})
    prefix=\$(echo \$input | cut -d '.' -f 1)


    samtools index ${bam}

    mkdir \${prefix}.sam2bam
    cp \${prefix}.transcripts* \${prefix}.sam2bam

    : '
    samtools view -bS ${sam} \
	| samtools sort - > \${prefix}.transcripts_mapped_unique_sorted_bam; \
	samtools index \${prefix}.transcripts_mapped_unique_sorted_bam; \
	2> \${prefix}_sam2bam_sort_and_index.log

    mkdir \${prefix}.sam2bam
    cp \${prefix}.transcripts* \${prefix}.sam2bam
    '
    """

}

process READ_LENGTH_HISTOGRAM {
    
    label "rcrunch_python"

    publishDir "${params.riboseq_process_data_outDir}/read_length_histogram", mode: 'copy'

    input:
    each(path(sam))
    path script_py

    output:
    path '*E14'

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

    publishDir "${params.riboseq_process_data_outDir}/determine_p_site_offset", mode: 'copy'

    input:
    //each(path(bam))
    each(path(bam_folder))
    path transcript_id_gene_id_CDS
    path script_py

    output:
    path '*.alignment_offset.json'

    script:
    """
    workd=\$(pwd)
    input=\$(basename ${bam_folder})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    echo \$workd
    echo \${prefix}

    cp ${bam_folder}/* .

    python ${script_py} \
	--bam \${prefix}.transcripts_mapped_unique_sorted_bam \
	--cds_coordinates ${transcript_id_gene_id_CDS} \
	--outdir \${prefix}_determine_p_site_offset \
	2> \${prefix}_determine_p_site_offset.log

    cp \${prefix}_determine_p_site_offset/* .
    mv alignment_offset.json \${prefix}.alignment_offset.json

    """

}

process COUNT_READS {

    label "pysam"

    publishDir "${params.riboseq_process_data_outDir}/count_reads", mode: 'copy'

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
                --bam \${prefix1}.transcripts_mapped_unique_sorted_bam \
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

    publishDir "${params.riboseq_process_data_outDir}/check_periodicity", mode: 'copy'

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
        	--bam \${prefix1}.transcripts_mapped_unique_sorted_bam \
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

process FILTER_READS_BASED_ON_READ_LENGTHS_AND_OFFSETS {
    
    label "pysam"

    publishDir "${params.riboseq_process_data_outDir}/filter_reads_based_on_read_lengths_and_offsets", mode: 'copy'

    input:
    each(path(bam_folder))
    each(path(p_site_offsets))
    path script_py

    output:
    path '*.transcripts.mapped.unique.a_site_profile.bam' optional true

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
                --bam \${prefix1}.transcripts_mapped_unique_sorted_bam \
                --p_site_offsets ${p_site_offsets} \
		--bam_out \${prefix1}.transcripts.mapped.unique.a_site_profile.bam \
                2> \${prefix1}_filter_reads_based_on_read_lengths_and_offsets.log
	
    fi

    """

}

process BAM_SORT_AND_INDEX {
    
    label "samtools"

    publishDir "${params.riboseq_process_data_outDir}/bam_sort_and_index", mode: 'copy'

    input:
    path bam

    output:
    path '*.transcripts.mapped.unique.a_site_profile.sorted.bam'
    path '*.transcripts.mapped.unique.a_site_profile.sorted.bam.bai'
    path '*.bam_sort_index'

    script:
    """
    input=\$(basename ${bam})
    prefix=\$(echo \$input | cut -d '.' -f 1)

    samtools sort ${bam} > \${prefix}.transcripts.mapped.unique.a_site_profile.sorted.bam; \
	samtools index \${prefix}.transcripts.mapped.unique.a_site_profile.sorted.bam; 2> \
	\${prefix}.bam_sort_and_index.log

    mkdir \${prefix}.bam_sort_index
    cp \${prefix}.transcripts* \${prefix}.bam_sort_index
    """

}




workflow RIBOSEQ_PROCESS_DATA_PIPE {

    take:
      pull_containers_ch
      reads_ch
      oligos_ch
      other_RNAs_sequence_ch
      count_oligos_script_ch
      other_RNAs_index_ch
      transcripts_index_ch
      find_overrepresented_sequences_script_ch
      transcripts_sequence_ch
      plot_read_lengths_script_ch
      transcript_id_gene_id_CDS_ch
      determine_p_site_offsets_script_ch
      count_reads_script_ch
      check_periodicity_script_ch
      filter_reads_based_on_read_lengths_and_offsets_script_ch
      gtf_ch

    main:
      if ( params.run_count_oligos == "true" ) {
	  oligos_counts = COUNT_OLIGOS(  reads_ch,
					 oligos_ch,
					 count_oligos_script_ch  )
      }

      trimmed_first_bases_fastq = TRIM_FIRST_BASES(reads_ch)
      pro_clipped_fastq = CLIP_READS(trimmed_first_bases_fastq)
      pro_trimmed_fastq = TRIM_READS(pro_clipped_fastq)
      pro_filtered_fastq = FILTER_READS(pro_trimmed_fastq)
      pro_filtered_fasta = FASTQ_TO_FASTA(pro_filtered_fastq)

      if ( params.aligner_other_genes == "segemehl" ) {
	  (other_genes_mapped_sam, other_genes_unmapped_fasta) = MAP_TO_OTHER_GENES_SEGEMEHL(  pro_filtered_fasta,
											     other_RNAs_index_ch,
											     other_RNAs_sequence_ch  )
}


      if ( params.aligner_other_genes == "star" ) {
	  (other_genes_mapped_sam, other_genes_unmapped_fasta) = MAP_TO_OTHER_GENES_STAR(  pro_filtered_fasta,
											   other_RNAs_index_ch  )
      }

      overrepresented_sequences_other = COUNT_OVERREPRESENTED_SEQUENCES_OTHER(  other_genes_mapped_sam,
										find_overrepresented_sequences_script_ch  )

      if ( params.aligner_genome == "segemehl" ) {
          (transcripts_mapped_sam, transcripts_unmapped_fasta) = MAP_TO_TRANSCRIPTS_SEGEMEHL(other_genes_unmapped_fasta, transcripts_index_ch, transcripts_sequence_ch)
      }

      if ( params.aligner_genome == "star" ) {
          (transcripts_mapped_sorted_bam, transcripts_unmapped_fasta) = MAP_TO_TRANSCRIPTS_STAR(other_genes_unmapped_fasta, transcripts_index_ch, gtf_ch)
      }

      if ( params.aligner_genome == "segemehl" ) {
	  transcripts_mapped_unique_sam = REMOVE_MULTIMAPPERS(transcripts_mapped_sam)
}

      if ( params.aligner_genome == "segemehl" ) {
	  (transcripts_mapped_unique_sorted_bam, transcripts_mapped_unique_sorted_bam_bai, bam_bai_folder) = SAM2BAM_SORT_AND_INDEX(transcripts_mapped_unique_sam)    
}

      if ( params.aligner_genome == "star" ) {
	  (transcripts_mapped_unique_sorted_bam, transcripts_mapped_unique_sorted_bam_bai, bam_sort_index_folder) = SAM2BAM_SORT_AND_INDEX(transcripts_mapped_sorted_bam)
}

      if ( params.aligner_genome == "segemehl" ) {
      read_length_histogram_pdf = READ_LENGTH_HISTOGRAM(  transcripts_mapped_unique_sam,
							  plot_read_lengths_script_ch  )
      
      alignment_offset_json = DETERMINE_P_SITE_OFFSET(  bam_bai_folder,
							transcript_id_gene_id_CDS_ch,
							determine_p_site_offsets_script_ch  )
      
      counts_tsv = COUNT_READS(  bam_bai_folder,
				 transcript_id_gene_id_CDS_ch,
				 alignment_offset_json,
				 count_reads_script_ch  )
      
      (periodicity_start_pdf, periodicity_stop_pdf, periodicity_analysis_start_ribo_seq) = CHECK_PERIODICITY(bam_bai_folder, transcript_id_gene_id_CDS_ch, alignment_offset_json, check_periodicity_script_ch)
      
      transcripts_mapped_unique_a_site_profile_bam = FILTER_READS_BASED_ON_READ_LENGTHS_AND_OFFSETS(bam_bai_folder, alignment_offset_json, filter_reads_based_on_read_lengths_and_offsets_script_ch)
      
      (transcripts_mapped_unique_a_site_profile_sorted_bam, transcripts_mapped_unique_a_site_profile_sorted_bam_bai, bam_sort_index_folder_filtered) = BAM_SORT_AND_INDEX(transcripts_mapped_unique_a_site_profile_bam)

}     

    emit:
      oligos_counts
      bam_sort_index_folder

}


