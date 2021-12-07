#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process COUNT_OLIGOS {
    
    label "htseq_biopython"

    publishDir "${params.riboseq_process_data_outDir}/count_oligos", mode: 'copy'

    input:
    path reads
    path oligos
    path py_script

    output:
    path 'oligos_counts'

    script:
    """
    python ${py_script} --fastq <(zcat ${reads}) \
        --oligos ${oligos} \
        --out oligos_counts 2> count_oligos.log

    """

}

process TRIM_FIRST_BASES {
    
    label "cutadapt"

    publishDir "${params.riboseq_process_data_outDir}/trim_first_bases", mode: 'copy'

    input:
    path reads

    output:
    path 'reads_out'

    script:
    """
    (cutadapt --cut ${params.cut} --minimum-length ${params.minimum_length} ${reads} | gzip > reads_out) 2> trim_first_bases.log

    """

}

process CLIP_READS {
    
    label "fastx"

    publishDir "${params.riboseq_process_data_outDir}/clip_reads", mode: 'copy'

    // "${sample_dict}[\$prefix]["minimum_quality"]

    input:
    path reads

    output:
    path 'reads_out'

    script:
    """
    (fastx_clipper ${params.clip_reads_v} ${params.clip_reads_n} -l ${params.clip_reads_l} ${params.clip_reads_c} ${params.clip_reads_z} -a ${params.clip_reads_adapter} -i <(zcat ${reads}) -o reads_out) 2> clip_reads.log

    """

}

/*
process TRIM_READS {
    
    label "fastx"

    publishDir "${params.riboseq_process_data_outDir}/trim_reads", mode: 'copy'

    input:
    path reads

    output:
    path 'reads_out'

    script:
    """
    fastq_quality_trimmer ${params.trim_reads_v} -l ${params.trim_reads_l} -t ${params.trim_reads_t} -Q ${params.trim_reads_Q} ${params.trim_reads_z} -i <(zcat ${reads}) -o reads_out 2> trim_reads.log

    """

}

process FILTER_READS {
    
    label "fastx"

    publishDir "${params.riboseq_process_data_outDir}/filter_reads", mode: 'copy'

    input:
    path reads

    output:
    path 'reads_out'

    script:
    """
    fastq_quality_filter ${params.filter_reads_v} -q ${params.filter_reads_q} -p ${params.filter_reads_p} -Q ${params.filter_reads_Q} ${params.filter_reads_z} -i <(zcat ${reads}) -o reads_out 2> filter_reads.log


    """

}

process FASTQ_TO_FASTA {
    
    label "fastx"

    publishDir "${params.riboseq_process_data_outDir}/filter_reads", mode: 'copy'

    input:
    path reads

    output:
    path 'reads_out'

    script:
    """
    fastq_to_fasta ${params.fq_to_fa_v} ${params.fq_to_fa_n} ${params.fq_to_fa_r} -i <(zcat ${reads}) -o reads_out 2> fastq_to_fasta.log


    """

}

process MAP_TO_OTHER_GENES {
    
    label "fastx"

    publishDir "${params.riboseq_process_data_outDir}/filter_reads", mode: 'copy'

    input:
    path reads
    path index
    path sequence

    //threads:    8

    output:
    path 'sam'
    path 'reads_out'

    script:
    """
    segemehl.x ${params.map_to_other_genes_silent} -i ${index} -d ${sequence} -q ${reads} --accuracy ${params.map_to_other_genes_accuracy} --threads ${threads} -o sam -u reads_out 2> map_to_other_genes.log

    : '
    segemehl.x \ ${params.map_to_other_genes_silent} \
        -i ${index} \
        -d ${sequence} \
        -q ${reads} \
        --accuracy ${params.map_to_other_genes_accuracy} \
        --threads {threads} \
        -o sam \
        -u reads_out 2> map_to_other_genes.log
    '
    """

}

process COUNT_OVERREPRESENTED_SEQUENCES_OTHER {
    
    label "pysam"

    publishDir "${params.riboseq_process_data_outDir}/count_overrepresented_sequences_other", mode: 'copy'

    input:
    path sam
    path script_py

    output:
    path 'overrepresented_sequences_counts'

    script:
    """
    python ${script_py} --sam ${sam} --out overrepresented_sequences_counts 2> count_overrepresented_sequences_other.log"
    
    # grep -P -v \"^@\" {input.sam} | cut -f 10 | sort | uniq -c | sort -n -r > overrepresented_sequences_counts 2> count_overrepresented_sequences_other.log

    """

}

process MAP_TO_TRANSCRIPTS {
    
    label "segemehl"

    publishDir "${params.riboseq_process_data_outDir}/map_to_transcripts", mode: 'copy'

    input:
    path reads
    path index
    path sequence

    //threads:    8

    output:
    path 'sam'
    path 'reads_out'

    script:
    """
    segemehl.x ${params.map_to_transcripts_silent} -i ${index} -d ${sequence} -q ${reads} --accuracy ${params.map_to_transcripts_accuracy} --threads ${threads} -o 
    sam -u reads_out 2> map_to_transcripts.log

    """

}

process REMOVE_MULTIMAPPERS {
    
    // basic container like alpine
    //label "segemehl"

    publishDir "${params.riboseq_process_data_outDir}/remove_multimappers", mode: 'copy'

    input:
    path sam

    // threads:     1  

    output:
    path 'sam_out'

    script:
    """
    (grep -P \"^@|\tNH:i:1\t\" ${sam} > sam_out) 2> remove_multimappers.log

    """

}

process SAM2BAM_SORT_AND_INDEX {
    
    label "samtools"

    publishDir "${params.riboseq_process_data_outDir}/sam2bam_sort_and_index", mode: 'copy'

    input:
    path sam

    // threads:     1  

    output:
    path 'bam_out'
    path 'bai' // possible error

    script:
    """
    samtools view -bS ${sam} | samtools sort - > bam_out; samtools index bam_out; 2> sam2bam_sort_and_index.log

    """

}

process READ_LENGTH_HISTOGRAM {
    
    label "rcrunch_python"

    publishDir "${params.riboseq_process_data_outDir}/read_length_histogram", mode: 'copy'

    input:
    path sam
    path script_py // scripts/plot_read_lengths.py

    // threads:     1 
    // attention: strange ouput you might have to fix

    output:
    path 'read_length_plot'

    script:
    """
    python ${script_py} --sam 
    ${sam} --outdir ${params.read_length_histogram_dir} 2> read_length_histogram.log

    """

}

process DETERMINE_P_SITE_OFFSET {
    
    label "pysam"

    publishDir "${params.riboseq_process_data_outDir}/determine_p_site_offset", mode: 'copy'

    input:
    path bam
    path transcript_id_gene_id_CDS
    path script_py // scripts/determine_p_site_offsets.py

    threads:     1  

    output:
    path 'p_site_offsets'
    path 'p_site_offset'

    script:
    """
    python ${script_py} --bam ${bam} --cds_coordinates ${transcript_id_gene_id_CDS} --outdir ${params.determine_p_site_offset_outdir} 2> determine_p_site_offset.log

    """

}

process COUNT_READS {
    
    label "pysam"

    publishDir "${params.riboseq_process_data_outDir}/count_reads", mode: 'copy'

    input:
    path bam
    path transcript_id_gene_id_CDS
    path p_site_offsets
    path script_py // scripts/count_reads.py

    threads:     1  

    output:
    path 'counts'
    path 'count'

    script:
    """
    python ${script_py} --bam ${bam} --tsv ${transcript_id_gene_id_CDS} --json ${p_site_offsets} --outdir ${params.count_reads_outdir} 2> count_reads.log

    """

}

process CHECK_PERIODICITY {
    
    label "rcrunch_python"

    publishDir "${params.riboseq_process_data_outDir}/check_periodicity", mode: 'copy'

    input:
    path bam
    path transcript_id_gene_id_CDS
    path p_site_offsets
    path script_py // scripts/check_periodicity.py

    threads:     1  

    output:
    path 'start_periodicity'
    path 'stop_periodicity'
    path 'dict_periodicity'

    script:
    """
    python ${script_py} --bam ${bam} --tsv ${transcript_id_gene_id_CDS} --json ${p_site_offsets} --outdir ${params.check_periodicity_outdir} --codnum ${params.check_peridocitiy_codnum} 2> check_periodicity.log

    """

}

process FILTER_READS_BASED_ON_READ_LENGTHS_AND_OFFSETS {
    
    label "pysam"

    publishDir "${params.riboseq_process_data_outDir}/filter_reads_based_on_read_lengths_and_offsets", mode: 'copy'

    input:
    path bam
    path p_site_offsets
    path script_py // scripts/filter_reads_based_on_read_lengths_and_offsets.py

    threads:     1  

    output:
    path 'bam_out'

    script:
    """
    python ${script_py} --bam ${bam} --p_site_offsets ${p_site_offsets} --bam_out bam_out 2> filter_reads_based_on_read_lengths_and_offsets.log

    """

}

process BAM_SORT_AND_INDEX {
    
    label "samtools"

    publishDir "${params.riboseq_process_data_outDir}/bam_sort_and_index", mode: 'copy'

    input:
    path bam

    threads:     1  

    output:
    path 'bam_out'
    path 'bai_out'

    script:
    """
    samtools sort ${bam} > bam_out; samtools index bai_out; 2> bam_sort_and_index.log

    """

}
*/


workflow RIBOSEQ_PROCESS_DATA_PIPE {

    take:
    reads_ch
    oligos_ch
    other_RNAs_sequence_ch
    count_oligos_script_ch
    other_RNAs_index_ch
    transcripts_index_ch
    find_overepresented_sequences_script_ch
    transcripts_sequence_ch
    plot_read_lengths_script_ch
    transcript_id_gene_id_CDS_ch
    determine_p_site_offsets_script_ch
    count_reads_script_ch
    check_periodicity_script_ch
    filter_reads_based_on_read_lengths_and_offsets_script_ch


    main:
      oligos_counts = COUNT_OLIGOS(reads_ch, oligos_ch, count_oligos_script_ch)
      
      trimmed_first_bases_fastq = TRIM_FIRST_BASES(reads_ch)
      
      pro_clipped_fastq = CLIP_READS(trimmed_first_bases_fastq)

      /*
      pro_trimmed_fastq = TRIM_READS(pro_clipped_fastq)

      pro_filtered_fastq = FILTER_READS(pro_trimmed_fastq)

      pro_filtered_fasta = FASTQ_TO_FASTA(pro_filtered_fastq)
      
      other_genes_mapped_sam, other_genes_unmapped_sam = MAP_TO_OTHER_GENES(pro_filtered_fasta, other_RNAs_index_ch, other_RNAs_sequence_ch)

      overepressented_sequences_other = COUNT_OVERREPRESENTED_SEQUENCES_OTHER(other_genes_mapped_sam, find_overepresented_sequences_script_ch)

      transcripts_mapped_sam, transcripts_unmapped = MAP_TO_TRANSCRIPTS(other_genes_unmapped_sam, transcripts_index_ch, transcripts_sequence_ch)

      transcripts_mapped_unique_sam = REMOVE_MULTIMAPPERS(transcripts_mapped_sam)

      transcripts_mapped_unique_sorted_bam, transcripts_mapped_unique_sorted_bam_bai = SAM2BAM_SORT_AND_INDEX(transcripts_mapped_unique_sam)

      read_length_histogram_pdf = READ_LENGTH_HISTOGRAM(transcripts_mapped_unique_sam, plot_read_lengths_script_ch)

      alignment_offset_json, p_site_offsets = DETERMINE_P_SITE_OFFSET(transcripts_mapped_unique_sorted_bam, transcript_id_gene_id_CDS_ch, determine_p_site_offsets_script_ch)

      counts_tsv, counts = COUNT_READS(transcripts_mapped_unique_sorted_bam, transcript_id_gene_id_CDS_ch, alignment_offset_json, count_reads_script_ch)

      periodicity_start_pdf, periodicity_stop_pdf, periodicity_analysis_start_ribo_seq = CHECK_PERIODICITY(transcripts_mapped_unique_sorted_bam, transcript_id_gene_id_CDS_ch, alignment_offset_json, check_periodicity_script_ch)

      transcripts_mapped_unique_a_site_profile_bam = FILTER_READS_BASED_ON_READ_LENGTHS_AND_OFFSETS(transcripts_mapped_unique_sorted_bam, alignment_offset_json, filter_reads_based_on_read_lengths_and_offsets_script_ch)

      transcripts_mapped_unique_a_site_profile_sorted_bam, transcripts_mapped_unique_a_site_profile_sorted_bam_bai = BAM_SORT_AND_INDEX(transcripts_mapped_unique_a_site_profile_bam)
      */

    emit:
      oligos_counts


}


