#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process SELECT_LONGEST_CODING_TRANSCRIPT {
    
    label "htseq"

    publishDir "${params.riboseq_annotate_outDir}/select_longest_coding_transcript", mode: 'copy'

    input:
    path input_gtf
    path select_longest_ct_py

    output:
    path 'longest_coding_transcript_per_gene.gtf'

    script:
    """
    python ${select_longest_ct_py} \
	--gtf ${input_gtf} \
	--out longest_coding_transcript_per_gene.gtf \
	2> select_longest_coding_transcript.log

    """

}

process GENERATE_SEGEMEHL_INDEX_OTHER_RNAS {
    
    label "segemehl"

    publishDir "${params.riboseq_annotate_outDir}/generate_segemehl_index_other_rnas", mode: 'copy'

    input:
    path sequence

    output:
    path 'segemehl_index_other_rnas.idx'

    script:
    """
    segemehl.x \
	-x segemehl_index_other_rnas.idx \
	-d ${sequence} \
	2> generate_segemehl_index_other_rnas.log

    """

}

process GENERATE_STAR_INDEX_OTHER_RNAS {

    label 'star'

    publishDir "${params.riboseq_annotate_outDir}/generate_star_index_other_rnas", mode: 'copy'

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

process EXTRACT_TRANSCRIPT_SEQUENCES {
    
    label "cufflinks"

    publishDir "${params.riboseq_annotate_outDir}/extract_transcript_sequences", mode: 'copy'

    input:
    path gtf
    path genome

    output:
    path 'transcripts_sequences.out'

    script:
    """
    gffread \
	${gtf} \
	-g ${genome} \
	-w transcripts_sequences.out \
	2> extract_transcript_sequences.log

    """

}

process CREATE_TAB_DELIMITED_CDS_FILE {
    
    echo true

    label "htseq_biopython"

    publishDir "${params.riboseq_annotate_outDir}/create_tab_delimited_CDS_file", mode: 'copy'

    input:
    path gtf
    path transcripts
    path td_CDS_script_py

    output:
    path 'CDS.tsv'

    script:
    """
    python ${td_CDS_script_py} \
	--gtf ${gtf} \
	--fasta ${transcripts} \
	--out CDS.tsv \
	2> create_tab_delimited_CDS_file.log

    """

}

process CREATE_BED_CDS_FILE {
    
    label "htseq_biopython"

    publishDir "${params.riboseq_annotate_outDir}/create_bed_cds_file", mode: 'copy'

    input:
    path tsv

    output:
    path 'CDS.bed'

    script:
    """
    tail -n+2 ${tsv} \
	| awk '{print \$1 "\t" \$3-1 "\t" \$4 "\t" \$2 }' > CDS.bed \
	2> create_bed_cds_file.log

    """

}

process GENERATE_SEGEMEHL_INDEX_TRANSCRIPTS {
    
    label "segemehl"

    publishDir "${params.riboseq_annotate_outDir}/generate_segemehl_index_transcripts", mode: 'copy'

    input:
    path sequence

    output:
    path 'segemehl_index_transcripts.idx'

    script:
    """
    segemehl.x \
	-x segemehl_index_transcripts.idx \
	-d ${sequence} \
	2> generate_segemehl_index_transcripts.log

    """

}

process GENERATE_STAR_INDEX_TRANSCRIPTS {

    label 'star'

    publishDir "${params.riboseq_annotate_outDir}/generate_star_index_transcripts", mode: 'copy'

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

workflow RIBOSEQ_ANNOTATE_PIPE {

    take:
      pull_containers_ch
      gtf_ch
      lct_script_ch
      other_RNAs_sequence_ch
      genome_ch
      ctdCDS_script_ch 

    main:
      longest_pc_transcript_per_gene_gtf = SELECT_LONGEST_CODING_TRANSCRIPT(gtf_ch, lct_script_ch)

      if ( params.other_genes_index_mode == "segemehl" ) {
	  other_RNAs_sequence_idx = GENERATE_SEGEMEHL_INDEX_OTHER_RNAS(other_RNAs_sequence_ch)
      }

      if ( params.other_genes_index_mode == "star" ) {
	  other_RNAs_sequence_idx = GENERATE_STAR_INDEX_OTHER_RNAS(other_RNAs_sequence_ch)
      }

      longest_pc_transcript_per_gene_fa = EXTRACT_TRANSCRIPT_SEQUENCES(  longest_pc_transcript_per_gene_gtf,
									 genome_ch  )

      transcript_id_gene_id_CDS_tsv = CREATE_TAB_DELIMITED_CDS_FILE(  longest_pc_transcript_per_gene_gtf,
								      longest_pc_transcript_per_gene_fa,
								      ctdCDS_script_ch  )

      transcript_id_gene_id_CDS_bed = CREATE_BED_CDS_FILE(transcript_id_gene_id_CDS_tsv)

      if ( params.genome_index_mode == "segemehl" ) {
          genome_idx = GENERATE_SEGEMEHL_INDEX_TRANSCRIPTS(longest_pc_transcript_per_gene_fa)
      }

      if ( params.genome_index_mode == "star" ) {
       	  genome_idx = GENERATE_STAR_INDEX_TRANSCRIPTS(genome_ch)
      }

    emit:
      longest_pc_transcript_per_gene_gtf
      other_RNAs_sequence_idx
      longest_pc_transcript_per_gene_fa
      transcript_id_gene_id_CDS_tsv
      transcript_id_gene_id_CDS_bed
      genome_idx

}













