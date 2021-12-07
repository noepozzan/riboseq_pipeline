#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process SELECT_LONGEST_CODING_TRANSCRIPT {
    
    label "htseq"

    publishDir "${params.riboseq_annotate_outDir}/select_longest_coding_transcript", mode: 'copy'

    input:
    path input_gtf
    path select_longest_ct_py

    output:
    path 'gtf_out'

    script:
    """
    python ${select_longest_ct_py} --gtf ${input_gtf} --out gtf_out 2> select_longest_coding_transcript.log

    """

}

process GENERATE_SEGEMEHL_INDEX_OTHER_RNAS {
    
    label "segemehl"

    publishDir "${params.riboseq_annotate_outDir}/generate_segemehl_index_other_rnas", mode: 'copy'

    input:
    path sequence

    output:
    path 'idx_out'

    script:
    """
    segemehl.x -x idx_out -d ${sequence} 2> generate_segemehl_index_other_rnas.log

    """

}

process EXTRACT_TRANSCRIPT_SEQUENCES {
    
    label "cufflinks"

    publishDir "${params.riboseq_annotate_outDir}/extract_transcript_sequences", mode: 'copy'

    input:
    path gtf
    path genome

    output:
    path 'transcripts_out'

    script:
    """
    gffread ${gtf} -g ${genome} -w transcripts_out 2> extract_transcript_sequences.log

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
    path 'tsv_out'

    script:
    """
    python ${td_CDS_script_py} --gtf ${gtf} --fasta ${transcripts} --out tsv_out 2> create_tab_delimited_CDS_file.log

    """

}

process CREATE_BED_CDS_FILE {
    
    label "htseq_biopython"

    publishDir "${params.riboseq_annotate_outDir}/create_bed_cds_file", mode: 'copy'

    input:
    path tsv

    output:
    path 'bed'

    script:
    """
    (tail -n+2 ${tsv} | awk '{print \$1 "\t" \$3-1 "\t" \$4 "\t" \$2 }' > bed) 2> create_bed_cds_file.log

    """

}

process GENERATE_SEGEMEHL_INDEX_TRANSCRIPTS {
    
    label "segemehl"

    publishDir "${params.riboseq_annotate_outDir}/generate_segemehl_index_transcripts", mode: 'copy'

    input:
    path sequence

    output:
    path 'idx'

    script:
    """
    segemehl.x -x idx -d ${sequence} 2> generate_segemehl_index_transcripts.log

    """

}


workflow RIBOSEQ_ANNOTATE_PIPE {

    take:
      gtf_ch
      lct_script_ch
      other_RNAs_sequence_ch
      genome_ch
      ctdCDS_script_ch 

    main:
      longest_pc_transcript_per_gene_gtf = SELECT_LONGEST_CODING_TRANSCRIPT(gtf_ch, lct_script_ch)
      
      other_RNAs_sequence_idx = GENERATE_SEGEMEHL_INDEX_OTHER_RNAS(other_RNAs_sequence_ch)
      
      longest_pc_transcript_per_gene_fa = EXTRACT_TRANSCRIPT_SEQUENCES(longest_pc_transcript_per_gene_gtf, genome_ch)
      
      transcript_id_gene_id_CDS_tsv = CREATE_TAB_DELIMITED_CDS_FILE(longest_pc_transcript_per_gene_gtf, longest_pc_transcript_per_gene_fa, ctdCDS_script_ch)

      transcript_id_gene_id_CDS_bed = CREATE_BED_CDS_FILE(transcript_id_gene_id_CDS_tsv)
      
      longest_pc_transcript_per_gene_idx = GENERATE_SEGEMEHL_INDEX_TRANSCRIPTS(longest_pc_transcript_per_gene_fa)

    emit:
      longest_pc_transcript_per_gene_gtf
      other_RNAs_sequence_idx
      longest_pc_transcript_per_gene_fa
      transcript_id_gene_id_CDS_tsv
      transcript_id_gene_id_CDS_bed
      longest_pc_transcript_per_gene_idx

}

