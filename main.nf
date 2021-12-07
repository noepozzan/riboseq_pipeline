nextflow.enable.dsl=2

include { PHILOSOPHER_PIPE } from './workspace/simple_data_analysis.nf'
include { RIBOSEQ_ANNOTATE_PIPE } from './riboseq/riboseq_annotate.nf'
include { RIBOSEQ_PROCESS_DATA_PIPE } from './riboseq/riboseq_process_data.nf'

// riboseq annotate channels
gtf_ch = Channel.fromPath(params.gtf)
lct_script_ch = Channel.fromPath(params.lct_script)
other_RNAs_sequence_ch = Channel.fromPath(params.other_RNAs_sequence)
genome_ch = Channel.fromPath(params.genome)
ctdCDS_script_ch = Channel.fromPath(params.ctdCDS_script)

// riboseq process data channels
reads_ch = Channel.fromPath(params.reads)
oligos_ch = Channel.fromPath(params.oligos)
other_RNAs_sequence_ch = Channel.fromPath(params.other_RNAs_sequence)
count_oligos_script_ch = Channel.fromPath(params.count_oligos_script)
// other_RNAs_index_ch = RIBOSEQ_ANNOTATE_PIPE.out.other_RNAs_sequence_idx
// transcripts_index_ch = RIBOSEQ_ANNOTATE_PIPE.outlongest_pc_transcript_per_gene_idx
find_overepresented_sequences_script_ch = Channel.fromPath(params.find_overepresented_sequences_script)
// transcripts_sequence_ch = RIBOSEQ_ANNOTATE_PIPE.out.longest_pc_transcript_per_gene_fa
plot_read_lengths_script_ch = Channel.fromPath(params.plot_read_lengths_script)
// transcript_id_gene_id_CDS_ch = RIBOSEQ_ANNOTATE_PIPE.out.transcript_id_gene_id_CDS_tsv
determine_p_site_offsets_script_ch = Channel.fromPath(params.determine_p_site_offsets_script)
count_reads_script_ch = Channel.fromPath(params.count_reads_script)
check_periodicity_script_ch = Channel.fromPath(params.check_periodicity_script)
filter_reads_based_on_read_lengths_and_offsets_script_ch = Channel.fromPath(params.filter_reads_based_on_read_lengths_and_offsets_script)


// philosopher channels
input_ch = Channel.fromPath(params.philosopher_input)
input_ch.view()
db_ch = Channel.fromPath(params.philosopher_db)
change_file_script_ch = Channel.fromPath(params.change_file_script)

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf -profile docker

        Mandatory arguments:
         --something                     something
         --something else                something else

        Optional arguments:
         --well                          still some else
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


workflow {

    RIBOSEQ_ANNOTATE_PIPE(  gtf_ch,
                            lct_script_ch,
                            other_RNAs_sequence_ch,
                            genome_ch,
                            ctdCDS_script_ch  )

    RIBOSEQ_PROCESS_DATA_PIPE(  reads_ch,
                                oligos_ch,
                                other_RNAs_sequence_ch,
                                count_oligos_script_ch,
                                RIBOSEQ_ANNOTATE_PIPE.out.other_RNAs_sequence_idx, // other_RNAs_index_ch
                                RIBOSEQ_ANNOTATE_PIPE.out.longest_pc_transcript_per_gene_idx, // transcripts_index_ch
                                find_overepresented_sequences_script_ch,
                                RIBOSEQ_ANNOTATE_PIPE.out.longest_pc_transcript_per_gene_fa, // transcripts_sequence_ch
                                plot_read_lengths_script_ch,
                                RIBOSEQ_ANNOTATE_PIPE.out.transcript_id_gene_id_CDS_tsv, // transcript_id_gene_id_CDS_ch
                                determine_p_site_offsets_script_ch,
                                count_reads_script_ch,
                                check_periodicity_script_ch,
                                filter_reads_based_on_read_lengths_and_offsets_script_ch  )

    //PHILOSOPHER_PIPE(RIBOSEQ_PROCESS_DATA_PIPE.out[1], input_ch, db_ch)      

}
