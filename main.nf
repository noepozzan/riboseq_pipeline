nextflow.enable.dsl=2

include { PHILOSOPHER_PIPE } from './subworkflows/simple_data_analysis.nf'
include { RIBOSEQ_ANNOTATE_PIPE } from './subworkflows/riboseq_annotate.nf'
include { RIBOSEQ_PROCESS_DATA_PIPE } from './subworkflows/riboseq_process_data.nf'
include { PULL_CONTAINERS } from './subworkflows/pull_containers.nf'
include { RIBOTISH_PIPE } from './subworkflows/ribotish_connect.nf'
include { CHECK_FILES_PIPE } from './subworkflows/check_files.nf'

// pull containers channels
config_file_ch = Channel.fromPath(params.slurm_config)
pull_file_ch = Channel.fromPath(params.pull_containers_file)

// check_files channels
check_files_script_ch = Channel.fromPath(params.check_files_script)

// riboseq annotate channels
gtf_ch = Channel.fromPath(params.gtf)
lct_script_ch = Channel.fromPath(params.lct_script)
other_RNAs_sequence_ch = Channel.fromPath(params.other_RNAs_sequence)
genome_ch = Channel.fromPath(params.genome)
ctdCDS_script_ch = Channel.fromPath(params.ctdCDS_script)

// riboseq process data channels
riboseq_reads_ch = Channel.fromPath(params.riboseq_reads)
oligos_ch = Channel.fromPath(params.oligos)
other_RNAs_sequence_ch = Channel.fromPath(params.other_RNAs_sequence)
count_oligos_script_ch = Channel.fromPath(params.count_oligos_script)
// other_RNAs_index_ch = RIBOSEQ_ANNOTATE_PIPE.out.other_RNAs_sequence_idx
// transcripts_index_ch = RIBOSEQ_ANNOTATE_PIPE.outlongest_pc_transcript_per_gene_idx
find_overrepresented_sequences_script_ch = Channel.fromPath(params.find_overrepresented_sequences_script)
// transcripts_sequence_ch = RIBOSEQ_ANNOTATE_PIPE.out.longest_pc_transcript_per_gene_fa
plot_read_lengths_script_ch = Channel.fromPath(params.plot_read_lengths_script)
// transcript_id_gene_id_CDS_ch = RIBOSEQ_ANNOTATE_PIPE.out.transcript_id_gene_id_CDS_tsv
determine_p_site_offsets_script_ch = Channel.fromPath(params.determine_p_site_offsets_script)
count_reads_script_ch = Channel.fromPath(params.count_reads_script)
check_periodicity_script_ch = Channel.fromPath(params.check_periodicity_script)
filter_reads_based_on_read_lengths_and_offsets_script_ch = Channel.fromPath(params.filter_reads_based_on_read_lengths_and_offsets_script)

// ribotish channels
genome_fai_ch = Channel.fromPath(params.genome_fai)

// philosopher channels
proteomics_reads_ch = Channel.fromPath(params.proteomics_reads)
philosopher_db_ch = Channel.fromPath(params.philosopher_db)
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

    PULL_CONTAINERS(  config_file_ch,
		      pull_file_ch  )

    CHECK_FILES_PIPE(  riboseq_reads_ch,
		       proteomics_reads_ch,
		       genome_ch,
		       genome_fai_ch,
		       gtf_ch,
		       other_RNAs_sequence_ch,
		       check_files_script_ch  )

    RIBOSEQ_ANNOTATE_PIPE(  gtf_ch,
                            lct_script_ch,
                            other_RNAs_sequence_ch,
                            genome_ch,
                            ctdCDS_script_ch  )

    RIBOSEQ_PROCESS_DATA_PIPE(  riboseq_reads_ch,
                                oligos_ch,
                                other_RNAs_sequence_ch,
                                count_oligos_script_ch,
                                RIBOSEQ_ANNOTATE_PIPE.out.other_RNAs_sequence_idx,
                                RIBOSEQ_ANNOTATE_PIPE.out.genome_idx,
                                find_overrepresented_sequences_script_ch,
                                RIBOSEQ_ANNOTATE_PIPE.out.longest_pc_transcript_per_gene_fa,
                                plot_read_lengths_script_ch,
                                RIBOSEQ_ANNOTATE_PIPE.out.transcript_id_gene_id_CDS_tsv,
                                determine_p_site_offsets_script_ch,
                                count_reads_script_ch,
                                check_periodicity_script_ch,
                                filter_reads_based_on_read_lengths_and_offsets_script_ch,
								gtf_ch  )
    
    RIBOTISH_PIPE(  gtf_ch,
					RIBOSEQ_PROCESS_DATA_PIPE.out.bam_sort_index_folder,
					genome_ch  )


    PHILOSOPHER_PIPE(  RIBOTISH_PIPE.out.ribo_pred,
					   proteomics_reads_ch,
					   philosopher_db_ch,
					   change_file_script_ch  )      
  
  

}
