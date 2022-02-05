nextflow.enable.dsl=2

include { PHILOSOPHER_PIPE } from './subworkflows/philosopher.nf'
include { ANNOTATE_PIPE } from './subworkflows/annotate.nf'
include { RIBOSEQ_PIPE } from './subworkflows/riboseq.nf'
include { PULL_CONTAINERS } from './subworkflows/pull_containers.nf'
include { RIBOTISH_PIPE } from './subworkflows/ribotish.nf'
include { CHECK_FILES_PIPE } from './subworkflows/check_files.nf'

// pull containers channels
config_file_ch = Channel.fromPath(params.slurm_config)
pull_file_ch = Channel.fromPath(params.pull_containers_file)

// check_files channels
check_files_script_ch = Channel.fromPath(params.check_files_script)

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

gtf_ch = channel.fromPath('./data/Mus_musculus.GRCm38.99.chr.gtf')
other_RNAs_sequence_ch = channel.fromPath('./data/mm10_rrnas.fa')
genome_ch = channel.fromPath('./data/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa')
genome_fai_ch = channel.fromPath('./data/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.fai')
riboseq_reads_ch = channel.fromPath('./data/riboseq_reads/*.gz')
oligos_ch = channel.fromPath('./data/oligos.txt')
proteomics_reads_ch = channel.fromPath('./data/proteomics_reads/*.mzML')
//only temporary, to be deleted as soon as ribotish outputs fasta file of predicted sORFs
philosopher_db_ch = channel.fromPath('./data/small_peptides_all_quad_samples.fasta')

workflow {

/*
    PULL_CONTAINERS(
		config_file_ch,
		pull_file_ch
	)

    CHECK_FILES_PIPE(
		riboseq_reads_ch,
		proteomics_reads_ch,
		genome_ch,
		genome_fai_ch,
		gtf_ch,
		other_RNAs_sequence_ch,
		check_files_script_ch
	)
*/

    ANNOTATE_PIPE(  
		gtf_ch,
		other_RNAs_sequence_ch,
		genome_ch,
	)

    RIBOSEQ_PIPE(  
		genome_ch,
		riboseq_reads_ch,
		oligos_ch,
		other_RNAs_sequence_ch,
		ANNOTATE_PIPE.out.longest_pc_transcript_per_gene_fa,
		ANNOTATE_PIPE.out.transcript_id_gene_id_CDS_tsv,
		gtf_ch
	)

    RIBOTISH_PIPE(
		gtf_ch,
		RIBOSEQ_PIPE.out.bam_sort_index_folder,
		genome_ch,
		genome_fai_ch,
		
	)


    PHILOSOPHER_PIPE(  
		RIBOTISH_PIPE.out.speptide,
		proteomics_reads_ch,
		philosopher_db_ch,
	)

}
