nextflow.enable.dsl=2

include { PHILOSOPHER_PIPE } from './subworkflows/philosopher.nf'
include { ANNOTATE_PIPE } from './subworkflows/annotate.nf'
include { RIBOSEQ_PIPE } from './subworkflows/riboseq.nf'
include { PULL_CONTAINERS } from './subworkflows/pull_containers.nf'
include { RIBOTISH_PIPE } from './subworkflows/ribotish.nf'
include { CHECK_FILES_PIPE } from './subworkflows/check_files.nf'
include { QC_PIPE } from './subworkflows/qc.nf'
include { READS_PIPE } from './subworkflows/prepare_reads.nf'
include { GENOME_PIPE } from './subworkflows/map_genome.nf'
include { rRNA_PIPE } from './subworkflows/map_rrna.nf'
include { TRANSCRIPTOME_PIPE } from './subworkflows/map_transcriptome.nf'
include { QC_ONLY_PIPE } from './subworkflows/qc_only.nf'

// pull containers channels
slurm_online_config = Channel.fromPath(params.slurm_config)
//slurm_offline_config = Channel.fromPath("${projectDir}/conf/slurm_offline.config")
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

gtf_ch = channel.fromPath(params.gtf)
other_RNAs_sequence_ch = channel.fromPath(params.other_RNAs_sequence)
genome_ch = channel.fromPath(params.genome)
genome_fai_ch = channel.fromPath(params.genome_fai)
oligos_ch = channel.fromPath('./data/oligos.txt')
riboseq_reads_ch = channel.fromPath(params.riboseq_reads)
proteomics_reads_ch = channel.fromPath(params.proteomics_reads)

workflow PULLING {

	PULL_CONTAINERS(
        slurm_online_config,
        pull_file_ch
    )

}

workflow PHILOSOPHER {

	speptides = Channel.fromPath("${projectDir}/data/small_peptides_all_quad_samples.fasta")

	PHILOSOPHER_PIPE(
        speptides,
        proteomics_reads_ch
    )

}

workflow CHECK_FILES {


    CHECK_FILES_PIPE(
		riboseq_reads_ch,
		proteomics_reads_ch,
		genome_ch,
		genome_fai_ch,
		gtf_ch,
		other_RNAs_sequence_ch,
		check_files_script_ch
	)

}

workflow RIBOTISH {

	verga = Channel.fromPath("${projectDir}/data/verga/*", type: 'dir')
	RIBOTISH_PIPE(
        gtf_ch,
        verga,
        genome_ch,
        genome_fai_ch,

    )

}

/*
workflow {

	if ( params.run_qc == true ){
    ANNOTATE_PIPE(  
		gtf_ch,
		other_RNAs_sequence_ch,
		genome_ch,
	)

	QC_PIPE(
		genome_ch,
        riboseq_reads_ch,
        oligos_ch,
        other_RNAs_sequence_ch,
        ANNOTATE_PIPE.out.longest_pc_transcript_per_gene_fa,
        ANNOTATE_PIPE.out.transcript_id_gene_id_CDS_tsv,
        gtf_ch
    )  
	}

    RIBOSEQ_PIPE(  
		genome_ch,
		riboseq_reads_ch,
		other_RNAs_sequence_ch,
		gtf_ch
	)

    RIBOTISH_PIPE(
		gtf_ch,
		RIBOSEQ_PIPE.out.bam_sort_index_folder,
		genome_ch,
		genome_fai_ch,
		
	)

    PHILOSOPHER_PIPE(  
		RIBOTISH_PIPE.out.speptide_combined,
		proteomics_reads_ch
	)

}

*/

workflow {

	READS_PIPE(
		riboseq_reads_ch
	)
	fasta_reads = READS_PIPE.out

	rRNA_PIPE(
		genome_ch,
		other_RNAs_sequence_ch,
		fasta_reads
	)
	rrna_unmapped = rRNA_PIPE.out.other_genes_unmapped_fasta
	rrna_mapped_sam = rRNA_PIPE.out.other_genes_mapped_sam

	GENOME_PIPE(
		genome_ch,
		gtf_ch,
		rrna_unmapped
	)
	genome_bam_folder = GENOME_PIPE.out

	if ( params.run_qc == true ){
	ANNOTATE_PIPE(
        gtf_ch,
        other_RNAs_sequence_ch,
        genome_ch,
    )

	TRANSCRIPTOME_PIPE(
		ANNOTATE_PIPE.out.longest_pc_transcript_per_gene_fa,
		rrna_unmapped
	)
	transcripts_mapped_unique_sam = TRANSCRIPTOME_PIPE.out.transcripts_mapped_unique_sam
	transcriptome_bam_folder = TRANSCRIPTOME_PIPE.out.bam_bai_folder
	
	QC_ONLY_PIPE(
		oligos_ch,
		ANNOTATE_PIPE.out.transcript_id_gene_id_CDS_tsv,
		transcripts_mapped_unique_sam,
		transcriptome_bam_folder,
		rrna_mapped_sam
	)
	}

	RIBOTISH_PIPE(
        gtf_ch,
        genome_bam_folder,
        genome_ch,
        genome_fai_ch,

    )

    PHILOSOPHER_PIPE(
        RIBOTISH_PIPE.out.speptide_combined,
        proteomics_reads_ch
    )
	
}



