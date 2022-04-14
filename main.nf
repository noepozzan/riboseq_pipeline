nextflow.enable.dsl=2

include { PHILOSOPHER } from './subworkflows/philosopher.nf'
include { PHILOSOPHER_PARALLEL } from './subworkflows/philosopher_parallel.nf'
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
include { FIX_NAMES } from './subworkflows/fix_names.nf'
include { PULLFILES } from './subworkflows/pull_files.nf'
include { RAWMZML } from './subworkflows/raw_to_mzml.nf'

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

workflow PULL_CONTAINERS {

	slurm_online_config = Channel.fromPath(params.slurm_config)
	//slurm_offline_config = Channel.fromPath("${projectDir}/conf/slurm_offline.config")
	pull_file_ch = Channel.fromPath(params.pull_containers_file)

	PULL_CONTAINERS(
        slurm_online_config,
        pull_file_ch
    )

}

workflow RAW_TO_MZML {

	RAWMZML()

}

workflow PULL_FILES {

	directory_path = channel.fromPath("${projectDir}/data")
	PULLFILES(directory_path)

}

workflow MAP_NAMES {

	csv_file = channel.fromPath("${projectDir}/data/experimental_conditions.csv")
	fix_script = channel.fromPath("${projectDir}/data/python_scripts/fix_names.py")

	FIX_NAMES(
		csv_file,
		fix_script
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

// main workflow that calls all processes in the subworkflows dir
workflow {

	READS_PIPE(
		riboseq_reads_ch
	)

	rRNA_PIPE(
		genome_ch,
		other_RNAs_sequence_ch,
		READS_PIPE.out
	)
	
	GENOME_PIPE(
		genome_ch,
		gtf_ch,
		rRNA_PIPE.out.other_genes_unmapped_fasta
	)

	if ( params.run_qc == true ){
	ANNOTATE_PIPE(
        gtf_ch,
        other_RNAs_sequence_ch,
        genome_ch,
    )

	TRANSCRIPTOME_PIPE(
		ANNOTATE_PIPE.out.longest_pc_transcript_per_gene_fa,
		rRNA_PIPE.out.other_genes_unmapped_fasta
	)
	
	QC_ONLY_PIPE(
		oligos_ch,
		ANNOTATE_PIPE.out.transcript_id_gene_id_CDS_tsv,
		TRANSCRIPTOME_PIPE.out.transcripts_mapped_unique_sam,
		TRANSCRIPTOME_PIPE.out.bam_bai_folder,
		rRNA_PIPE.out.other_genes_mapped_sam
	)
	}

	RIBOTISH_PIPE(
        gtf_ch,
        GENOME_PIPE.out,
        genome_ch,
        genome_fai_ch,

    )

    PHILOSOPHER(
        RIBOTISH_PIPE.out.speptide_combined,
        proteomics_reads_ch
    )
/*
	PHILOSOPHER_PARALLEL(
		//PHILOSOPHER.out.ionquant,
		PHILOSOPHER.out.report,
		PHILOSOPHER.out.msfragger_params,
		RIBOTISH_PIPE.out.speptide_combined,
        proteomics_reads_ch		
	)
*/
}

workflow MOCK {

	speptide_fasta = channel.fromPath("${projectDir}/data/small_peptides_all_quad_samples.fasta")
	PHILOSOPHER(
        //RIBOTISH_PIPE.out.speptide_combined,
        speptide_fasta,
		proteomics_reads_ch
    )
	PHILOSOPHER_PARALLEL(
        PHILOSOPHER.out.report,
        PHILOSOPHER.out.msfragger_params,
        //RIBOTISH_PIPE.out.speptide_combined,
        speptide_fasta,
		proteomics_reads_ch
    )

}


