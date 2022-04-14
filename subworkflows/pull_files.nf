nextflow.enable.dsl=2

process PULL {

	echo true
	label "philosopher"

	publishDir "${projectDir}/results/pull_files/pull", mode: 'copy', pattern: '{*.fai, *.gtf, *.gff3}'
	publishDir "${params.log_dir}/pull_files/pull", mode: 'copy', pattern: '*.log'

	input:
	path dir

	output:
	path '*.fa', emit: genome_fasta
	path '*.gtf', emit: gtf
	path '*.gff3', emit: gff3
	path '*.log', emit: log

	script:
	"""
	echo \$(pwd)
	# whole genome DNA fasta file
	wget http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
		&> genome_fasta.log
	# gff3 annotation files of noncoding RNA and mRNA
	wget http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3/mus_musculus.GRCm39.gff3.gz \
		&> ncRNA_annotation.log
	wget http://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.105.chr.gtf.gz \
		&> mRNA_annotation.log

	# decompress files
	zcat Mus_musculus.GRCm39.dna.primary_assembly.fa.gz > Mus_musculus.GRCm39.dna.primary_assembly.fa
	zcat mus_musculus.GRCm39.gff3.gz > mus_musculus.GRCm39.gff3
	zcat Mus_musculus.GRCm39.105.chr.gtf.gz > Mus_musculus.GRCm39.105.chr.gtf

	cp *.fa *.gtf *.gff3 ${dir}
	"""


}

process FASTA_INDEX {

	label "samtools"

	publishDir "${projectDir}/results/pull_files/fasta_index", mode: 'copy', pattern: '*.fai'
	publishDir "${params.log_dir}/fasta_index", mode: 'copy', pattern: '*.log'

	input:
	path genome_fasta
	path dir

	output:
	path "*.fai", emit: fai
	path "*.log", emit: log

	script:
	"""
	samtools faidx \
		${genome_fasta} \
		&> fasta_index.log

	cp *.fai ${dir}
	"""

}

process GFFREAD {

	label "gffread"

	publishDir "${projectDir}/results/pull_files/gffread", mode: 'copy', pattern: ''
	publishDir "${params.log_dir}/gffread", mode: 'copy', pattern: '*.log'

	input:
	path gff3

	output:
	path "", emit: gtf
	path "*.log", emit: log

	script:
	"""
	input=\$(basename ${gff3})
    prefix=\$(echo \$input | cut -d '.' -f 1,2)

	gffread ${gff3} \
		-T \
		-o \${prefix}.gtf \
		&> gff3_to_gtf.log

	"""
}

process CONCAT {

	label "gffread"

	publishDir "${projectDir}/results/pull_files/concat", mode: 'copy', pattern: ''
    publishDir "${params.log_dir}/concat", mode: 'copy', pattern: '*.log' 

	input:
	path gene_gtf
	path noncoding_gtf

	output:
	path "combined.gtf", emit: gtf
	path "*.log", emit: log, optional: true

	script:
	"""
	cat ${gene_gtf} ${noncoding_gtf} >> combined.gtf
	"""
}

process MOVE {

	label "philosopher"

	publishDir "${projectDir}/results/pull_files/move", mode: 'copy', pattern: ''
	publishDir "${params.log_dir}/move", mode: 'copy', pattern: ''

	input:
	path genome_fasta
	path genome_fai
	path gene_gtf
	path noncoding_gtf
	path combined_gtf

	output:

	script:
	"""
	cp ${genome_fasta} ${params.genome}
	cp ${genome_fai} ${params.genome_fai}
	cp ${gene_gtf} ${params.gtf}
	cp ${noncoding_gtf} ${params.rnacentral_gtf}
	cp ${combined_gtf} ${params.combined_gtf}
	"""

}

workflow PULLFILES {

	take:
	directory_path

	main:
	// this subworkflow pulls gtf, gff3 and fasta files from the web
	// and converts gff3 files to gtf, then concatenates the gtf files
	// into 1 larger gtf file which will be used to map reads to.
	// the workflow also puts the files in the right dir
	// and writes this to the nextflow.config file

	PULL(
		directory_path
	)
	PULL.out.gff3.view()
	FASTA_INDEX(
		PULL.out.genome_fasta,
		directory_path
	)

	GFFREAD(
		PULL.out.gff3
	)

	CONCAT(
		PULL.out.gtf,
		GFFREAD.out.gtf
	)
	combined_gtf = CONCAT.out.gtf

	MOVE(
		PULL.out.genome_fasta,
        FASTA_INDEX.out.fai,
        PULL.out.gtf,
        GFFREAD.out.gtf,
		combined_gtf
	)

	emit:
	combined_gtf

}




