nextflow.enable.dsl=2

process MAKE {

	label "samtools"

	publishDir "${projectDir}/results/make_test_files", mode: 'copy', pattern: 'TEST*'
	publishDir "${projectDir}/data/test_riboseq", mode: 'copy', pattern: 'TEST*.fastq.gz'
	publishDir "${projectDir}/data/", mode: 'copy', pattern: 'TEST*.gtf'
	publishDir "${projectDir}/data/", mode: 'copy', pattern: 'TEST*.fa'
	publishDir "${projectDir}/data/", mode: 'copy', pattern: 'TEST*.fai'

	input:
	val lines
	path riboseq_file
	path out
	path gtf
	path genome

	output:
	path 'TEST*.fastq.gz', emit: test_riboseq
	path 'TEST*.gtf', emit: test_gtf
	path 'TEST*.fa', emit: test_genome
	path '*.fai', emit: test_fai

	script:
	"""
	echo \$(pwd)

	zcat ${riboseq_file} \
		| head -n ${lines} \
		| gzip -c > TEST_${riboseq_file}

	head -n ${lines} ${gtf} > TEST_${gtf}
	
	head -n ${lines} ${genome} > TEST_${genome}
	samtools faidx TEST_${genome}

	"""

}

workflow MAKE_TEST_FILES {

	take:
	n_test_lines
	riboseq_reads_ch
	out_path
	gtf
	genome

	main:
	MAKE(
		n_test_lines,
		riboseq_reads_ch,
		out_path,
		gtf,
		genome
	)

	riboseq_test = MAKE.out.test_riboseq
	gtf_test = MAKE.out.test_gtf
	genome_test = MAKE.out.test_genome
	fai_test = MAKE.out.test_fai

	emit:
	riboseq_test
	gtf_test
	genome_test
	fai_test
}


