#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//https://github.com/Nesvilab/philosopher/wiki/Simple-Data-Analysis


process WORKSPACE {

    label 'philosopher'

    input:
    path ribotish_speptide

	output:
	val 'tubel'

    script:
    """
	# create new working directory (remove old if applicable,
	# to be sure that there is a fresh new workspace)
    rm -rf ${params.workspace}
    mkdir ${params.workspace}

	# move input files to working directory

	for VAR in ${params.proteomics_reads}
	do
    	mv \$VAR ${params.workspace}/
	done

	# move the predicted sPeptides file to the workspace
    mv ${ribotish_speptide} ${params.workspace}/

    # initialize philosopher in the philosopher
    # working directory which is called workspace
    cd ${params.workspace}
    philosopher workspace --clean
    philosopher workspace --init

	# copy input files back to data folder as a backup
    dir=${params.proteomics_reads}
	parentdir="\$(dirname "\$dir")"
	cp *.mzML \$parentdir
	
	"""

}

process DATABASE {

    label 'philosopher'
    
    publishDir "${params.philosopher_dir}/database", mode: 'copy'

    input:
    val workspace
    path db

    output:
    path '*.fas'

    script:
    """
    workd=\$(pwd)

	# same pattern always:
	# change into working directory, execute command,
	# then copy the generated output files back to
	# nextflow directory to allow nextflow to track files
    cd ${params.workspace}
    
	philosopher database \
		--custom ${db} \
		--contam
    
    cp ${params.workspace}/*.fas \$workd

    """

}

process GENERATE_CHANGE_PARAMS {

    label 'msfragger'

    publishDir "${params.philosopher_dir}/generate_change_params", mode: 'copy'

    input:
    path db
    path change_file_script

    output:
    path 'msfragger.params'

    script:
    """
	# generate MSFRAGGER parameter file
    java -jar /MSFragger.jar --config closed

	# python script to change some parameters
	python ${change_file_script} \
		--database ${db} \
		--param closed_fragger.params \
		--out msfragger.params \
		&> generate_change_params.log

	# copy params file to the working directory
	# since this process took place in some subdirectory
    cp msfragger.params ${params.workspace}

    """

}

process MSFRAGGER {

    label "msfragger"

    publishDir "${params.philosopher_dir}/msfragger/", mode: 'copy'

    input:
    path mzML_file
    path closed_fragger
    path db_file

    output:
    path '*.pepXML'
    
    script:
	"""
	workd=\$(pwd)	

	# change to working directory
	cd ${params.workspace}

	# search for matches of the predicted peptides
	# in the peptiomics data
	for VAR in ${mzML_file}
	do

	java -Xmx${params.fragger_ram}g \
		-jar /MSFragger.jar \
		${closed_fragger} \
		\${VAR} 

	done

	cp *.pepXML \$workd

	"""

}

process PEPTIDEPROPHET {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/peptideprophet", mode: 'copy'
       
    input:
    path db
    path pepXML
    
    output:
	path "interact*.pep.xml"
    
    script:
	"""
	workd=\$(pwd)

	cd ${params.workspace}

	# peptide assignment validation
	philosopher peptideprophet \
		${params.peptideprophet_args} \
		--database ${db} \
		${pepXML}

	cp ${params.workspace}/interact*.pep.xml \$workd

	"""

}

process PROTEINPROPHET {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/proteinprophet", mode: 'copy'
       
    input:
    path pepxml
    
    output:
    path "interact.prot.xml"
     
    script:
    if( params.skip_proteinprophet == true )
        """
		workd=\$(pwd)

		cd ${params.workspace}

        touch interact.prot.xml

		cp ${params.workspace}/interact.prot.xml \$workd
        """

    else if( params.skip_proteinprophet == false )
        """
        workd=\$(pwd)

        cd ${params.workspace}

		# group peptides by their corresponding protein(s)
		# to compute probabilities that those proteins were
		# present in in the original sample
		philosopher proteinprophet \
			${pepxml} \
			2> proteinprophet.out

        cp ${params.workspace}/interact.prot.xml \$workd

        """

}

process FILTER_FDR {
    
    label "philosopher"

	echo true

    publishDir "${params.philosopher_dir}/filter_fdr", mode: 'copy'

    input:
    path pepxml
    path protxml

    output:
    val "filtering_done_pseudo"
 
    script:
    if( params.skip_proteinprophet == true )
        """
        cd ${params.workspace}

        philosopher filter \
			${params.philosopher_filter_args} \
			--pepxml ${pepxml}

        """
	else if( params.skip_proteinprophet == false && params.combine_peptideprophet == true )
		"""
		cd ${params.workspace}

		philosopher filter \
			${params.philosopher_filter_args} \
			--pepxml ${pepxml} \
			--protxml ${protxml}

		"""
    else if( params.skip_proteinprophet == false && params.combine_peptideprophet == false )
        """
        cd ${params.workspace}

		# filter matches and estimate FDR
		# skip the --sequential parameter
		philosopher filter \
			--psm 0.05 \
			--ion 0.05 \
			--pep 0.05 \
			--prot 1 \
			--razor \
			--picked \
			--tag rev_ \
			--pepxml interact-0*.pep.xml \
			--protxml interact.prot.xml

        """

}

process FREEQUANT {

    label "philosopher"

    publishDir "${params.philosopher_dir}/freequant", mode: 'copy'

    input:
    val filter_fdr

    output:
    val 'freequant_done_pseudo'

    script:
    """
    cd ${params.workspace}

	# Perform label-free quantification via 
	# precursor (MS1) abundances and spectral counting
    philosopher freequant \
		--dir . \
		2> freequant.out

    """

}

process REPORT {
    
    label "philosopher"

	echo true

    publishDir "${params.philosopher_dir}/report", mode: 'copy'

    input:
    val freequant

    output:
    path '*'
 
    script:
    """
	workd=\$(pwd)
    
	cd ${params.workspace}

	# reports about the findings for easy interpretation
	philosopher report \
		--msstats \
		2> report.out

	cp msstats.csv peptide.tsv psm.tsv ion.tsv  \$workd
	"""

}


process IONQUANT {

    label "ionquant"

    publishDir "${params.philosopher_dir}/ionquant", mode: 'copy'

    input:
    val filter_fdr
	path pepXML

    output:
    path "*_quant.csv"

    script:
    """
	workd=\$(pwd)
    cd ${params.workspace}

	# extract parent dir
	dir=${params.proteomics_reads}
    parentdir="\$(dirname "\$dir")"
	# Perform label-free quantification via 
	# precursor (MS1) abundances and spectral counting
	java -jar /IonQuant.jar \
		--specdir \${parentdir} \
		--multidir multidir \
		*.pepXML \
		&> ionquant.log

	cp *_quant.csv \$workd

    """

}

process CLEAN_UP_WORKSPACE {

    label "philosopher"

    input:
    path report

    script:
    """

    cd ${params.workspace}
    cd ..
    rm -rf ${params.workspace}

    """

}

workflow PHILOSOPHER_PIPE {

    take:
    ribotish_predict_ch
    proteomics_reads

    main:
    WORKSPACE(ribotish_predict_ch)

    DATABASE(WORKSPACE.out, ribotish_predict_ch)
    
	GENERATE_CHANGE_PARAMS(DATABASE.out, params.change_params_script)
    
	MSFRAGGER(
		proteomics_reads.collect(),
		GENERATE_CHANGE_PARAMS.out,
		DATABASE.out
	)
    
	PEPTIDEPROPHET(DATABASE.out, MSFRAGGER.out)
    
	PROTEINPROPHET(PEPTIDEPROPHET.out)

	FILTER_FDR(
		PEPTIDEPROPHET.out,
		PROTEINPROPHET.out
	)
	
	FREEQUANT(FILTER_FDR.out)
    
	REPORT(FREEQUANT.out)
    report = REPORT.out

	IONQUANT(
        FILTER_FDR.out,
        MSFRAGGER.out
    )

	//CLEAN_UP_WORKSPACE(REPORT.out)

    emit:
	report

}
