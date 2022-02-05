#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
https://github.com/Nesvilab/philosopher/wiki/Simple-Data-Analysis
*/


process WORKSPACE {

    label 'philosopher'

    input:
    path 'riboseq'

    output:
    val 'tubel'

    script:
    """
    rm -rf ${params.workspace}
    mkdir ${params.workspace}

    mzML_dir=${params.proteomics_reads}
    parentdir_mzML="\$(dirname "\$mzML_dir")"

    mv \${parentdir_mzML}/* ${params.workspace}/

    db_dir=${params.philosopher_db}
    parentdir_db="\$(dirname "\$db_dir")"

    mv ${params.philosopher_db} ${params.workspace}/

    # initialize philosopher workspace in workspace
    cd ${params.workspace}
    philosopher workspace --clean
    philosopher workspace --init

    cp *.mzML \$parentdir_mzML
    cp *.fasta \$parentdir_db
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

    cd ${params.workspace}
    philosopher database --custom ${db} #--contam
    
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
    path 'closed_fragger.params'

    script:
    """
    java -jar /MSFragger.jar --config

    python ${change_file_script} ${db} closed_fragger.params

    cp closed_fragger.params ${params.workspace}

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
    if (params.fragger_mode == "docker" )
	"""
	workd=\$(pwd)	

	cd ${params.workspace}
	
	for VAR in ${mzML_file}
	do

	java -Xmx${params.fragger_ram}g \
		-jar /MSFragger.jar \
		${closed_fragger} \${VAR} 

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
    path "*.pep.xml"
    
    script:
    """
    workd=\$(pwd)

    cd ${params.workspace}
    
    philosopher peptideprophet \
	--combine --database ${db} \
	--decoy rev_ --ppm --accmass \
	--expectscore --decoyprobs \
	--nonparam ${pepXML} \
	2> peptideprophet.out
    
    cp ${params.workspace}/interact.pep.xml \$workd

    """

}

process PROTEINPROPHET {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/proteinprophet", mode: 'copy'
       
    input:
    path pepxml
    
    output:
    path "*.prot.xml"
     
    script:
    if( params.skip_proteinprophet == true )
        """
        touch interact.prot.xml

        """

    else if( params.skip_proteinprophet == false )
        """
        workd=\$(pwd)

        cd ${params.workspace}
        philosopher proteinprophet ${pepxml} 2> proteinprophet.out

        cp ${params.workspace}/interact.prot.xml \$workd

        """

}

process FILTERANDFDR {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/filterandfdr", mode: 'copy'

    input:
    path pepxml
    path protXML

    output:
    val 'filterandfdr'
 
    script:
    if( params.skip_proteinprophet == true )
        """
        cd ${params.workspace}

        philosopher filter --picked --tag rev_ --pepxml ${pepxml}

        """

    else if( params.skip_proteinprophet == false )
        """
        cd ${params.workspace}

        philosopher filter --psm 0.05 --ion 0.05 --pep 0.05 --prot 1 --sequential --razor --picked --tag rev_ --pepxml interact.pep.xml --protxml interact.prot.xml

        """

}

process QUANTIFY {

    label "philosopher"

    publishDir "${params.philosopher_dir}/quantify", mode: 'copy'

    input:
    val filterandfdr

    output:
    val 'freequant'

    script:
    """
    cd ${params.workspace}
    philosopher freequant --dir . 2> freequant.out

    """

}

process REPORT {
    
    label "philosopher"

    publishDir "${params.philosopher_dir}/report", mode: 'copy'

    input:
    val freequant

    output:
    path '*'
 
    script:
    """
	workd=\$(pwd)
    cd ${params.workspace}
    philosopher report --msstats 2> report.out

	cp msstats.csv peptide.tsv psm.tsv ion.tsv  \$workd
    """

}

process CLEAN_UP_WORKSPACE {

    echo true

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
      db_ch

    main:
      workspace_obj = WORKSPACE(ribotish_predict_ch)
      db_obj = DATABASE(workspace_obj, db_ch)
      change_obj = GENERATE_CHANGE_PARAMS(db_obj, params.change_file_script)
      pepXML_obj = MSFRAGGER(proteomics_reads.collect(), change_obj, db_obj)
      interact_pep_xml = PEPTIDEPROPHET(db_obj, pepXML_obj)
      protXML_obj = PROTEINPROPHET(interact_pep_xml)
      filter_obj = FILTERANDFDR(interact_pep_xml, protXML_obj)
      quant_obj = QUANTIFY(filter_obj)
      report_obj = REPORT(quant_obj)
      CLEAN_UP_WORKSPACE(report_obj)

    emit:
      report_obj

}
