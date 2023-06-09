process FASTP {
    tag "Fastp on $sample_id"

    label 'medium_resource_req'

    publishDir "${params.output_dir}/fastp", 
        mode: 'copy',
        pattern: "*.fastp.*"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple(val(sample_id),  path("*.trim.fastq.gz"), emit: trimmed_fastqs)
        path("*.fastp.*")
        path("*.fail.fastq.gz")

    script:
    """
    [ ! -f  ${sample_id}_1.fastq.gz ] && ln -s ${reads[0]} ${sample_id}_1.fastq.gz
    [ ! -f  ${sample_id}_2.fastq.gz ] && ln -s ${reads[1]} ${sample_id}_2.fastq.gz
    INPUT_READS='--in1 ${sample_id}_1.fastq.gz --in2 ${sample_id}_2.fastq.gz'
    OUTPUT_READS='--out1 ${sample_id}_1.trim.fastq.gz --out2 ${sample_id}_2.trim.fastq.gz --unpaired1 ${sample_id}_1.fail.fastq.gz --unpaired2 ${sample_id}_2.fail.fastq.gz'

    # run fastp
    fastp \$INPUT_READS \$OUTPUT_READS \\
        --cut_front --cut_tail \\
        --cut_mean_quality 30 \\
        --length_required  50 \\
        --thread $task.cpus \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html
    """
}

process FASTQC {
    tag "FastQC on $sample_id"

    label 'medium_resource_req'

    publishDir "${params.output_dir}/fastqc/",
        pattern: "${sample_id}_fastqc_outputs",
        mode: 'copy'

    input:
    tuple(val(sample_id), path(reads))

    output:
    path("*_fastqc_outputs"), emit: fastqc_directories
    tuple(val(sample_id), path("${sample_id}_*_fastqc.txt"), emit: fastqc_report)


    script:
    """
    mkdir ${sample_id}_fastqc_outputs
    fastqc --quiet --extract --threads $task.cpus *.trim.fastq.gz --outdir ${sample_id}_fastqc_outputs
    cp ${sample_id}_fastqc_outputs/${sample_id}_1.trim_fastqc/summary.txt ${sample_id}_1_fastqc.txt
    cp ${sample_id}_fastqc_outputs/${sample_id}_2.trim_fastqc/summary.txt ${sample_id}_2_fastqc.txt
    """
}

process MULTIQC {
    tag { 'MultiQC for fastqc' }

    label 'medium_resource_req'

    publishDir "${params.output_dir}/multiqc",
        mode: 'copy',
        pattern: "multiqc_report.html",
        saveAs: { "fastqc_multiqc_report.html" }

    publishDir "${params.output_dir}/multiqc",
        mode: 'copy',
        pattern: "multiqc_data/multiqc_fastqc.txt",
        saveAs: { "fastqc_multiqc_report.tsv" }

    input:
    path(fastqc_directories) 

    output:
    path("multiqc_report.html"), emit: multiqc_html_report
    path("multiqc_data/multiqc_fastqc.txt"), emit: multiqc_tsv_report

    script:
    """
    multiqc --interactive .
    """
}

process CONTAMINATION_CHECK {
    tag { "ConFindr on ${sample_id}" }

    label 'big_resource_req'

    publishDir "${params.output_dir}/confindr",
        mode: 'copy',
        saveAs: { file -> "${sample_id}_${file}" }

    input:
    tuple( val(sample_id), path(trimmed_reads) )
    each(confindr_db_path)

    output:
    tuple( val(sample_id), path("confindr_report.csv") )

    script:
    """
    confindr.py -i . -o . -d ${confindr_db_path} -t $task.cpus \\
        -bf 0.025 -b 2 --cross_detail \\
        -fid _1 -rid _2
    """
}

process QUALIFYR {
    tag "Qualifyr on $sample_id"

    label 'small_resource_req'

    publishDir "${params.output_dir}/trimmed_reads/pass",
        mode: 'copy',
        pattern: 'trimmed_reads/pass/*',
        saveAs: { file -> file.split('\\/')[-1] }

    publishDir "${params.output_dir}/trimmed_reads/warning",
        mode: 'copy',
        pattern: 'trimmed_reads/warning/*',
        saveAs: { file -> file.split('\\/')[-1] }
    
    publishDir "${params.output_dir}/trimmed_reads/failure",
        mode: 'copy',
        pattern: 'trimmed_reads/failure/*',
        saveAs: { file -> file.split('\\/')[-1] }

    publishDir "${params.output_dir}/qualifyr/warn_and_fail_reports/",
        mode: 'copy',
        pattern: 'qualifyr/warn_and_fail_reports/*',
        saveAs: { file -> file.split('\\/')[-1] }

    input:
    each(qc_conditions_yml)
    tuple(val(sample_id), path(fastqc_reports), path(confindr_report), path(reads))

    output:
    path('trimmed_reads/**/*')
    path('qualifyr/warn_and_fail_reports/*'), optional: true
    path("${sample_id}.qualifyr.json"), emit: json_files

    """
    mkdir -p qualifyr/warn_and_fail_reports
    result=`qualifyr check -y ${qc_conditions_yml} -f ${fastqc_reports} -c ${confindr_report} -s ${sample_id} 2> ERR`
    return_code=\$?
    if [[ \$return_code -ne 0 ]]; then
        exit 1;
    else
        if [[ \$result == "PASS" ]]; then
            qc_level="pass"
        elif [[ \$result == "WARNING" ]]; then
            qc_level="warning"
        elif [[ \$result == "FAILURE" ]]; then
            qc_level="failure"
        fi

        mkdir -p trimmed_reads/\${qc_level}
        mv ${reads} trimmed_reads/\${qc_level}/

        if [[ \$result != "PASS" ]]; then
            mv ERR qualifyr/warn_and_fail_reports/${sample_id}_\${qc_level}_qc_result.tsv
        fi
    fi

    # make json file
    qualifyr check -y ${qc_conditions_yml} -f ${fastqc_reports} -c ${confindr_report} -s ${sample_id} -j -o .
    """

}

process QUALIFYR_REPORT {
    tag { 'qualifyr report' }

    label 'small_resource_req'

    publishDir "${params.output_dir}/qualifyr",
        mode: 'copy',
        pattern: "qualifyr_report.*"

    input:
    path(json_files)
    val(version)

    output:
    path("qualifyr_report.*")

    script:
    """
    qualifyr report -i . -c 'confindr.contam_status' -t "Analysis with readQC Pipeline version ${version}"
    """

}

