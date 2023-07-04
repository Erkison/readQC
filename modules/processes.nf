process FASTP {
    tag "Fastp on $sample_id"

    label 'medium_resource_req'

    publishDir "${params.output_dir}/fastp", 
        mode: 'copy',
        pattern: "*.fastp.html"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple(val(sample_id),  path("*_trimmed_*.fastq.gz"), emit: trimmed_fastqs)
        path("*.fastp.*")
        path("*failed_*.fastq.gz")

    script:
    """
    [ ! -f  ${sample_id}_1.fastq.gz ] && ln -s ${reads[0]} ${sample_id}_1.fastq.gz
    [ ! -f  ${sample_id}_2.fastq.gz ] && ln -s ${reads[1]} ${sample_id}_2.fastq.gz
    INPUT_READS='--in1 ${sample_id}_1.fastq.gz --in2 ${sample_id}_2.fastq.gz'
    OUTPUT_READS='--out1 ${sample_id}_trimmed_1.fastq.gz --out2 ${sample_id}_trimmed_2.fastq.gz --unpaired1 ${sample_id}_failed_1.fastq.gz --unpaired2 ${sample_id}_failed_2.fastq.gz'

    # run fastp
    fastp \$INPUT_READS \$OUTPUT_READS \\
        --cut_front --cut_tail \\
        --cut_mean_quality 20 \\
        --length_required  50 \\
        --thread $task.cpus \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html
    """
}

process CHECK_TRIMMED_FASTQ_SIZE {
    tag { sample_id }
    
    input:
    tuple(val(sample_id), path(reads))

    output:
    tuple(val(sample_id), stdout)

    script:
    """
    size=\$(stat -Lc %s ${reads[0]} 2>/dev/null |  awk '{printf "%f", \$1/1000000}')
    # if stat failed (on MacOS), use MacOS stat
    if [ "\$size" == "" ]; then
      size=\$(stat -Lf%z ${reads[0]} 2>/dev/null |  awk '{printf "%f", \$1/1000000}')
    fi
    echo \$size
    """
}

process WRITE_OUT_EXCLUDED_READS {
    tag "$sample_id"

    publishDir "${params.output_dir}/trimmed_reads/failure/",
        mode: 'move'

    input:
    tuple(val(sample_id), path(reads))

    output:
    path("too_small_size_post_trimming/*_trimmed_*.fastq.gz")


    script:
    """
    mkdir too_small_size_post_trimming
    cp ${reads} too_small_size_post_trimming
    """
}

process FASTQC {
    tag "FastQC on $sample_id"

    label 'medium_resource_req'

    // publishDir "${params.output_dir}/fastqc/",
    //     pattern: "${sample_id}_fastqc_outputs",
    //     mode: 'copy'

    input:
    tuple(val(sample_id), path(reads))

    output:
    path("*_fastqc_outputs"), emit: fastqc_directories
    tuple(val(sample_id), path("${sample_id}_*_fastqc.txt"), emit: fastqc_report)


    script:
    """
    mkdir ${sample_id}_fastqc_outputs
    fastqc --quiet --extract --threads $task.cpus *_trimmed_*.fastq.gz --outdir ${sample_id}_fastqc_outputs
    cp ${sample_id}_fastqc_outputs/${sample_id}_trimmed_1_fastqc/summary.txt ${sample_id}_1_fastqc.txt
    cp ${sample_id}_fastqc_outputs/${sample_id}_trimmed_2_fastqc/summary.txt ${sample_id}_2_fastqc.txt
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

    publishDir "${params.output_dir}/qualifyr/",
        mode: 'copy',
        pattern: '*_sample_ids.txt'

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

    input:
    each(qc_conditions_yml)
    tuple(val(sample_id), path(fastqc_reports), path(confindr_report), path(reads))

    output:
    path('trimmed_reads/**/*')
    path('*_sample_ids.txt')
    path("${sample_id}.qualifyr.json"), emit: json_files

    """
    result=`qualifyr check -y ${qc_conditions_yml} -f ${fastqc_reports} -c ${confindr_report} -s ${sample_id} 2> ERR`
    return_code=\$?
    if [[ \$return_code -ne 0 ]]; then
        exit 1;
    else
        if [[ \$result == "PASS" ]]; then
            qc_level="pass"
            echo ${sample_id} >> pass_sample_ids.txt
        elif [[ \$result == "WARNING" ]]; then
            qc_level="warning"
            echo ${sample_id} >> pass_sample_ids.txt
        elif [[ \$result == "FAILURE" ]]; then
            qc_level="failure"
            echo ${sample_id} >> fail_sample_ids.txt
        fi

        mkdir -p trimmed_reads/\${qc_level}
        mv ${reads} trimmed_reads/\${qc_level}/

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

