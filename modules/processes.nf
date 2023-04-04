process FASTQC {
    tag "FastQC on $sample_id"
    time { 1.hour * task.attempt }
    memory { 4.GB * task.attempt }

    publishDir "${params.output_dir}/fastqc/pre_trimming",
    mode: 'copy',
    pattern: "*.html"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.html")
    tuple(val(sample_id), path("${sample_id}_R*_fastqc.txt"),   emit: qc_post_trimming_files)
    path("*_fastqc_data"),                                      emit: fastqc_directories


    script:
    r1_prefix = reads[0].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
    r2_prefix = reads[1].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
    """
    fastqc ${reads[0]} ${reads[1]} --extract
    mv ${r1_prefix}_fastqc/summary.txt ${sample_id}_R1_fastqc.txt
    mv ${r2_prefix}_fastqc/summary.txt ${sample_id}_R2_fastqc.txt

    # move files for fastqc
    mkdir ${r1_prefix}_fastqc_data
    mkdir ${r2_prefix}_fastqc_data
    mv ${r1_prefix}_fastqc/fastqc_data.txt ${r1_prefix}_fastqc_data
    mv ${r2_prefix}_fastqc/fastqc_data.txt ${r2_prefix}_fastqc_data
    """
}

process MULTIQC {
    tag { 'multiqc for fastqc' }
    time { 1.hour * task.attempt }
    memory { 4.GB * task.attempt }

    publishDir "${params.output_dir}/",
    mode: 'copy',
    pattern: "multiqc_report.html",
    saveAs: { "fastqc_multiqc_report.html" }

    input:
    path(fastqc_directories) 

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc --interactive .
    """
}
