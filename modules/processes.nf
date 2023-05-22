process FASTP {
    tag "Fastp on $sample_id"

    label 'small_resource_req'

    publishDir "${params.output_dir}/fastp", 
        mode: 'copy',
        pattern: "*.fastp.*"
    
    publishDir "${params.output_dir}/trimmed_fastqs", 
        mode: 'copy',
        pattern: "*.trim.fastq.gz"

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

    label 'small_resource_req'

    publishDir "${params.output_dir}/fastqc/",
        mode: 'copy'

    input:
    tuple(val(sample_id), path(reads))

    output:
    path("*_fastqc_outputs"), emit: fastqc_directories


    script:
    """
    mkdir ${sample_id}_fastqc_outputs
    fastqc --quiet --extract --threads $task.cpus *.trim.fastq.gz --outdir ${sample_id}_fastqc_outputs
    """
}

process MULTIQC {
    tag { 'MultiQC for fastqc' }

    label 'medium_resource_req'

    publishDir "${params.output_dir}/multiqc",
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
