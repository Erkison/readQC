include { FASTP; CHECK_TRIMMED_FASTQ_SIZE; WRITE_OUT_EXCLUDED_READS; FASTQC; MULTIQC; 
        CONTAMINATION_CHECK; QUALIFYR; QUALIFYR_REPORT } from '../modules/processes.nf' 

workflow QC {
    take:
        reads_ch
        confindr_db_path
        qc_conditions
        posttrim_file_size_check
        version

    main: 
        FASTP(reads_ch)

        // exclude empty fastqs after trimming
        CHECK_TRIMMED_FASTQ_SIZE(FASTP.out.trimmed_fastqs)
        included_reads_based_on_size = CHECK_TRIMMED_FASTQ_SIZE.out.filter { it[1].toFloat() >= posttrim_file_size_check }
        excluded_reads_based_on_size = CHECK_TRIMMED_FASTQ_SIZE.out.filter { it[1].toFloat() < posttrim_file_size_check }
        included_reads = FASTP.out.trimmed_fastqs
            .join(included_reads_based_on_size)
            .map { items -> [items[0], items[1]] }
        excluded_reads = FASTP.out.trimmed_fastqs
            .join(excluded_reads_based_on_size)
            .map { items -> [items[0], items[1]] }

        WRITE_OUT_EXCLUDED_READS(excluded_reads)

        // Run rest of workflow for included reads
        FASTQC(included_reads)

        MULTIQC(FASTQC.out.fastqc_directories.collect())

        CONTAMINATION_CHECK(included_reads, confindr_db_path)

        // summarise quality
        quality_files = FASTQC.out.fastqc_report.join(CONTAMINATION_CHECK.out).join(included_reads)

        // quality_files.view()

        QUALIFYR(qc_conditions, quality_files)
        combined_qualifyr_json_files = QUALIFYR.out.json_files.collect()
        // combined_qualifyr_json_files.view()
        QUALIFYR_REPORT(combined_qualifyr_json_files, version)


    emit:
        included_reads
}