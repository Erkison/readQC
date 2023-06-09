include { FASTP; FASTQC; MULTIQC; CONTAMINATION_CHECK; QUALIFYR; QUALIFYR_REPORT } from '../modules/processes.nf' 

workflow QC {
    take:
        reads_ch
    take: 
        confindr_db_path
    take: 
        qc_conditions
    take: 
        version

    main: 
        // reads_ch.view()
        FASTP(reads_ch)

        FASTQC(FASTP.out.trimmed_fastqs)

        MULTIQC(FASTQC.out.fastqc_directories.collect())

        CONTAMINATION_CHECK(FASTP.out.trimmed_fastqs, confindr_db_path)

        // summarise quality
        quality_files = FASTQC.out.fastqc_report.join(CONTAMINATION_CHECK.out).join(FASTP.out.trimmed_fastqs)

        // quality_files.view()

        QUALIFYR(qc_conditions, quality_files)
        combined_qualifyr_json_files = QUALIFYR.out.json_files.collect()
        combined_qualifyr_json_files.view()
        QUALIFYR_REPORT(combined_qualifyr_json_files, version)


    emit:
        FASTP.out.trimmed_fastqs
}