include { FASTP; FASTQC; MULTIQC; CONTAMINATION_CHECK } from '../modules/processes.nf' 

workflow QC {
    take:
        reads_ch
    take: 
        confindr_db_path

    main: 
        // reads_ch.view()
        FASTP(reads_ch)

        FASTQC(FASTP.out.trimmed_fastqs)

        MULTIQC(FASTQC.out.fastqc_directories.collect())

        CONTAMINATION_CHECK(FASTP.out.trimmed_fastqs, confindr_db_path)

    emit:
        MULTIQC.out
}