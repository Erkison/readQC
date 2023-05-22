nextflow.enable.dsl=2

// include non-process modules
include { help_message; version_message; complete_message; error_message; pipeline_start_message } from './modules/messages.nf'
include { default_params; check_params } from './modules/params_parser.nf'
include { help_or_version } from './modules/params_utils.nf'

version = '1.0dev'

// setup default params
default_params = default_params()
// merge defaults with user params
merged_params = default_params + params
// help and version messages
help_or_version(merged_params, version)
final_params = check_params(merged_params)
// starting pipeline
pipeline_start_message(version, final_params)

// include processes and workflows
include { QC } from './modules/workflows.nf'

workflow {
    // Setup input Channel from Read path
    reads_ch = Channel
        .fromFilePairs( final_params.reads, checkIfExists: true )
        .ifEmpty { error "Cannot find any reads matching: ${final_params.reads}" }
    confindr_db_path = Channel
        .fromPath( final_params.confindr_db_path, checkIfExists: true )
        .ifEmpty { error "Cannot find specified confindr db at: ${final_params.confindr_db_path}" }


    QC(reads_ch, confindr_db_path)
}

// Messages on completion 
workflow.onComplete {
    complete_message(final_params, workflow, version)
}

workflow.onError {
    error_message(workflow)
}