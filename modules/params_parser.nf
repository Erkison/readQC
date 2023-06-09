include { check_mandatory_parameter; check_optional_parameters; check_parameter_value } from './params_utils.nf'

def default_params(){
    /***************** Setup inputs and channels ************************/
    def params = [:] as nextflow.script.ScriptBinding$ParamsMap
    // Defaults for configurable variables
    params.help = false
    params.version = false
    params.reads = false
    params.confindr_db_path = false
    params.output_dir = false
    params.qc_conditions = false

    return params
}

def check_params(Map params) { 
    final_params = params
    
    // set up reads files
    final_params.reads = check_mandatory_parameter(params, 'reads')
     
    // set up output directory
    final_params.output_dir = check_mandatory_parameter(params, 'output_dir') - ~/\/$/

    // check other mandatory params
    final_params.confindr_db_path = check_mandatory_parameter(params, 'confindr_db_path')

    final_params.qc_conditions = check_mandatory_parameter(params, 'qc_conditions')

      
    return final_params
}