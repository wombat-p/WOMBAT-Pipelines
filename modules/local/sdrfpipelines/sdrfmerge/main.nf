//import org.yaml.snakeyaml.Yaml
process SDRFMERGE { 
      
    label 'process_medium'
    conda (params.enable_conda ? "bioconda::sdrf-pipelines=0.0.21" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://wombatp/maxquant-pipeline:dev"
    } else {
        container "wombatp/maxquant-pipeline:dev"
    }

    publishDir "${params.outdir}/sdrf_merge", mode:'copy'


    input:
      path sdrf
      file parameters
      path map
     

    output:
      path "sdrf_local.tsv"         , emit: sdrf_local
      path "params.yml"          , emit: parameters_out
      //value parameters_out, emit: parameters


    script:
    // load parameters
    //def par_filename = file( ["${task.workDir}", "${parameters}"].join(File.separator) )
    //def par_filename = parameters
    //println(par_filename)
    //def yaml = (Map)new Yaml().load((par_filename).text)    

    """
    if [[ "$sdrf" != "sdrf.tsv" ]]
    then
	cp "${sdrf}" sdrf.tsv
    fi
    if [[ "$parameters" != "params.yml" ]] 
    then
        cp "${parameters}" params.yml
    fi
    if [[ "$map" != "params2sdrf.yml" ]]
    then
        cp "${map}" params2sdrf.yml
    fi
    # TODO change to package when available
    python $projectDir/scripts/add_data_analysis_param.py
    python $projectDir/scripts/sdrf2params.py
    echo "preliminary version" > sdrf-merge.version.txt
    """
}
