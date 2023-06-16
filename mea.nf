params.geneColName = 'markname'
params.pvalColName = 'meta_p'
params.pvalFileName = "/app/data/pvals/cma/ABI.csv"
params.moduleFileDir = "/app/data/modules/cherryPickModules_noCoexpression/"
params.pipeline = "cma"
params.trait = "fhshdl"
params.numRP = 2

nextflow.enable.dsl=2

process RandomPermutation {
    container 'edkang0925/mea-m1'

    output:
    path("outputs/RP/*.csv")

    script:
    """
    python3 /app/scripts/randomPermutation.py ${params.pvalFileName} "outputs/RP/" ${params.geneColName} ${params.numRP}
    """
}

process PreProcessForPascal{
    container 'edkang0925/mea-m1'

    input:
    path(geneScoreFile)

    output:
    path("pascalInput/GS_*")
    path("pascalInput/Module_*")
    path("pascalInput/GO_*")

    script:
    """
    python3 /app/scripts/preProcessForPascal.py \
        ${geneScoreFile} \
        ${params.moduleFileDir} \
        "pascalInput/" \
        ${params.pipeline} \
        ${params.trait} \
        ${params.geneColName} \
        ${params.pvalColName}
    """

}

process RunPascal{
    container 'acharyas/pascalx'

    input:
    path(geneScoreFile)
    path(moduleFile)

    output:
    path("pascalOutput/*")

    script:
    """
    python3 /app/scripts/runPascal.py \
        ${geneScoreFile} \
        ${moduleFile} \
        "pascalOutput/" \
        ${params.pipeline} \
        ${params.trait}
        
    """
}

process ProcessPascalOutput{
    container 'edkang0925/mea-m1'

    input:
    path(pascalOutputFile)
    path(geneScoreFilePascalInput) // used to decide number of tests

    output:
    path("masterSummaryPiece/*")
    path("significantModules/")

    """
    python3 /app/scripts/processPascalOutput.py \
        ${pascalOutputFile} \
        0.05 \
        "masterSummaryPiece/" \
        ${geneScoreFilePascalInput} \
        "significantModules/"
    """
}

process GoAnalysis{
    container 'edkang0925/webgestalt-m1'

    input:
    path(sigModuleDir)
    path(backGroundGenesDir)

    output:
    path("GO_summaries/")

    """
    Rscript /app/scripts/ORA_cmd.R --sigModuleDir ${sigModuleDir} --backGroundGenesFile ${backGroundGenesDir} \
        --summaryRoot "GO_summaries/" --reportRoot "GO_reports/"

    """
}

process MergeORAsummaryAndMasterSummary{
    container 'edkang0925/mea-m1'

    input:
    path(oraSummaryDir)
    path(masterSummaryPiece)

    """
    python3 /app/scripts/mergeORAandSummary.py \
        ${oraSummaryDir} \
        ${masterSummaryPiece} 
    """

}



workflow {
    // For each module file in the module directory, preprocess the data for pascal.
    preProcessedFiles = PreProcessForPascal(RandomPermutation()|flatten)
    pascalOut = RunPascal(preProcessedFiles[0]|flatten, preProcessedFiles[1]|flatten)
    processedPascalOutput = ProcessPascalOutput(pascalOut[0]|flatten, preProcessedFiles[0]|flatten)
    goAnalysisOut = GoAnalysis(processedPascalOutput[1]|flatten, preProcessedFiles[2]|flatten)
    MergeORAsummaryAndMasterSummary(goAnalysisOut[0]|flatten, processedPascalOutput[0]|flatten)

}
