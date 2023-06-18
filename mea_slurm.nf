params.geneColName = 'markname'
params.pvalColName = 'meta_p'
params.pvalFileName = "/app/data/pvals/cma/fhshdl.csv"
params.moduleFileDir = "/app/data/modules/cherryPickModules_noCoexpression/"
params.pipeline = "cma"
params.trait = "fhshdl"
params.numRP = 50

nextflow.enable.dsl=2

process RandomPermutation {
    container 'mea_latest.sif'
    label "process_low"

    output:
    path("outputs/RP/*.csv")

    script:
    """
    python3 /app/scripts/randomPermutation.py ${params.pvalFileName} "outputs/RP/" ${params.geneColName} ${params.numRP}
    """
}

process PreProcessForPascal{
    container 'mea_latest.sif'
    label "process_low"

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
    container 'pascalx_latest.sif'
    label "process_low"

    input:
    path(geneScoreFile)
    path(moduleFile)
    path(goFile)

    output:
    path("pascalOutput/*")
    path(geneScoreFile)
    path(goFile)

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
    container 'mea_latest.sif'
    label "process_low"

    input:
    path(pascalOutputFile)
    path(geneScoreFilePascalInput) // used to decide number of tests
    path(goFile)

    output:
    path("masterSummaryPiece/master_summary_slice_*")
    path("significantModules/")
    path(geneScoreFilePascalInput)
    path(goFile)

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
    container 'webgestalt_latest.sif'
    label "process_low"

    input:
    path(masterSummarySlice)
    path(sigModuleDir)
    path(geneScoreFilePascalInput) // used to decide number of tests
    path(goFile)

    output:
    path(masterSummarySlice)
    path("GO_summaries/")
    path(geneScoreFilePascalInput) // used to decide number of tests
    path(goFile)

    """
    Rscript /app/scripts/ORA_cmd.R --sigModuleDir ${sigModuleDir} --backGroundGenesFile ${goFile} \
        --summaryRoot "GO_summaries/" --reportRoot "GO_reports/"

    """
}

process MergeORAsummaryAndMasterSummary{
    container 'mea_latest.sif'
    label "process_low"
    input:
    path(oraSummaryDir)
    path(masterSummaryPiece)

    output:
    path("mergedSummary/*")

    """
    python3 /app/scripts/mergeORAandSummary.py \
        ${oraSummaryDir} \
        ${masterSummaryPiece} \
        "mergedSummary/"
    """

}

process VerticalMergeMasterSummaryPieces{
    container 'mea_latest.sif'
    publishDir "./masterSummaries/", mode: 'copy'
    label "process_medium"
    input:
    path(mergedSummaryFiles)

    output:
    path("master_summary_*")

    """
    python3 /app/scripts/verticalMerge.py \
        ${mergedSummaryFiles} 
    """
    
}



workflow {
    // For each module file in the module directory, preprocess the data for pascal.
    preProcessedFiles = PreProcessForPascal(RandomPermutation()|flatten)
    pascalOut = RunPascal(preProcessedFiles[0]|flatten, preProcessedFiles[1]|flatten, preProcessedFiles[2]|flatten)
    processedPascalOutput = ProcessPascalOutput(pascalOut[0]|flatten, pascalOut[1]|flatten, pascalOut[2]|flatten)
    goAnalysisOut = GoAnalysis(processedPascalOutput[0]|flatten, processedPascalOutput[1]|flatten, processedPascalOutput[2]|flatten,processedPascalOutput[3]|flatten)
    horizontallyMergedOut = MergeORAsummaryAndMasterSummary(goAnalysisOut[0]|flatten, goAnalysisOut[1]|flatten)
    VerticalMergeMasterSummaryPieces(horizontallyMergedOut.collect())
}
