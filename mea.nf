nextflow.enable.dsl=2
params.geneColName = 'markname'
params.pvalColName = 'meta_p'
params.moduleFileDir = "/app/data/modules/cherryPickModules_noCoexpression/"
params.numRP = 3
// FIX BELOW PARAMS before running.
params.pipeline = "cma"
params.pvalFileName = "/app/data/pvals/cma/fhshdl.csv"
params.trait = "fhshdl"
params.numTests = 178106


process RandomPermutation {
    container 'edkang0925/mea-m1'
    publishDir ".", mode: 'copy'

    output:
    path("RPscores/${params.trait}/*.csv")

    script:
    """
    python3 /app/scripts/randomPermutation.py ${params.pvalFileName} "RPscores/${params.trait}/" ${params.geneColName} ${params.numRP}
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
    container 'edkang0925/mea-m1'

    input:
    path(pascalOutputFile)
    path(geneScoreFilePascalInput) // used to decide number of tests
    path(goFile)

    output:
    path("masterSummaryPiece/master_summary_slice_*")
    path("significantModules/")
    path(goFile)

    """
    python3 /app/scripts/processPascalOutput.py \
        ${pascalOutputFile} \
        0.05 \
        "masterSummaryPiece/" \
        ${geneScoreFilePascalInput} \
        "significantModules/" \
        ${params.numTests} 
    """
}

process GoAnalysis{
    container 'edkang0925/webgestalt-m1'
    publishDir ".", pattern:"GO_summaries/${params.trait}/*", mode: 'copy' // copy ORA results

    input:
    path(masterSummarySlice)
    path(sigModuleDir)
    path(goFile)


    output:
    path(masterSummarySlice)
    path("GO_summaries/${params.trait}/GO_summaries_${goFile.baseName.split('_')[2]}_${goFile.baseName.split('_')[3]}/")
    path(goFile)

    script:
    def oraSummaryDir = "GO_summaries/${params.trait}/GO_summaries_${goFile.baseName.split('_')[2]}_${goFile.baseName.split('_')[3]}/"
    """
    Rscript /app/scripts/ORA_cmd.R --sigModuleDir ${sigModuleDir} --backGroundGenesFile ${goFile} \
        --summaryRoot "${oraSummaryDir}" --reportRoot "GO_reports/"
    """

}

process MergeORAsummaryAndMasterSummary{
    container 'edkang0925/mea-m1'

    input:
    path(masterSummaryPiece)
    path(oraSummaryDir)
    path(goFile)

    output:
    path("mergedSummary/*")

    script:
    """
    python3 /app/scripts/mergeORAandSummary.py \
        ${masterSummaryPiece} \
        ${oraSummaryDir} \
        "mergedSummary/" \
        ${goFile}
    """

}
process VerticalMergeMasterSummaryPieces{
    container 'edkang0925/mea-m1'
    publishDir "./masterSummaries/", mode: 'copy'

    input:
    path(mergedSummaryFiles)

    output:
    path("master_summary_*.csv")

    script:
    """
    # Use the header from the first file
    head -n 1 \$(echo ${mergedSummaryFiles[0]}) > master_summary_${params.trait}.csv

    # Append the contents of each file excluding the header
    for file in ${mergedSummaryFiles.join(' ')}; do
        awk 'NR > 1' \$file >> master_summary_${params.trait}.csv
    done
    """
}

// process VerticalMergeMasterSummaryPieces{
//     container 'edkang0925/mea-m1'
//     publishDir "./masterSummaries/", mode: 'copy'

//     input:
//     path(mergedSummaryFiles)

//     output:
//     path("master_summary_*")

//     script:
//     def filesToMerge = task.workDir.resolve("merged_file_paths.txt")
//     file(filesToMerge).text = mergedSummaryFiles.collect { it.toString() }.join('\n')

//     """
//     python3 /app/scripts/verticalMerge.py $filesToMerge
//     """
// }



workflow {
    // For each module file in the module directory, preprocess the data for pascal.
    preProcessedFiles = PreProcessForPascal(RandomPermutation()|flatten)
    pascalOut = RunPascal(preProcessedFiles[0]|flatten, preProcessedFiles[1]|flatten, preProcessedFiles[2]|flatten)
    processedPascalOutput = ProcessPascalOutput(pascalOut[0]|flatten, pascalOut[1]|flatten, pascalOut[2]|flatten)
    goAnalysisOut = GoAnalysis(processedPascalOutput[0]|flatten, processedPascalOutput[1]|flatten, processedPascalOutput[2]|flatten)
    horizontallyMergedOut = MergeORAsummaryAndMasterSummary(goAnalysisOut[0]|flatten, goAnalysisOut[1]|flatten, goAnalysisOut[2]|flatten)
    VerticalMergeMasterSummaryPieces(horizontallyMergedOut.collect())
}


// workflow {

//     main_ch = RandomPermutation()
//              .flatten()
//              .PreProcessForPascal()
//              .RunPascal()
//              .ProcessPascalOutput()
//              .GoAnalysis()
//              .MergeORAsummaryAndMasterSummary()
    
//     VerticalMergeMasterSummaryPieces(main_ch.collect())
// }