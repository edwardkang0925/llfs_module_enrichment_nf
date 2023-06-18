#!/usr/bin/env Rscript

## TO USE:
## need one folder containing the gene lists of interest
## need another folder containing their reference gene lists
## all gene list files should be single-column .txt
## input-reference list pairs must share the first '.'-separated word of their names
## output files are named based on the first '_'-separated word of the input
## change FIELDS to desired name, method, database, gene id format, input folder, and gene universe folder
## data on the terms with significant overlap will be in SUMMARIES_PATH if such terms exist

library("WebGestaltR")
library(stringr)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("--sigModuleDir"), type="character",
                     help="directory having significant modules")
parser <- add_option(parser, c("--backGroundGenesFile"), type="character",
                     help="file having background genes")
parser <- add_option(parser, c("--summaryRoot"), type="character",
                     help="directory to save summary")
parser <- add_option(parser, c("--reportRoot"), type="character",
                     help="directory to save report")


opt = parse_args(parser)
print(opt$sigModuleDir)
METHOD = "ORA" # ORA | GSEA | NTA
DATABASES=c("geneontology_Biological_Process") 
GENE_ID="genesymbol" # see options with listIdType()

INPUT_PATH = file.path(opt$sigModuleDir)
# reports are more in-depth than summaries - advisable to keep reports FALSE if not needed
REPORTS_PATH= file.path(opt$reportRoot) # only used if GENERATE_REPORT=TRUE
SUMMARIES_PATH=file.path(opt$summaryRoot) # will be created if does not exist
GENERATE_REPORT=FALSE

# path must exist even if GENERATE_REPORT=FALSE
if (!dir.exists(REPORTS_PATH)) {
    dir.create(REPORTS_PATH, recursive=TRUE)
}

if (!dir.exists(SUMMARIES_PATH)) {
    dir.create(SUMMARIES_PATH, recursive=TRUE)
}

for(fileName in list.files(INPUT_PATH)){
    enrich_df <- NULL
    name <- ""
    if(grepl("sig_", fileName)){
        name <- sub("^(sig_.*)\\.txt$", "\\1", fileName) 
        tf_method = paste0(name, '_', METHOD)
        tryCatch(
            # perform enrichment analysis
            enrich_df <- WebGestaltR(
                enrichMethod = METHOD,
                organism = "hsapiens",
                enrichDatabase = DATABASE,
                interestGeneFile = significantModule,
                interestGeneType = GENE_ID,
                referenceGeneFile = opt$backGroundGenesFile,
                referenceGeneType = GENE_ID,
                minNum = 10, # default 10
                maxNum = 500, # default 500
                reportNum = 20, # default 20
                isOutput = GENERATE_REPORT,
                outputDirectory = REPORTS_PATH,
                projectName = tf_method
            ),
            error = function(e){
                print(paste0("ERROR while running WebGestalt for ",tf_method))
                enrich_df = NULL
            }
        )    
    }else{
        name <- sub("^(dummy_.*)\\.txt$", "\\1", fileName) 
    }
    # save summary as a .csv file
    if (!is.null(enrich_df)) {
    # remove link column
    sig_df <- subset(enrich_df, select = -c(link))
    # affinity propagation 
    idsInSet <- sapply(sig_df$overlapId, strsplit, split=";")
    names(idsInSet) <- sig_df$geneSet
    minusLogP <- -log(sig_df$pValue)
    minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
    apRes <- affinityPropagation(idsInSet, minusLogP)
    # subset GO terms for exemplar terms
    apGO_full <- sig_df[sig_df$geneSet %in% apRes$representatives,]
    if (nrow(apGO_full) > 0) {
        apGO_full['database'] <- rep(DATABASE, nrow(apGO_full))
        write.csv(apGO_full,file.path(SUMMARIES_PATH,paste0(name,".csv")),row.names = FALSE)
    } else {
        print("NO SIGNIFICANT OVERLAPS")
        write.csv(NULL,file.path(SUMMARIES_PATH,paste0(name,".csv")),row.names = FALSE)
    }
    } 
    else {
        print("NO SIGNIFICANT OVERLAPS")
        write.csv(NULL,file.path(SUMMARIES_PATH,paste0(name,".csv")),row.names = FALSE)
    }
        
    
}


