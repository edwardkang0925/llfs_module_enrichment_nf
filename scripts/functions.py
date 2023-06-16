import os
import re
import shutil
from typing import List
import functools as ft
import ast
from cmath import nan
import pandas as pd
import statsmodels.stats.multitest as smt
from statsmodels.sandbox.stats.multicomp import multipletests

def createOrCleanDir(DIRPATH:str, clean:bool=True) -> None:
    """
    Create or clean a directory

    Args:
        DIRPATH (str): path to the direcotry which needs to be created or cleaned
        clean (bool): flag indicating whether the pre-existing files under the directory should be deleted
    """
    if not os.path.exists(DIRPATH):
        print(f"{DIRPATH} does not exist. Creating ...")
        os.makedirs(DIRPATH)
    else:
        if clean:
            print(f"{DIRPATH} exists. Cleaning up previous files ...")
            for file in os.listdir(DIRPATH):
                os.remove(file)


def queryDirectories(DIRPATH:str) -> List[str]:
    """_summary_
        From a directory, return a list of paths to subdirectories
    
    Args:
        DIRPATH (str): path to the directory where the list of subdirectories are of interest

    Returns:
        List[str]: list of path of subdirectories under DIRPATH 
    """
    result = []
    if os.path.exists(DIRPATH):
        for item in os.listdir(DIRPATH):
            fullPath = os.path.join(DIRPATH, item)
            if os.path.isdir(fullPath):
                result.append(fullPath)
    return result


def querySpecificFiles(DIRPATH:str, endswith='.csv') -> List[str]:
    """
    Given a path to a directory, return list of full path to all csv files under the directory.

    Args:
        DIRPATH (str): Path to the directory containing csv files

    Returns:
        List[str]: List of full path to each of the csv files located under DIRPATH
    """
    result = []
    if os.path.exists(DIRPATH):
        for item in os.listdir(DIRPATH):
            if item.endswith(endswith):
                fullPath = os.path.join(DIRPATH, item)
                result.append(fullPath)
    else:
        print(f"{DIRPATH} does not exist but is queried for CSV files")
    return result
    
    

def extractCategoryFromFileName(DIRPATH:str, regex:str) -> str:
    """
    Given a path to a csv file, extract the category substring using the given regular expression

    Args:
        DIRPATH (str): path to a csv file containing gene score. 
        regex (str): regular expression to use to extract the category from the file name

    Returns:
        str: _description_
    """
    filename = DIRPATH.split("/")[-1]
    result = re.findall(regex, filename)
    return result[0]

def extractGeneSetFromModuleFile(DIRPATH:str):
    """
    Read a module file and extract a set of genes in the file. 
    It assumes the input file is tsv format where the gene name starts to appear from the thrid column

    Args:
        DIRPATH (str): path to the module file

    Returns:
        _type_: set of genes appear in the module file
    """
    ret = set()
    with open(DIRPATH, "r") as f:
        lines = f.readlines()
        for line in lines:
            columns = line.split()
            for column in columns[2:]:
                ret.add(column)
    return ret
            

def pairwiseProcessGeneScoreAndModule(GSPATH:str, MODULEPATH:str, OUTPUTPATH:str, GOPATH:str, pipeline:str, trait:str, geneNameCol:str, pvalCol:str, sep:str=',') -> None:
    """
    Given a pair of Gene score file and a module file, drop genes which does not exist in either of the files. 
    Write a pair of processed file with the same name. Pascal will take this pair as an input to proceed module enrichment.

    Args:
        GSPATH (str): path to the gene score file
        MODULEPATH (str): path to a pre-defined module file
        OUTPUTPATH (str): path to start building nested subdir for outputs. [pipeline > trait > output]
        pipeline (str): name of the pipeline. e.g. twas, gwas, staar, or cma
        trait (str): name of the trait
        geneNameCol (str): column name for gene name in the GS file
        pvalCol (str): column name for the pvalue in the GS file
        sep (str): if input gene score file is tab separated, pass in '\t' otherwise default should work.
    """
    
    df_gs = pd.read_csv(GSPATH, sep=sep) 
    genesWithScore = set(df_gs[geneNameCol])
    genesInModule = extractGeneSetFromModuleFile(MODULEPATH)
    intersectingGenes = genesWithScore.intersection(genesInModule)
    
    
    # create output directory
    outPvalDir = os.path.join(OUTPUTPATH, pipeline, trait, "pvals")
    createOrCleanDir(outPvalDir, clean=False) # create dir to output processed gs file
    GODir = os.path.join(GOPATH,pipeline, trait)
    createOrCleanDir(GODir, clean=False)
    outModuleDir = os.path.join(OUTPUTPATH, pipeline, trait, "modules")
    createOrCleanDir(outModuleDir, clean=False) # create dir to output processed gs file

    # output files
    gsFileName = GSPATH.split("/")[-1]
    moduleFileName = MODULEPATH.split("/")[-1]
    # output GS file to be used for PASCAL. 
    with open(os.path.join(outPvalDir, f"{pipeline}_{trait}_{moduleFileName[:-4]}.tsv"), "w") as f:
        for index, row in df_gs.iterrows():
            f.write(row[geneNameCol] + "\t" + str(row[pvalCol]) + "\n")
            
    # ouput GO background set 
    with open(os.path.join(GODir, f"{pipeline}_{trait}_{moduleFileName[:-4]}.txt"),"w") as f:
        for index, row in df_gs.iterrows():
            if row[geneNameCol] in intersectingGenes:
                f.write(f"{row[geneNameCol]}\n")
    
    # output module file after intersecting with Gene Score file
    with open(os.path.join(outModuleDir, f"{pipeline}_{trait}_{moduleFileName[:-4]}.tsv"), "w") as f:
        with open(MODULEPATH, "r") as g:
            droppedGeneCounter = 0
            lines = g.readlines()
            for line in lines:
                columns = line.split()
                f.write(columns[0]) # write the module index 
                for column in columns[2:]: # column[1] is always 1.0, so dropped
                    if column in intersectingGenes:
                        f.write("\t" + column)
                    else:
                        droppedGeneCounter += 1
                f.write("\n")
            print(f"Total {droppedGeneCounter} genes were dropped out of {len(genesInModule)} from {moduleFileName}")
            
    


def processOnePascalOutput(DIRPATH:str, alpha:float, outputPATH:str):
    """
    Given a path to a pascal output file, extract module index, module genes, and BH-corrected module pvalue

    Args:
        DIRPATH (str): path to a pascal output file
        alpha (float): significance threshold for modules pvalue after BH correction
        outputPATH (str): path to save the whole pascal result

    Returns:
        _type_: list of processed module info, total number of significant pathways
    """
    with open(DIRPATH, "r") as f:
        results = f.read()
        # flatten pval parts
        results = results.replace(",\n", ",")
        results = results.replace(" ", "")
        # 0: module index, 1: module genes, 2: gene uniform pval, 3: module uncorrected pval
        parsedResults = re.findall("\[(.+?),(.*?),array\((.*)\),(.*?)\]", results)
        
        pathwayIndexList = []
        pathwayGenesList = []
        pathwayPvalList = []
        
        for i in range(len(parsedResults)):
            # if a module lost all genes due to missing gene score, exclude it from FDR
            if parsedResults[i][3] != "nan": 
                # ex) "'5'" -> 5
                pathwayIndexList.append(int(parsedResults[i][0].replace("'", "")))
                # string list to list conversion via ast.literal_eval
                pathwayGenesList.append(ast.literal_eval(parsedResults[i][1]))
                pathwayPvalList.append(float(parsedResults[i][3]))
        
        # FDR correction BH or Bonferroni
        #correctedPathwayPvalList = smt.fdrcorrection(pathwayPvalList, alpha) # BH
        correctedPathwayPvalList = multipletests(pathwayPvalList, alpha, method='bonferroni') #Bonferroni
        
        # output csv file 
        df = pd.DataFrame(list(zip(pathwayIndexList, pathwayGenesList,
                           pathwayPvalList, correctedPathwayPvalList[1])),
                          columns=['moduleIndex', 'moduleGenes', 'modulePval', 'correctedModulePval'])
        df.to_csv(outputPATH)
        
        result = []
        
        numSigPathway = sum(correctedPathwayPvalList[0])
        for tup in zip(pathwayIndexList, pathwayGenesList,
                       correctedPathwayPvalList[0], correctedPathwayPvalList[1], pathwayPvalList):
            result.append(tup)
        # sort by corrected module pvalue
        result.sort(key=lambda x: x[-1])
        return result, numSigPathway
    
def extractGenesBasedOnPval(DIRPATH:str, pval:float):
    """
    Given a GeneScore file in tsv file format without header and a pvalue threshold, 
    output list of genes with pval less than the threshold

    Args:
        DIRPATH (str): Path to GS file in tsv format. ASSUMPTION: the first col is gene name and the second col is pval
        pval (float): pvalue threshold

    Returns:
        _type_: list of significant genes
    """
    df = pd.read_table(DIRPATH, header=None)
    return list(df[df[1] < pval][0])

            
def recordModulesFromPascalResult(result, OUTPUTPATH:str, sigGenesList, almostSigGenesList, sig4GenesList,
                                             sig3GenesList, sig2GenesList):
    moduleIndexToSize = {}
    moduleIndexToModulePval = {}
    moduleIndexToCorrectedModulePval = {}
    moduleIndexToSigFlag = {}
    moduleIndexSigGenes = {}
    moduleIndexAlmostSigGenes = {}
    moduleIndexToSig4Genes = {}
    moduleIndexToSig3Genes = {}
    moduleIndexToSig2Genes = {}
    # each item represents a module
    for item in result:
        sigGenes = []
        almostSigGenes = []
        sig4Genes = []
        sig3Genes = []
        sig2Genes = []
        moduleSizeCounter = 0
        for gene in item[1]:
            moduleSizeCounter += 1
            if gene in sigGenesList:
                sigGenes.append(gene)
            if gene in almostSigGenesList:
                almostSigGenes.append(gene)
            if gene in sig4GenesList:
                sig4Genes.append(gene)
            if gene in sig3GenesList:
                sig3Genes.append(gene)
            if gene in sig2GenesList:
                sig2Genes.append(gene)        
                
        # assumes index of 2 represents bool indicating significance of the module
        if item[2]:
            moduleIndexToSigFlag[item[0]] = True
            # assumes item[1] is list of genes
            saveSignificantModules(OUTPUTPATH.replace(".txt", f"_{item[0]}.txt"), item[1])
        else:
            moduleIndexToSigFlag[item[0]] = False
        moduleIndexToSize[item[0]] = moduleSizeCounter
        moduleIndexToModulePval[item[0]] = item[4]
        moduleIndexToCorrectedModulePval[item[0]] = item[3]
        moduleIndexSigGenes[item[0]] = sigGenes
        moduleIndexAlmostSigGenes[item[0]] = almostSigGenes
        moduleIndexToSig4Genes[item[0]] = sig4Genes
        moduleIndexToSig3Genes[item[0]] = sig3Genes
        moduleIndexToSig2Genes[item[0]] = sig2Genes

    return moduleIndexToSize, moduleIndexToModulePval, moduleIndexToCorrectedModulePval, moduleIndexToSigFlag, moduleIndexSigGenes, moduleIndexAlmostSigGenes, moduleIndexToSig4Genes, moduleIndexToSig3Genes, moduleIndexToSig2Genes
                    
                

def outputMergableORADF(ORAPATH:str, studies:List[str]):
    dict_ora = {'study':[], 'trait':[], 'network':[], 'moduleIndex':[],
                'geneontology_Biological_Process':[], 'BPminCorrectedPval':[],
                'BPminFDREnrichmentRatio':[], 'BPmaxEnrichmentRatio':[],
                'geneontology_Molecular_Function':[], "MFminCorrectedPval":[],
                'MFminFDREnrichmentRatio':[], 'MFmaxEnrichmentRatio':[]}
    GOtypes = ['geneontology_Biological_Process', 'geneontology_Molecular_Function']
    for study in studies:
        ora_trait_dirs = queryDirectories(os.path.join(ORAPATH, study))
        for ora_trait_dir in ora_trait_dirs:
            trait = ora_trait_dir.split("/")[-1]
            for ora_type in GOtypes:
                module_ora_files = querySpecificFiles(os.path.join(ora_trait_dir, ora_type), endswith=".csv")
                for module_ora_file in module_ora_files:
                    networkType = module_ora_file.split("_")[-2]
                    moduleIndex, GOcount, minPval, enrichmentRatio_minFDR, enrichmentRatio_max = countGOterms(module_ora_file)
                    # only append network and moduleIndex once for (study, trait) 
                    if ora_type == GOtypes[0]:
                        dict_ora['study'].append(study)
                        dict_ora["trait"].append(trait)
                        dict_ora["network"].append(networkType)
                        dict_ora["moduleIndex"].append(moduleIndex)
                        dict_ora['BPminCorrectedPval'].append(minPval)
                        dict_ora["BPminFDREnrichmentRatio"].append(enrichmentRatio_minFDR)
                        dict_ora["BPmaxEnrichmentRatio"].append(enrichmentRatio_max)
                    else:
                        dict_ora['MFminCorrectedPval'].append(minPval)
                        dict_ora['MFminFDREnrichmentRatio'].append(enrichmentRatio_minFDR)
                        dict_ora['MFmaxEnrichmentRatio'].append(enrichmentRatio_max)
                    dict_ora[ora_type].append(GOcount)

    return pd.DataFrame(dict_ora)
