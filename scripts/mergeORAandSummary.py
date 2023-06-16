import argparse
import os
import re
from typing import List
from .preprocess import *
import pandas as pd

def countGOterms(DIRPATH:str)-> int:
    df = pd.read_csv(DIRPATH)
    moduleIndex = DIRPATH.split("/")[-1].split("_")[-1].replace(".csv","")
    maxEnrichmentRatio = -1
    enrichMentRatio = -1
    if len(df) > 0:
        minTermPval = df["FDR"].min()
        maxEnrichmentRatio = df["enrichmentRatio"].max()
        enrichMentRatio = df[df["FDR"] == df["FDR"].min()]["enrichmentRatio"].max()
    else:
        minTermPval = -1
        maxEnrichmentRatio = -1
        enrichMentRatio = -1
    return int(moduleIndex), len(df), minTermPval, enrichMentRatio, maxEnrichmentRatio


def outputMergableORA_df(module_ora_file:str, study:str, trait:str, network:str, moduleIndex:int,module_ora_file:str):
    if "dummy" in module_ora_file:
        return pd.DataFrame({'study':[study], 'trait':[trait], 'network':[network], 'moduleIndex':[moduleIndex],
                             'geneontology_Biological_Process':["NA"], 'BPminCorrectedPval':["NA"],
                             'geneontology_Molecular_Function':["NA"], "MFminCorrectedPval":["NA"],
                                'MFminFDREnrichmentRatio':["NA"], 'MFmaxEnrichmentRatio':["NA"]})
                             
    dict_ora = {'study':[], 'trait':[], 'network':[], 'moduleIndex':[],
                'geneontology_Biological_Process':[], 'BPminCorrectedPval':[],
                'BPminFDREnrichmentRatio':[], 'BPmaxEnrichmentRatio':[],
                'geneontology_Molecular_Function':[], "MFminCorrectedPval":[],
                'MFminFDREnrichmentRatio':[], 'MFmaxEnrichmentRatio':[]}
    GOtypes = 'geneontology_Biological_Process'
    moduleIndex, GOcount, minPval, enrichmentRatio_minFDR, enrichmentRatio_max = countGOterms(module_ora_file)
    dict_ora['study'].append(study)
    dict_ora["trait"].append(trait)
    dict_ora["network"].append(network)
    dict_ora["moduleIndex"].append(moduleIndex)
    dict_ora['BPminCorrectedPval'].append(minPval)
    dict_ora["BPminFDREnrichmentRatio"].append(enrichmentRatio_minFDR)
    dict_ora["BPmaxEnrichmentRatio"].append(enrichmentRatio_max)
    dict_ora[ora_type].append(GOcount)

    return pd.DataFrame(dict_ora)
    

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Randomly permute the first column of a CSV file.")
    # Add arguments to parser
    parser.add_argument("oraResult", help="pathToORAResult")
    parser.add_argument("masterSummaryPiece", help="pathToSummaryResult")
    parser.add_argument("outputMergedFilePath", help="pathToOutputMergedFile")

    # Parse the arguments
    args = parser.parse_args()
        
    print(args.oraResult)
    print(args.masterSummaryPiece)
    
    study, trait, network, moduleIndex = os.path.basename(args.oraResult).split("_")[0:4]
    df_ora = outputMergableORADF(args.oraResult, study, trait, network, moduleIndex)
    df_merge = pd.merge(args.masterSummaryPiece, df_ora, how='left', on=['study','trait','network', 'moduleIndex'])
    df_merge.fillna("NA", inplace=True)
    df_merge.to_csv(args.outputMergedFilePath, index=False)