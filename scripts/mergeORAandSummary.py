import argparse
import os
import re
from typing import List
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


def outputMergableORA_df(module_ora_file:str, study:str, trait:str, network:str, moduleIndex:int):
    # if "dummy" in module_ora_file:
    #     return pd.DataFrame({'study':[study], 'trait':[trait], 'network':[network], 'moduleIndex':[moduleIndex],
    #                          'geneontology_Biological_Process':["NA"], 'BPminCorrectedPval':["NA"],
    #                           "BPminFDREnrichmentRatio":["NA"], 'BPmaxEnrichmentRatio':["NA"]})
                             
    dict_ora = {'study':[], 'trait':[], 'network':[], 'moduleIndex':[],
                'geneontology_Biological_Process':[], 'BPminCorrectedPval':[],
                'BPminFDREnrichmentRatio':[], 'BPmaxEnrichmentRatio':[]}
    GOtype = 'geneontology_Biological_Process'
    moduleIndex, GOcount, minPval, enrichmentRatio_minFDR, enrichmentRatio_max = countGOterms(module_ora_file)
    dict_ora['study'].append(study)
    dict_ora["trait"].append(trait)
    dict_ora["network"].append(network)
    dict_ora["moduleIndex"].append(moduleIndex)
    dict_ora['BPminCorrectedPval'].append(minPval)
    dict_ora["BPminFDREnrichmentRatio"].append(enrichmentRatio_minFDR)
    dict_ora["BPmaxEnrichmentRatio"].append(enrichmentRatio_max)
    dict_ora[GOtype].append(GOcount)

    return pd.DataFrame(dict_ora)
    

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="merge horizontally summary + ora summary")
    # Add arguments to parser
    parser.add_argument("masterSummaryPiece", help="pathToSummaryResult")
    parser.add_argument("oraResultsDir", help="path To Dir containing ORAResults")
    parser.add_argument("output_directory", help="path To save OutputMergedFile")
    parser.add_argument("goFile", help="path to GO background file")

    # Parse the arguments
    args = parser.parse_args()
        
    print(args.oraResultsDir)
    print(args.masterSummaryPiece)
    
    # Check if the output directory exists, if not create it
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
        
    df_summary_piece = pd.read_csv(args.masterSummaryPiece)
    df_summary_piece[['study', 'trait', 'network', 'moduleIndex']] = df_summary_piece[['study', 'trait', 'network', 'moduleIndex']].astype(str)
    study, trait, network = os.path.basename(args.goFile).split(".")[0].split("_")[1:4]

    
    ora_dfs = []
    if len(os.listdir(args.oraResultsDir)) == 0: # if there are no significant module from ORA result
        ora_dfs.append(pd.DataFrame({'study':[study], 'trait':[trait], 'network':[network], 'moduleIndex':[0],
                                    'geneontology_Biological_Process':["NA"], 'BPminCorrectedPval':["NA"],
                                    "BPminFDREnrichmentRatio":["NA"], 'BPmaxEnrichmentRatio':["NA"]}))
    else:
        for file in os.listdir(args.oraResultsDir):
            file = os.path.join(args.oraResultsDir, file)
            moduleIndex = os.path.basename(file).split(".")[0].split("_")[4]
            df_ora = outputMergableORA_df(file, study, trait, network, moduleIndex)
            ora_dfs.append(df_ora)
            
    df_ora_merged = pd.concat(ora_dfs, ignore_index=True)
    df_ora_merged[['study', 'trait', 'network', 'moduleIndex']] = df_ora_merged[['study', 'trait', 'network', 'moduleIndex']].astype(str)
    df_merge = pd.merge(df_summary_piece, df_ora_merged, how='left', on=['study','trait','network', 'moduleIndex'])
    df_merge.fillna("NA", inplace=True)
    mergedFileName = f"{study}_{trait}_{network}.csv"
    df_merge.to_csv(os.path.join(args.output_directory, mergedFileName), index=False)
    

if __name__ == "__main__":
    main()
