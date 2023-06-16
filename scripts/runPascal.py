import argparse
import os
import glob
from typing import List

from PascalX import pathway
from PascalX import genescorer

        
def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Randomly permute the first column of a CSV file.")
    
    # Add arguments to parser
    parser.add_argument("scoreFile", help="Path to the scoreFile.")
    parser.add_argument("moduleFile", help="Path to the moduleFile.")
    parser.add_argument("outputPath", help="Path to the output directory.")
    parser.add_argument("pipelineName", help="Name of the pipeline.")
    parser.add_argument("traitName", help="Name of the trait.")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Check if the output directory exists, if not create it
    if not os.path.exists(args.outputPath):
        os.makedirs(args.outputPath)
    
    #for moduleFile, scoreFile in zip(moduleFiles, scoreFiles):
    Scorer = genescorer.chi2sum()
    Scorer.load_scores(args.scoreFile)
    Pscorer = pathway.chi2rank(Scorer, fuse=False)
    M = Pscorer.load_modules(args.moduleFile, ncol=0, fcol=1)
    RESULT = Pscorer.score(M)
    fileName = os.path.basename(args.scoreFile).replace("tsv", "txt").replace("GS_", "")
    f = open(os.path.join(args.outputPath, fileName), "w")
    for r in RESULT[0]:
        f.write(str(r)+"\n")
if __name__ == "__main__":
    main()