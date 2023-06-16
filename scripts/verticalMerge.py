import pandas as pd
import argparse

def concatenate_csv(file_paths, outputFileName):
    # Read the CSV files and store them in a list
    data_frames = [pd.read_csv(file_path) for file_path in file_paths]
    
    # Concatenate the data frames vertically
    concatenated_df = pd.concat(data_frames, ignore_index=True)
    
    # Save the concatenated data frame to an output CSV file
    concatenated_df.to_csv(outputFileName, index=False)

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description='Concatenate CSV files vertically.')
    parser.add_argument('file_paths', metavar='N', type=str, nargs='+',
                        help='paths to the CSV files to concatenate')
    args = parser.parse_args()
    print(args.file_paths)
    suffix = args.file_paths[0].split("_")[1].split("-")[1]
    outputFileName = f"master_summary_{suffix}_RP.csv"
    # Call the concatenate_csv function with the provided file paths
    concatenate_csv(args.file_paths, outputFileName)
