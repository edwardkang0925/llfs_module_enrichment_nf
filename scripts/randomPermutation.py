import pandas as pd
import argparse
import os

"""
Given input_file_path, output_directory, columnToPermute, and seed, this function will:
generate 1000 RP files in the output_directory. The RP files will have the same name as the input_file_path with the addition of the seed number. (ABI-1.csv, ABI-2.csv, ABI-3.csv, etc.)
Seeds are integers from 1 to 1000.

Usage:
python3 randomPermutation.py input_file_path output_directory columnToPermute 
"""

def permute_first_column(input_file_path, output_directory, columnToPermute, seed=None):
    # Step 1: Read the CSV file into a DataFrame
    df = pd.read_csv(input_file_path)
    
    # Step 2: Randomly permute the values in the specified column
    df[columnToPermute] = df[columnToPermute].sample(frac=1, random_state=seed).values
    
    # Step 3: Drop the 'Unnamed: 0' column if it exists
    if 'Unnamed: 0' in df.columns:
        df = df.drop(columns=['Unnamed: 0'])
    
    # Step 4: Save the modified DataFrame to a new CSV file
    # Extract the base file name and append "_{seed}.csv"
    base_name = os.path.splitext(os.path.basename(input_file_path))[0]
    output_file_name = f"{seed}-{base_name}.csv"
    output_file_path = os.path.join(output_directory, output_file_name)
    df.to_csv(output_file_path, index=False)

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Randomly permute the first column of a CSV file.")
    
    # Add arguments to parser
    parser.add_argument("input_file_path", help="Path to the input CSV file.")
    parser.add_argument("output_directory", help="Directory to save the permuted CSV files.")
    parser.add_argument("column_name_to_permute", help="Name of the column to permute.")
    parser.add_argument("numRP", type=int, help="number of permutations")

    
    # Parse the arguments
    args = parser.parse_args()
    
    # Check if the output directory exists, if not create it
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    
    # Loop to generate 1000 permuted DataFrames
    for i in range(1, args.numRP+1):
        # Set seed for reproducibility
        seed = i
        # Call the function
        permute_first_column(args.input_file_path, args.output_directory, args.column_name_to_permute,seed)

if __name__ == "__main__":
    main()
