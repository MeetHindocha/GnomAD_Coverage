import os
import pandas as pd

# Function to calculate coverage for each transcript
def calculate_transcript_coverage(df):
    # Filter rows where 'biotype' column equals 'transcript'
    transcript_mask = df['biotype'] == 'transcript'
    transcripts_df = df[transcript_mask]
    transcripts_df.index = range(len(transcripts_df))

    # Function to calculate coverage for a single transcript
    def calculate_coverage(row):
        transcript_start = row['start']
        transcript_end = row['end']
        # Filter rows outside of the current transcript's range and calculate mean coverage
        transcript_coverage = df.loc[~transcript_mask & (df['start'] >= transcript_start) & (df['end'] <= transcript_end), 'mean']
        return transcript_coverage.mean()

    # Apply the calculate_coverage function to each row and store the result in a new column
    transcripts_df['transcript_coverage'] = transcripts_df.apply(calculate_coverage, axis=1)

    # Return a subset of the DataFrame containing only relevant columns
    return transcripts_df[['biotype', 'start', 'end', 'transcript_coverage',"gene_id","transcript_id"]]

# Folder containing input files
folder_path = "C:/Users/meeth/Desktop/RP2()/GenomeCoverageData"

# List files in the folder
file_list = os.listdir(folder_path)

# Loop through each file in the folder
for file in file_list:
    # Check if the file ends with '.tsv.gz'
    if file.endswith('.tsv.gz'):
        file_path = os.path.join(folder_path, file)

        # Read the compressed CSV file into a DataFrame
        df1 = pd.read_csv(file_path, sep='\t', compression='gzip')
        df = pd.DataFrame(df1)

        # Calculate transcript coverage using the provided function
        result_df = calculate_transcript_coverage(df)

        # Define output file path with '_transcriptcoverage.csv' appended to the original filename
        output_file_path = os.path.join(folder_path, f"{os.path.splitext(file)[0]}_transcriptcoverage.csv")
        
        # Write the result DataFrame to a new CSV file
        result_df.to_csv(output_file_path, index=False)

        print(f"Processed file: {file}")

    else:
        # Skip files that do not end with '.tsv.gz'
        print(f"Skipping non-CSV file: {file}")
