import pandas as pd

# Path to the first CSV file containing gnomAD data
file_path = "C:/Users/meeth/Desktop/RP2()/gnomAD_canonical build37.csv"

# Read the first CSV file into a DataFrame
df1 = pd.read_csv(file_path, sep=',', header=0)

# Path to the second CSV file containing genome coverage data for canonical transcripts
fpath = "C:/Users/meeth/Desktop/Genome Coverage/CanonicalTranscriptsGcover.csv"

# Read the second CSV file into a DataFrame
df2 = pd.read_csv(fpath, sep=',', header=0)

# Merge the two DataFrames based on 'transcript_id', keeping all rows from the first DataFrame
# and adding matching columns from the second DataFrame
merged_df = df1.merge(df2, on='transcript_id', how='left')

# Path to the output file where the merged DataFrame will be saved
output_file = "C:/Users/meeth/Desktop/Exome Coverage/Canonicalgenoemcoverage.csv"

# Write the merged DataFrame to a new CSV file without including the index
merged_df.to_csv(output_file, index=False)




