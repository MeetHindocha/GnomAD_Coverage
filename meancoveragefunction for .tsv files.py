import os
import pandas as pd
import gzip

# Directory containing .gz files
gz_folder = "C:/Users/meeth/Desktop/Coverage data/ExomeCoverageData"

# Initialize an empty list to store DataFrames
dataframes = []

# List all the .gz files in the folder
gz_files = [f for f in os.listdir(gz_folder) if f.endswith('.gz')]

# Loop through each .gz file and extract data
for gz_file in gz_files:
    with gzip.open(os.path.join(gz_folder, gz_file), 'rt') as f:
        # Read the content of the .gz file into a DataFrame
        data_df = pd.read_csv(f, sep='\t')

        # Append the DataFrame to the list
        dataframes.append(data_df)

# Concatenate the list of DataFrames into a single DataFrame
master_df = pd.concat(dataframes, ignore_index=True)

# Save the merged DataFrame to an output file



##############################
#Filtering canonical data
################################

fpath="C:/Users/meeth/Desktop/RP2Files/Canonical data/gnomAD_canonical build37.csv"
canondata = pd.read_csv(fpath, sep=",", header=0)



# Merge the DataFrames to filter rows in df1
master_df = master_df.merge(canondata, on='transcript_id')
merged_df=master_df.drop(["gene_id_y","gene_name"], axis=1)
output_file = 'merged_exomecoverage_data.csv'
master_df.to_csv(output_file, index=True)
print(f"Merged data saved to '{output_file}'.")


file_path = "C:/Users/meeth/Desktop/Scripts/filteredDGVdata.csv"
DGV_data = pd.read_csv(file_path, sep=",", header=0)

from intervaltree import IntervalTree

# Assuming you have DGV_data and merged_df DataFrames

# Filter rows in merged_df for biotypes 'exon' and 'CDS'
exon_cds_rows = merged_df[merged_df['biotype'].isin(['exon', 'CDS'])]

# Create an interval tree for exon and CDS ranges
interval_tree = IntervalTree()
for idx, row in exon_cds_rows.iterrows():
    interval_tree[row['start']:row['end']] = idx  # Store the index for reference

# Filter DGV_data based on interval tree
filtered_rows = []

for idx, row in DGV_data.iterrows():
    overlapping_intervals = interval_tree[row['start']:row['end']]
    if overlapping_intervals:
        # Check if any overlapping interval belongs to exon or CDS
        overlapping_biotypes = set(exon_cds_rows.loc[iv.data, 'biotype'] for iv in overlapping_intervals)
        if any(biotype in overlapping_biotypes for biotype in ['exon', 'CDS']):
            filtered_rows.append(row)

# Create a new DataFrame from filtered rows
filtered_DGV_data = pd.DataFrame(filtered_rows)
d6=filtered_DGV_data
d6=d6.drop(["variantaccession","varianttype","variantsubtype","genes"], axis=1)

# Print the filtered DGV_data DataFrame
print(len(filtered_DGV_data))
filtered_DGV_data.to_csv("Filtered_DGV_data.csv",sep=",", index=False)
# Define the path to the large CSV file
csv_file_path = "C:/Users/meeth/Desktop/coverage/gnomad.exomes.coverage.summary1.tsv"
import dask.dataframe as dd
# Define data types for columns (adjust as needed)
dtype_mapping = {'chrom': str, 'pos': int, 'mean': float}

csv_ddf = dd.read_csv(csv_file_path, dtype=dtype_mapping)
#changing name of column
d6 = d6.rename(columns={'chr': 'chrom'})

dask_ddf = csv_ddf.compute().reset_index(drop=True)

# Function to calculate mean coverage using Dask
def calculate_mean_coverage_dask(row):
    coverage_slice = dask_ddf[
        (dask_ddf['chrom'] == row['chrom']) & 
        (dask_ddf['pos'] >= row['start']) & 
        (dask_ddf['pos'] <= row['end'])
    ]
    return coverage_slice['mean'].mean().compute()

# Apply the function to the Pandas DataFrame
d6['mean_coverage'] = d6.apply(calculate_mean_coverage_dask, axis=1)

# Print the updated Pandas DataFrame
print(d6)