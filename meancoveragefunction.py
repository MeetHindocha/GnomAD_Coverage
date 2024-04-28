import dask.dataframe as dd
import pandas as pd
import time
# Load Dask dataframes from CSV files


range_df = pd.read_csv("Filtered_DGV_data.csv",sep=",")
# Convert multiple columns into lists using list comprehensions
start_positions  = range_df['start'].tolist()
end_positions = range_df['end'].tolist()
chrom_list=[21]

print(start_positions)
print(end_positions)

fpath="D:\ExomeandGenomeCoveragedataZip\gnomad.exomes.coverage.summary1.tsv"
df = dd.read_csv("21.tsv",sep="\t")
df = df.compute()


def calculate_mean_coverage(df, chrom_list, start_positions, end_positions):
    """
    Calculate the mean coverage for specified ranges and chromosomes using a DataFrame.

    Args:
    df (pd.DataFrame): The DataFrame containing numeric 'chrom', 'pos', and 'mean' columns.
    chrom_list (list): List of chromosome numbers to consider.
    start_positions (list): List of start positions.
    end_positions (list): List of end positions.

    Returns:
    list: List of mean coverages for each specified range and chromosome.
    """
    mean_coverages = []
    
    for chrom in chrom_list:
        # Filter DataFrame by chromosome
        chrom_df = df[df['chrom'] == chrom]
        
        # Calculate mean coverage for each range
        range_mean_coverages = []
        for start, end in zip(start_positions, end_positions):
            # Filter DataFrame to include only rows within the specified range
            filtered_df = chrom_df[(chrom_df['pos'] >= start) & (chrom_df['pos'] <= end)]
            
            # Calculate the mean coverage for the filtered rows
            mean_coverage = filtered_df['mean'].mean()
            
            range_mean_coverages.append(mean_coverage)
        
        mean_coverages.extend(range_mean_coverages)
    
    return mean_coverages

# Measure execution time
start_time = time.time()

# Calculate the mean coverage for the specified ranges and chromosomes
mean_coverages = calculate_mean_coverage(df, chrom_list, start_positions, end_positions)
mean_coverages
# Measure end time
end_time = time.time()

# Calculate and print the execution time
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")

# Print the results
for i, (chrom, start, end) in enumerate(zip(chrom_list, start_positions, end_positions)):
    print(f"Mean coverage for Range {i+1} (Chromosome {chrom}, {start}-{end}): {mean_coverages[i]}")