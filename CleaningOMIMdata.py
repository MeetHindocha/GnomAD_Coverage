import pandas as pd

# Disable chained assignment warning
pd.options.mode.chained_assignment = None  # default='warn'

# Read data from the first CSV file and drop the 'Entrez Gene ID (NCBI)' column
file_path = "C:/Users/meeth/Desktop/OMIM data/Geneidentifierofdgenesexcel.csv"
data = pd.read_csv(file_path, sep=",", header=0)
d = data.drop(["Entrez Gene ID (NCBI)"], axis=1)

# Save the modified DataFrame to a new CSV file
d.to_csv("GeneIdentifier.csv", sep=",", index=False)

# Filter rows where 'MIM Entry Type' is 'gene'
d = d[d['MIM Entry Type (see FAQ 1.3 at https://omim.org/help/faq)'] == 'gene']

# Display the filtered DataFrame
print(d)

# Read data from three CSV files
f1path = "C:/Users/meeth/Desktop/RP2()/OMIMDiseases1.csv"
f2path = "C:/Users/meeth/Desktop/RP2()/OMIMDiseases2.csv"
f3path = "C:/Users/meeth/Desktop/RP2()/OMIMDiseases3.csv"
d1 = pd.read_csv(f1path, sep=",", header=0)
d2 = pd.read_csv(f2path, sep=",", header=0)
d3 = pd.read_csv(f3path, sep=",", header=0)

# Merge the three DataFrames on 'MIM Number'
merged_df = pd.merge(d1, d2, on='MIM Number', how='outer')
merged_df = pd.merge(merged_df, d3, on='MIM Number', how='outer')

# Display the merged DataFrame
print(merged_df)

# Save the merged DataFrame to a new CSV file
merged_df.to_csv("OMIM.csv", sep=",", index=False)

# Read data from the fourth CSV file
f4path = "C:/Users/meeth/Desktop/RP2()/GeneIdentifier.csv"
d4 = pd.read_csv(f4path, sep=",", header=0)

# Drop rows with NaN values in the 'Approved Gene Symbol' column
merged_df = merged_df.dropna(subset=['Approved Gene Symbol'])

# Merge the fourth DataFrame with the filtered 'merged_df' DataFrame on 'Approved Gene Symbol'
merged_df1 = pd.merge(d4, merged_df, on='Approved Gene Symbol', how='inner')

# Drop unnecessary columns from the merged DataFrame
merged_df1 = merged_df1.drop(["MIM Entry Type (see FAQ 1.3 at https://omim.org/help/faq)",
                              "Ensembl Gene ID (Ensembl)", "Gene Symbols_x", "Entrez Gene ID"], axis=1)

# Save the final merged DataFrame to a new CSV file
merged_df1.to_csv("FinalOMIM.csv", sep=",", index=False)

# Read data from the fifth and sixth CSV files
f6path = "C:/Users/meeth/Desktop/RP2()/FinalOMIM.csv"
f7path = "C:/Users/meeth/Desktop/RP2Files/Canonicalexomecoverage.csv"
d6 = pd.read_csv(f6path, sep=",", header=0)
d7 = pd.read_csv(f7path, sep=",", header=0)

# Merge the fifth and sixth DataFrames on 'gene_name'
merged_df2 = pd.merge(d7, d6, on='gene_name', how='left')

# Save the final merged DataFrame to a new CSV file
merged_df2.to_csv("ExomeCanoncoverdisease.csv", sep=",", index=False)
