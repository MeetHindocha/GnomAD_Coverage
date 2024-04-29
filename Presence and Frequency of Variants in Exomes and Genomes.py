import functools
import pandas as pd
import numpy as np
#Exome Coverage data 
#
f1path="C:/Users/meeth/Desktop/RP2Files/Canonical data/Canonicalexomecoverage.csv"
df1=pd.read_csv(f1path, sep=",", header=0)
df1
df1 = df1.dropna(subset=['exon_coverage', 'CDS_coverage'])

# Convert decimal numbers to the nearest integer value
df1['exon_coverage'] = df1['exon_coverage'].round().astype(int)
df1['CDS_coverage'] = df1['CDS_coverage'].round().astype(int)


# Drop rows with value 0 in both 'exon_coverage' and 'CDS_coverage'
df1 = df1[(df1['exon_coverage'] != 0) & (df1['CDS_coverage'] != 0)]
# Print the updated DataFrame
print(len(df1))
df1

f2path="C:/Users/meeth/Desktop/Scripts/merged_exomecoverage_data.csv"
df2=pd.read_csv(f2path, sep=",", header=0)
df2=df2.drop(["mean","strand"], axis=1)
df2 = df2[df2['biotype'] == 'transcript']
df2


#Creating a merged dataframe
# Merge the start and end columns from df2 into df1 based on transcript_id
merged_df1 = pd.merge(df1, df2[['transcript_id', 'start', 'end']], on='transcript_id', how='left')
merged_df1 = merged_df1.drop(["biotype","gene_id_y",], axis=1)
merged_df1# Calculate nucleotide count for each row and add it as a new column
merged_df1['nucleotide_count'] = merged_df1['end'] - merged_df1['start'] + 1
merged_df1


variantpath="C:/Users/meeth/Desktop/Scripts/FilteredDGVdata.csv"
variantdf=pd.read_csv(variantpath, sep=",", header=0)
variantdf



def count_variants(row):
    start = row['start']
    end = row['end']
    variant_count = variantdf[(variantdf['start'] >= start) & (variantdf['end'] <= end)].shape[0]
    return variant_count

# Apply the function to calculate variant count and add it as a new column
merged_df1['variant_count'] = merged_df1.apply(count_variants, axis=1)

# Print the updated DataFrame
print(merged_df1)



#################################
# Create a list of coverage values from 1 to 100
coverage_values = list(range(1, 101))

# Initialize empty lists to store the number of variants for each coverage value
variant_counts_exon = []
variant_counts_CDS = []

# Calculate the number of variants for each coverage value in 'exon_coverage' and 'CDS_coverage'
for coverage in coverage_values:
    variants_at_coverage_exon = merged_df1[merged_df1['exon_coverage'] == coverage]['variant_count'].sum()
    variant_counts_exon.append(variants_at_coverage_exon)
    
    variants_at_coverage_CDS = merged_df1[merged_df1['CDS_coverage'] == coverage]['variant_count'].sum()
    variant_counts_CDS.append(variants_at_coverage_CDS)

# Create a new DataFrame with coverage values, variant counts for 'exon_coverage', and variant counts for 'CDS_coverage'
new_df = pd.DataFrame({
    'Coverage': coverage_values,
    'Variant_Count_Exon': variant_counts_exon,
    'Variant_Count_CDS': variant_counts_CDS
})

# Print the new DataFrame
print(new_df)



#####################################################################
# Create a list of coverage values from 1 to 100
coverage_values = list(range(1, 101))

# Initialize empty lists to store the number of nucleotides for each coverage value in 'exon_coverage' and 'CDS_coverage'
nucleotide_counts_exon = []
nucleotide_counts_CDS = []

# Calculate the total number of nucleotides for each coverage value in 'exon_coverage' and 'CDS_coverage'
for coverage in coverage_values:
    nucleotides_at_coverage_exon = merged_df1[merged_df1['exon_coverage'] == coverage]['nucleotide_count'].sum()
    nucleotide_counts_exon.append(nucleotides_at_coverage_exon)
    
    nucleotides_at_coverage_CDS = merged_df1[merged_df1['CDS_coverage'] == coverage]['nucleotide_count'].sum()
    nucleotide_counts_CDS.append(nucleotides_at_coverage_CDS)

# Add new columns to the existing new_df DataFrame
new_df['Nucleotide_Count_Exon'] = nucleotide_counts_exon
new_df['Nucleotide_Count_CDS'] = nucleotide_counts_CDS

# Print the updated new_df DataFrame
print(new_df)
Exome_df = pd.DataFrame(new_df)
new_df.to_csv("Avg_Variants_Exomecoverage.csv",sep=",", index=False)

new_df
Exome_df=pd.DataFrame(new_df)
Exome_df
##############################################################################
import pandas as pd
import matplotlib.pyplot as plt

# Assuming you have a DataFrame named 'new_df' with 'Coverage', 'Variant_Count_Exon', and 'Variant_Count_CDS' columns

# Calculate the average number of variants per coverage value
new_df['Avg_Variants_Exon'] = new_df['Variant_Count_Exon'] / new_df['Nucleotide_Count_Exon']
new_df['Avg_Variants_CDS'] = new_df['Variant_Count_CDS'] / new_df['Nucleotide_Count_CDS']

# Plotting the figures
plt.figure(figsize=(12, 6))

# Plot for Exon data
plt.subplot(1, 2, 1)
plt.bar(new_df['Coverage'], new_df['Avg_Variants_Exon'], color='blue')
plt.xlabel('Coverage')
plt.ylabel('Average Variants per Nucleotide')
plt.title('Exon Data')
plt.grid()

# Plot for CDS data
plt.subplot(1, 2, 2)
plt.bar(new_df['Coverage'], new_df['Avg_Variants_CDS'], color='green')
plt.xlabel('Coverage')
plt.ylabel('Average Variants per Nucleotide')
plt.title('CDS Data')
plt.grid()

plt.tight_layout()
plt.show()

#########################################################################################################
#Genome Coverage data 


f1path="C:/Users/meeth/Desktop/RP2Files/Canonical data/Canonicalgenomecoverage.csv"
df1=pd.read_csv(f1path, sep=",", header=0)
df1
df1 = df1.dropna(subset=['exon_coverage', 'CDS_coverage'])
print(len(df1))
# Convert decimal numbers to the nearest integer value
df1['exon_coverage'] = df1['exon_coverage'].round().astype(int)
df1['CDS_coverage'] = df1['CDS_coverage'].round().astype(int)


# Drop rows with value 0 in both 'exon_coverage' and 'CDS_coverage'
df1 = df1[(df1['exon_coverage'] != 0) & (df1['CDS_coverage'] != 0)]
# Print the updated DataFrame
print(len(df1))
df1

f2path="C:/Users/meeth/Desktop/Scripts/merged_genomecoverage_data.csv"
df2=pd.read_csv(f2path, sep=",", header=0)
df2=df2.drop(["Unnamed: 0","strand","mean","exon_id"], axis=1)
df2 = df2[df2['biotype'] == 'transcript']
df2


#Creating a merged dataframe
# Merge the start and end columns from df2 into df1 based on transcript_id
merged_df = pd.merge(df1, df2[['transcript_id', 'start', 'end']], on='transcript_id', how='left')
merged_df = merged_df.drop(["biotype","gene_id_y"], axis=1)
merged_df# Calculate nucleotide count for each row and add it as a new column
merged_df['nucleotide_count'] = merged_df['end'] - merged_df['start'] + 1
merged_df
merged_df

variantpath="C:/Users/meeth/Desktop/Scripts/Filtered_DGV_data.csv"
variantdf=pd.read_csv(variantpath, sep=",", header=0)
variantdf




def count_variants(row):
    start = row['start']
    end = row['end']
    variant_count = variantdf[(variantdf['start'] >= start) & (variantdf['end'] <= end)].shape[0]
    return variant_count

# Apply the function to calculate variant count and add it as a new column
merged_df['variant_count'] = merged_df.apply(count_variants, axis=1)

# Print the updated DataFrame
print(merged_df)



#################################
coverage_values = list(range(1, 101))

# Initialize empty lists to store the number of variants for each coverage value
variant_counts_exon = []
variant_counts_CDS = []
# Calculate the number of variants for each coverage value in 'exon_coverage' and 'CDS_coverage'
for coverage in coverage_values:
    variants_at_coverage_exon = merged_df[merged_df['exon_coverage'] == coverage]['variant_count'].sum()
    variant_counts_exon.append(variants_at_coverage_exon)
    
    variants_at_coverage_CDS = merged_df[merged_df['CDS_coverage'] == coverage]['variant_count'].sum()
    variant_counts_CDS.append(variants_at_coverage_CDS)

# Ensure that both lists have the same length
while len(variant_counts_exon) < len(coverage_values):
    variant_counts_exon.append(0)

while len(variant_counts_CDS) < len(coverage_values):
    variant_counts_CDS.append(0)

# Create a new DataFrame with coverage values, variant counts for 'exon_coverage', and variant counts for 'CDS_coverage'
new_df = pd.DataFrame({
    'Coverage': coverage_values,
    'Variant_Count_Exon': variant_counts_exon,
    'Variant_Count_CDS': variant_counts_CDS
})

# Print the new DataFrame
print(new_df)
new_df.to_csv("file.csv",sep=",", index=False)
Genome_df=pd.DataFrame(new_df)
Exome_df
Genome_df
#####################################################################
# Create a list of coverage values from 1 to 100
coverage_values = list(range(1, 101))

# Initialize empty lists to store the number of nucleotides for each coverage value in 'exon_coverage' and 'CDS_coverage'
nucleotide_counts_exon = []
nucleotide_counts_CDS = []

# Calculate the total number of nucleotides for each coverage value in 'exon_coverage' and 'CDS_coverage'
for coverage in coverage_values:
    nucleotides_at_coverage_exon = merged_df[merged_df['exon_coverage'] == coverage]['nucleotide_count'].sum()
    nucleotide_counts_exon.append(nucleotides_at_coverage_exon)
    
    nucleotides_at_coverage_CDS = merged_df[merged_df['CDS_coverage'] == coverage]['nucleotide_count'].sum()
    nucleotide_counts_CDS.append(nucleotides_at_coverage_CDS)

# Add new columns to the existing new_df DataFrame
new_df['Nucleotide_Count_Exon'] = nucleotide_counts_exon
new_df['Nucleotide_Count_CDS'] = nucleotide_counts_CDS

# Print the updated new_df DataFrame
print(new_df)


# Calculate the average number of variants per coverage value
new_df['Avg_Variants_Exon'] = new_df['Variant_Count_Exon'] / new_df['Nucleotide_Count_Exon']
new_df['Avg_Variants_CDS'] = new_df['Variant_Count_CDS'] / new_df['Nucleotide_Count_CDS']





import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Assuming you have two dataframes: Exome_df and Genome_df
# Extract the relevant columns for both dataframes
exome_avg_variants_exon = Exome_df['Avg_Variants_Exon']
exome_avg_variants_cds = Exome_df['Avg_Variants_CDS']

genome_avg_variants_exon = Genome_df['Avg_Variants_Exon']
genome_avg_variants_cds = Genome_df['Avg_Variants_CDS']

# Create a figure with subplots
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10), sharex=True)

# Plot 1: Avg_Variants_Exon for Exome_df
sns.barplot(x=np.arange(1, 101), y=exome_avg_variants_exon, ax=axes[0, 0])
axes[0, 0].set_title('a) Avg_Variants_Exon (Exome)')
axes[0, 0].set_ylabel('Average Variants')

# Plot 2: Avg_Variants_CDS for Exome_df
sns.barplot(x=np.arange(1, 101), y=exome_avg_variants_cds, ax=axes[0, 1])
axes[0, 1].set_title('b) Avg_Variants_CDS (Exome)')
axes[0, 1].set_ylabel('Average Variants')

# Plot 3: Avg_Variants_Exon for Genome_df
sns.barplot(x=np.arange(1, 101), y=genome_avg_variants_exon, ax=axes[1, 0])
axes[1, 0].set_title('c) Avg_Variants_Exon (Genome)')
axes[1, 0].set_xlabel('Coverage')
axes[1, 0].set_ylabel('Average Variants')

# Plot 4: Avg_Variants_CDS for Genome_df
sns.barplot(x=np.arange(1, 101), y=genome_avg_variants_cds, ax=axes[1, 1])
axes[1, 1].set_title('d) Avg_Variants_CDS (Genome)')
axes[1, 1].set_xlabel('Coverage')
axes[1, 1].set_ylabel('Average Variants')


# Adjust x-axis tick labels to show only '1'
for ax in axes.flat:
    ax.set_xticks([0,10,20,30,40,50,60,70,80,90,99])

# Adjust layout
plt.tight_layout()
plt.show()
