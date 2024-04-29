import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
file_path = "C:/Users/meeth/Desktop/RP2Files/Canonical data/ExomeCanoncoverdisease.csv"
data = pd.read_csv(file_path, sep=",", header=0)
data
# Assuming you have a DataFrame called 'df' with 'exon_coverage' and 'CDS_coverage' columns

# Define the threshold
threshold = 20

# Create a boolean mask to identify rows where 'CDS_coverage' is above the threshold
mask = data['CDS_coverage'] > threshold

# Use the mask to drop rows where 'CDS_coverage' is above the threshold
Lowcovexomecandis = data[~mask]
Lowcovexomecandis = Lowcovexomecandis.drop_duplicates(subset='gene_name', keep='first')
Lowcovexomecandis.to_csv("Lowcovexomecandis.csv",sep=",", index=False)
print(len(Lowcovexomecandis))







f1 = "C:/Users/meeth/Desktop/RP2Files/Canonical data/GeneCanoncoverdisease.csv"
data1 = pd.read_csv(f1, sep=",", header=0)
print((data1))
# Assuming you have a DataFrame called 'df' with 'exon_coverage' and 'CDS_coverage' columns

# Define the threshold
threshold = 20

# Create a boolean mask to identify rows where 'CDS_coverage' is above the threshold
mask = data1['CDS_coverage'] > threshold

# Use the mask to drop rows where 'CDS_coverage' is above the threshold
Lowcovgencandis = data1[~mask]
Lowcovgencandis = Lowcovgencandis.drop_duplicates(subset='gene_name', keep='first')


Lowcovgencandis.to_csv("Lowcovgencandis.csv",sep=",", index=False)
print(len(Lowcovgencandis))





f2="C:/Users/meeth/Desktop/RP2Files/Gevirscores.csv"
d2= pd.read_csv(f2, sep=",", header=0)
d2





# Merge exome_df and genome_df on their respective gene_name columns for exome data
merged_exome_df = pd.merge(Lowcovexomecandis, d2, on='gene_name', how='inner')
merged_exome_df
# Merge exome_df and genome_df on their respective gene_name columns for genome data
merged_genome_df = pd.merge(Lowcovgencandis, d2, on='gene_name', how='inner')
merged_genome_df
# Calculate the correlations
exome_exon_correlation = merged_exome_df['exon_coverage'].corr(merged_exome_df['GeVIR %'])
exome_cds_correlation = merged_exome_df['CDS_coverage'].corr(merged_exome_df['GeVIR %'])
genome_exon_correlation = merged_genome_df['exon_coverage'].corr(merged_genome_df['GeVIR %'])
genome_cds_correlation = merged_genome_df['CDS_coverage'].corr(merged_genome_df['GeVIR %'])

# Print the correlation scores
print(f'Exome Exon Coverage vs. GeVIR % Correlation: {exome_exon_correlation:.2f}')
print(f'Exome CDS Coverage vs. GeVIR % Correlation: {exome_cds_correlation:.2f}')
print(f'Genome Exon Coverage vs. GeVIR % Correlation: {genome_exon_correlation:.2f}')
print(f'Genome CDS Coverage vs. GeVIR % Correlation: {genome_cds_correlation:.2f}')








# Create a 2x2 grid of subplots with shared x-axes
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10), sharex='col')
fig.subplots_adjust(hspace=0.6)  # Adjust the vertical spacing between subplots

# Plot 1: Exon_coverage vs GeVIR% for merged_exome_df
sns.scatterplot(data=merged_exome_df, x='exon_coverage', y='GeVIR %', ax=axes[0, 0])
axes[0, 0].set_title('a) Exon Coverage vs GeVIR% (Exome)')

# Plot 2: CDS_coverage vs GeVIR% for merged_exome_df
sns.scatterplot(data=merged_exome_df, x='CDS_coverage', y='GeVIR %', ax=axes[0, 1])
axes[0, 1].set_title('b) CDS Coverage vs GeVIR% (Exome)')

# Plot 3: Exon_coverage vs GeVIR% for merged_genome_df
sns.scatterplot(data=merged_genome_df, x='exon_coverage', y='GeVIR %', ax=axes[1, 0])
axes[1, 0].set_title('c) Exon Coverage vs GeVIR% (Genome)')

# Plot 4: CDS_coverage vs GeVIR% for merged_genome_df
sns.scatterplot(data=merged_genome_df, x='CDS_coverage', y='GeVIR %', ax=axes[1, 1])
axes[1, 1].set_title('d) CDS Coverage vs GeVIR% (Genome)')

# Set individual x-axis labels for the lower row
axes[1, 0].set_xlabel('Exon Coverage')
axes[1, 1].set_xlabel('CDS Coverage')

# Show the plots
plt.tight_layout()
plt.show()

