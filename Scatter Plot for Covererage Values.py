import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
f1path="C:/Users/meeth/Desktop/RP2Files/Canonical data\gnomAD_canonical build37.csv"
data = pd.read_csv(f1path, sep=",", header=0)
data
f2path="C:/Users/meeth/Desktop/RP2Files/MergedGenexomdata/merged_exomecoverage_data.csv"
d1= pd.read_csv(f2path, sep=",", header=0)
d1=d1.drop(["Unnamed: 0","strand"], axis=1)
d1



######################Filtering with canonical data
# Create a list of unique transcript_ids from the first dataframe
# Create a list of unique transcript_ids from the first dataframe

d1 = d1.merge(data[['transcript_id']], on='transcript_id', how='inner')

f3path="C:/Users/meeth/Desktop/RP2Files/Canonical data/Canonicalexomecoverage.csv"
f4path="C:/Users/meeth/Desktop/RP2Files/Canonical data/Canonicalgenomecoverage.csv"
d3=pd.read_csv(f3path, sep=",", header=0)
d4=pd.read_csv(f4path, sep=",", header=0)
d3
d4

#Scatter plot

total_transcript_ids1 = len(d3)  # Assuming both dataframes have the same length
total_transcript_ids2 = len(d4)
# Create subplots for four scatter plots (2x2 grid)
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle("Scatter Plots of Exon and CDS Coverages (Threshold = 20)")

# Define the threshold for low coverage
threshold = 20

# Scatter plot for exon_coverage vs. row index in exome data (labeled as 'a')
axes[0, 0].scatter(np.arange(total_transcript_ids1), d3['exon_coverage'], s=10, color='blue', alpha=0.7)
axes[0, 0].axhline(y=threshold, color='red', linestyle='--', label=f'Threshold ({threshold})', linewidth=1)
axes[0, 0].axvline(x=np.where(d3['exon_coverage'] < threshold)[0], color='red', linestyle='--', alpha=0.5, label=f'Below Threshold', linewidth=1)
axes[0, 0].set_title("a) Exon Coverage (Exome)")
axes[0, 0].set_xlabel("Row Index")
axes[0, 0].set_ylabel("Coverage")
axes[0, 0].legend()

# Scatter plot for CDS_coverage vs. row index in exome data (labeled as 'b')
axes[0, 1].scatter(np.arange(total_transcript_ids1), d3['CDS_coverage'], s=10, color='green', alpha=0.7)
axes[0, 1].axhline(y=threshold, color='red', linestyle='--', label=f'Threshold ({threshold})', linewidth=1)
axes[0, 1].axvline(x=np.where(d3['CDS_coverage'] < threshold)[0], color='red', linestyle='--', alpha=0.5, label=f'Below Threshold', linewidth=1)
axes[0, 1].set_title("b) CDS Coverage (Exome)")
axes[0, 1].set_xlabel("Row Index")
axes[0, 1].set_ylabel("Coverage")
axes[0, 1].legend()

# Scatter plot for exon_coverage vs. row index in genome data (labeled as 'c')
axes[1, 0].scatter(np.arange(total_transcript_ids2), d4['exon_coverage'], s=10, color='red', alpha=0.7)
axes[1, 0].axhline(y=threshold, color='blue', linestyle='--', label=f'Threshold ({threshold})', linewidth=1)
axes[1, 0].axvline(x=np.where(d4['exon_coverage'] < threshold)[0], color='blue', linestyle='--', alpha=0.5, label=f'Below Threshold', linewidth=1)
axes[1, 0].set_title("c) Exon Coverage (Genome)")
axes[1, 0].set_xlabel("Row Index")
axes[1, 0].set_ylabel("Coverage")
axes[1, 0].legend()

# Scatter plot for CDS_coverage vs. row index in genome data (labeled as 'd')
axes[1, 1].scatter(np.arange(total_transcript_ids2), d4['CDS_coverage'], s=10, color='purple', alpha=0.7)
axes[1, 1].axhline(y=threshold, color='blue', linestyle='--', label=f'Threshold ({threshold})', linewidth=1)
axes[1, 1].axvline(x=np.where(d4['CDS_coverage'] < threshold)[20], color='blue', linestyle='--', alpha=0.5, label=f'Below Threshold', linewidth=1)
axes[1, 1].set_title("d) CDS Coverage (Genome)")
axes[1, 1].set_xlabel("Row Index")
axes[1, 1].set_ylabel("Coverage")
axes[1, 1].legend()

# Adjust layout for subplots
plt.tight_layout(rect=[0, 0, 1, 0.95])

# Show the scatter plots
plt.show()


lowcoverexomepath="C:/Users/meeth/Desktop/RP2Files/Canonical data/ExomeCanoncoverdisease.csv"
coverage_df= pd.read_csv(lowcoverexomepath, sep=",", header=0)



coverage_df['CDS_coverage'] = coverage_df['CDS_coverage'].astype(float)
coverage_df['exon_coverage'] = coverage_df['exon_coverage'].astype(float)

# Filter rows where both 'CDS_coverage' and 'exon_coverage' are >= 20
filtered_df = coverage_df[(coverage_df['CDS_coverage'] <= 20.0)]

# Reset the index of the filtered DataFrame
filtered_df.reset_index(drop=True, inplace=True)
filtered_df


d1 = d1.merge(filtered_df[['transcript_id']], on='transcript_id', how='inner')
d1
d1.to_csv("Canonexomedata.csv",sep=",", index=False)

df = pd.DataFrame(d1)
df
