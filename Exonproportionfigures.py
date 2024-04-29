import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt

# Define the file path
file_path = "C:/Users/meeth/Desktop/RP2Files/MergedGenexomdata/merged_exomecoverage_data.csv"  # Replace with your file path

# Read the data into a DataFrame
df = pd.read_csv(file_path, sep=',')  # Assuming tab-delimited file
df=df.drop(["Unnamed: 0","strand"], axis=1)
df



f1 = "C:/Users/meeth/Desktop/RP2Files/MergedGenexomdata/merged_genomecoverage_data.csv"  # Replace with your file path

# Read the data into a DataFrame
df1 = pd.read_csv(f1, sep=',')  # Assuming tab-delimited file
df1=df1.drop(["Unnamed: 0","strand"], axis=1)
df1





# Define the low coverage thresholds for both DataFrames
low_coverage_threshold_df = 20
low_coverage_threshold_df1 = 20

# Create masks for exons and CDS regions below the low coverage thresholds for df
low_coverage_mask_df = (df['mean'] < low_coverage_threshold_df)

# Create masks for exons and CDS regions below the low coverage thresholds for df1
low_coverage_mask_df1 = (df1['mean'] < low_coverage_threshold_df1)

# Calculate the number of exons and CDS regions below the low coverage thresholds for df
exons_low_coverage_df = len(df[(df['biotype'] == 'exon') & low_coverage_mask_df])
cds_low_coverage_df = len(df[(df['biotype'] == 'CDS') & low_coverage_mask_df])

# Calculate the total number of exons for df
total_exons_df = len(df[df['biotype'] == 'exon'])

# Calculate the number of exons and CDS regions below the low coverage thresholds for df1
exons_low_coverage_df1 = len(df1[(df1['biotype'] == 'exon') & low_coverage_mask_df1])
cds_low_coverage_df1 = len(df1[(df1['biotype'] == 'CDS') & low_coverage_mask_df1])

# Calculate the total number of exons for df1
total_exons_df1 = len(df1[df1['biotype'] == 'exon'])

# Create data for the pie charts
labels = ['Exons Below Threshold', 'Exons Above Threshold', 'CDS Below Threshold', 'CDS Above Threshold']
sizes_df = [exons_low_coverage_df, total_exons_df - exons_low_coverage_df, cds_low_coverage_df, len(df[df['biotype'] == 'CDS']) - cds_low_coverage_df]
sizes_df1 = [exons_low_coverage_df1, total_exons_df1 - exons_low_coverage_df1, cds_low_coverage_df1, len(df1[df1['biotype'] == 'CDS']) - cds_low_coverage_df1]
colors = ['red', 'lightcoral', 'blue', 'lightskyblue']
explode = (0.1, 0, 0.1, 0)  # Explode slices for emphasis (optional)

# Create a 2x2 grid for pie charts
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Plot the pie chart for DataFrame df
axes[0].pie(sizes_df, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140, explode=explode)
axes[0].set_title(f'Exome Coverage (Lower Coverage Threshold = {low_coverage_threshold_df})')
axes[0].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

# Plot the pie chart for DataFrame df1
axes[1].pie(sizes_df1, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140, explode=explode)
axes[1].set_title(f'Genome Coverage (Lower Coverage Threshold = {low_coverage_threshold_df1})')
axes[1].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

# Adjust spacing between subplots
plt.tight_layout()

# Show the chart
plt.show()



