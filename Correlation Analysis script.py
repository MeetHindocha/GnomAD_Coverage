import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load data from the file
data = {
    'Gene_name': ['RPS17', 'SMN2', 'OPN1MW', 'GRAP', 'C4A', 'C4B', 'IKBKG', 'CCL3L1', 'CFC1', 'SMN1'],
    'Exome_Exon_Coverage_data': [2, 4, 14, 11, 9, 8, 13, 7, 17, 4],
    'Exome_CDS_Coverage': [3, 5, 12, 11, 12, 12, 13, 13, 17, 6],
    'Genome_Exon_Coverage': [1, 3, 6, 13, 5, 3, 9, 12, 11, 8],
    'Genome_CDS_Coverage': [1, 3, 5, 13, 4, 5, 7, 12, 12, 8],
    'GeVIR_score': [90.93, 70.56, 77.51, 29.51, 52.04, 53.67, 75.03, 89.93, 88.36, 71.22],
    'Transcript_stable_ID': ['ENST00000330339', 'ENST00000380743', 'ENST00000369929', 'ENST00000284154',
                             'ENST00000428956', 'ENST00000435363', 'ENST00000369609', 'ENST00000422211',
                             'ENST00000259216', 'ENST00000380707']
}

df = pd.DataFrame(data)

# Calculate correlation coefficients
corr_matrix = df[['Exome_Exon_Coverage_data', 'Exome_CDS_Coverage', 'Genome_Exon_Coverage', 
                  'Genome_CDS_Coverage', 'GeVIR_score']].corr()

# Plotting the correlation heatmap
plt.figure(figsize=(8, 6))
heatmap = sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1, fmt=".2f")

# Rotate x-axis tick labels
heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')

plt.title('Correlation Heatmap')
plt.tight_layout()
plt.show()
