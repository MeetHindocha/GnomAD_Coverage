from scipy.stats import fisher_exact
import pandas as pd

import pandas as pd

# Assuming you have loaded your two dataframes, wes_data and wgs_data, which contain 'gene_name', 'exon_coverage', and 'CDS_coverage' columns
f1path= "C:/Users/meeth/Desktop/Scripts/Lowcovexomecandis1.csv"
wes_data= pd.read_csv(f1path, sep=",", header=0)
f2path="C:/Users/meeth/Desktop/Scripts/Lowcovgencandis1.csv"
wgs_data=pd.read_csv(f2path, sep=",", header=0)
# Create a new column to indicate low coverage in each dataset
wes_data['low_coverage_wes'] = True
wgs_data['low_coverage_wgs'] = True

# Merge the two dataframes based on 'gene_name' to align the data
merged_data = pd.merge(wes_data[['gene_name', 'low_coverage_wes']], wgs_data[['gene_name', 'low_coverage_wgs']], on='gene_name', how='outer')

# Fill missing values with False (genes that don't exist in both datasets)
merged_data = merged_data.fillna(False)

# Create a contingency table
contingency_table = pd.crosstab(merged_data['low_coverage_wgs'], merged_data['low_coverage_wes'])

# Print the contingency table
print(contingency_table)



# Perform Fisher's exact test
odds_ratio, p_value = fisher_exact(contingency_table)

# Print the results
print("Odds Ratio:", odds_ratio)
print("P-Value:", p_value)

# Interpret the results based on the p-value
alpha = 0.05  # Significance level
if p_value < alpha:
    print("Reject null hypothesis: There is a significant association between low coverage in exomes and low coverage in WGS.")
else:
    print("Fail to reject null hypothesis: There is no significant association between low coverage in exomes and low coverage in WGS.")
    
    
    
    
    
    
import pandas as pd
import statsmodels.api as sm
wes_data['Phenotypes'].fillna(0, inplace=True)
df=pd.DataFrame(wgs_data)
from scipy.stats import ttest_ind
import pandas as pd
from scipy.stats import ttest_ind


import pandas as pd
from scipy.stats import ttest_ind

# Load your dataset into a DataFrame (assuming it's named df)
# Make sure your dataframe includes 'exon_coverage' and 'CDS_coverage' columns

# Perform a t-test for 'exon_coverage'
t_stat_exon, p_value_exon = ttest_ind(df['exon_coverage'], df['CDS_coverage'], equal_var=False)

# Print the results
print(f"T-Test for Exon Coverage vs. CDS Coverage: t-statistic = {t_stat_exon}, p-value = {p_value_exon}")







