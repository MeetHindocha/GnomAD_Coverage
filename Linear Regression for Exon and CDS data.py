from scipy.stats import fisher_exact
import pandas as pd

f1path= "C:/Users/meeth/Desktop/Scripts/Lowcovexomecandis1.csv"
df1= pd.read_csv(f1path, sep=",", header=0)
f2path= "C:/Users/meeth/Desktop/RP2Files/MergedGenexomdata/merged_exomecoverage_data.csv"
df2=pd.read_csv(f2path, sep=",", header=0)
df1


# Filter rows from df2 where 'biotype' is 'transcript'
filtered_df2 = df2[df2['biotype'] == 'transcript']

# Merge the filtered_df2 with df1 based on 'transcript_id'
df1 = df1.merge(filtered_df2[['transcript_id', 'Chr']], on='transcript_id', how='left')
df1.to_csv("df1.csv",sep=",", index=False)
df1


f3path="C:/Users/meeth/Desktop/Scripts/Lowcovgencandis1.csv"
f4path="C:/Users/meeth/Desktop/RP2Files/MergedGenexomdata/merged_genomecoverage_data.csv"
df3= pd.read_csv(f3path, sep=",", header=0)
df4= pd.read_csv(f4path, sep=",", header=0)


# Filter rows from df2 where 'biotype' is 'transcript'
filtered_df4 = df4[df4['biotype'] == 'transcript']

# Merge the filtered_df2 with df1 based on 'transcript_id'
df3 = df3.merge(filtered_df4[['transcript_id', 'Chr']], on='transcript_id', how='left')
df3.to_csv("df3.csv",sep=",", index=False)



list=df3["Chr"]




import pandas as pd
import statsmodels.api as sm

# Load your dataset into a DataFrame (assuming it's named df)
# Make sure your dataframe includes 'gene_name', 'exon_coverage', 'CDS_coverage', and 'Chr' columns

# Create dummy variables for the 'Chr' column (categorical predictor)
df1 = pd.get_dummies(df1, columns=['Chr'], drop_first=True)

# Define your predictors (independent variables)
# Here, you include the dummy variables for 'Chr'
X = df1[[list]]  # Add more dummy variables as needed

# Add a constant term to the predictors for the intercept
X = sm.add_constant(X)

# Define the dependent variables (outcomes)
Y_exon = df1['exon_coverage']
Y_cds = df1['CDS_coverage']

# Fit linear regression models
model_exon = sm.OLS(Y_exon, X).fit()
model_cds = sm.OLS(Y_cds, X).fit()

# Print model summaries
print("Linear Regression for Exon Coverage:")
print(model_exon.summary())

print("\nLinear Regression for CDS Coverage:")
print(model_cds.summary())