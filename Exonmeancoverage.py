import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
#The first dataframe
file_path = "C:/Users/meeth/Desktop/RP2()/gnomad.exomes.coverage.summary.tsv.bgz"
data = pd.read_csv(file_path, sep='\t', compression='gzip',chunksize=6008666,usecols=["mean","pos","chrom"])
dataframe1=data.get_chunk(6008666)



#The second dataframe
fpath="C:/Users/meeth/Desktop/RP2()/hg19.ensGene.gtf"
df=pd.read_csv(fpath, sep="\t", header=0, chunksize=235052)
df.columns =['Chromosome name',"Transcript", 'Exon', 'SP',"EP","Gene Information","1","2","3"]

df1=df.get_chunk(235052)
df1.columns =['Chromosome name', 'Transcript',"Exon", 'SP', 'EP',"Gene Information","1","2","3"]
df2 = df1[df1["Exon"].str.contains("transcript") == False]
dataframe2 = df2.drop(df2.columns[[1,2,5,6,7]],axis = 1)



# Function to get the average value of positions from dataframe2
def get_average_value(start, end, df2):
    positions = df2[(df2['pos'] >= start) & (df2['pos'] <= end)]
    return positions['mean'].mean()

# Calculate average values for each row in dataframe1
average_values = df2.apply(lambda row: get_average_value(row['SP'], row['EP'], d1), axis=1)

# Add the average_values to the first dataframe as a new column
df2['average_value'] = average_values

df2.to_csv('Exoncoverage.csv')
f = pd.read_csv("file_name3.csv")
f