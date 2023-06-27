import pandas as pd
import sys

### This is a script by Marta Alvarez Presas from the Institute of Evolutionary Biology (IBE, UPF-CSIC) generated on 
### the 27th of june 2023 with the aim of prepare the input file for the GO-Figure! software in order to get figures
### from a GO term analysis. If you use it, please, cite.

### You need to give the name of the input file as an argument and also the name of the output file

# Generate a dataframe from a file with GO results
df = pd.read_csv(sys.argv[1], sep=',')

# Delete the first column, that contains no relevant information
df = df.drop('Unnamed: 0', axis=1)

# Rename columns
new_names = ['% GOterm', 'enrichment_P-value']

df.columns = new_names

# Delete rows with 0 value in the evalue column
df = df[(df['enrichment_P-value'] != 0)]

# Print dataframe on screen
print(df)

# Save the dataframe in a new file with the tsv format tab sepparated
#df.to_csv("BilateriaNovelSC1_GO_clean.tsv", sep='\t', index=False)
df.to_csv(sys.argv[2], sep='\t', index=False)
