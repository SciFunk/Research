import pandas as pd
import os

#combined file with all chr: 'allchrSampleswithPos.csv'
dfpos = pd.read_csv('allchrSampleswithPos.csv', dtype="str", sep=",", header=None, names = ['rs', 'sample', 'allelepos', 'chr'])

#split allelepos line, which has both the allele and chr pos in one cell
df2 = dfpos['allelepos'].str.split(',', expand=True).rename(columns={0:'allele', 1:'pos'})
df3 = pd.concat([dfpos, df2], axis=1)
df3.drop(['allelepos'], axis=1, inplace=True)

#strip extra characters from columns
df3['chr'] = df3['chr'].map(lambda x: x.strip('chr'))
df3['chr'] = df3['chr'].map(lambda x: x.strip('.txt'))
df3['allele'] = df3['allele'].map(lambda x: x.strip('('))
df3['allele'] = df3['allele'].map(lambda x: x.strip("'"))
df3['pos'] = df3['pos'].map(lambda x: x.strip(')'))

#reorder columns
df3 = df3[['chr', 'pos', 'rs', 'sample', 'allele']]

df3.to_csv(path_or_buf = "allchrSampleswithPos_CLEAN.csv")
