import pandas as pd
import gzip
import re
import sys
from collections import defaultdict

pop_info = defaultdict(list)  #if the key doesn't exist, it makes an empty list for key and then add key
pop_infoFile = open(r"C:\Users\SciFunk\Downloads\Working\pop_locations.txt", "r")
for line in pop_infoFile:
    if line.startswith("sample"):
         continue
    spline = line.strip().split()
    pop_info[spline[1]].append(spline[0])
    
#creates a list of all sample numbers that are in the African continent
afrPops = ["POS", "REF", "ALT"]

for x in ["LWK", "MSL", "GWD", "ESN", "YRI"]:
    afrPops.extend(pop_info[x]) 
    
#fdat = gzip.open(sys.argv[1], "rt")
fdat = gzip.open(r"C:\Users\SciFunk\Downloads\Working\chr1small.vcf.gz","rt")
#dframe = pandas.read_table(fdat, sep='\t', comment = "#")
header = False
data = []
AADict = {}
for line in fdat:
    if line.startswith("#CHR"):
        headerline = line.strip().split()
        header = True
    elif not line.startswith("##") and header == True:
        spline = line.strip().split()
        data.append(spline)
        match = re.search(r'AA=(.)', spline[7]) #\w means match a single character, match = gives T/F
        if match:
            if match.group(1) == '|':
                AADict[spline[2]] = "."
            else:
                AADict[spline[2]] = match.group(1).upper()    #1 will take whatever's in . position
        else:
            AADict[spline[2]] = "."
df = pd.DataFrame.from_records(data, columns=headerline)
df.drop(['#CHROM', 'QUAL', "FILTER", "INFO", "FORMAT"], axis = 1, inplace=True) #removes useless columns
df2 = df[(df['REF'].map(len)==1) & (df['ALT'].map(len)==1)].copy() #makes sure we're only looking at bi-allelic
df2.set_index('ID', inplace=True)
dfAfr = df2[afrPops]

sumOnes = ((dfAfr.apply(lambda row: row.str.contains(r"1\|1").sum()*2, axis = 1) +
dfAfr.apply(lambda row: row.str.contains(r"1\|0").sum(), axis = 1)+
dfAfr.apply(lambda row: row.str.contains(r"0\|1").sum(), axis = 1))/(len(afrPops)*2-6))

sumZeroes = ((dfAfr.apply(lambda row: row.str.contains(r"0\|0").sum()*2, axis = 1) +
dfAfr.apply(lambda row: row.str.contains(r"1\|0").sum(), axis = 1)+
dfAfr.apply(lambda row: row.str.contains(r"0\|1").sum(), axis = 1))/(len(afrPops)*2-6))

df2.loc[:,"sumOnes"] = sumOnes
df2.loc[:,"sumZeroes"] = sumZeroes
df3 = df2[(df2.sumOnes < 0.01) | (df2.sumZeroes < 0.01)]

df3['CountThis'] = df3['sumOnes'].apply(lambda row: 1 if row < 0.01 else 0)

df3.drop(["sumOnes", "sumZeroes"], axis = 1, inplace = True)
#now df3 is only SNPS that occur <1% in Afr

finalDict = defaultdict(lambda: defaultdict(str))
for pop in pop_info:
    if pop in ["LWK", "MSL", "GWD", "ESN", "YRI"]:
        continue
    nonAfrPops = ["POS", "REF", "ALT", "CountThis"]
    nonAfrPops.extend(pop_info[pop])
    dftemp = df3[nonAfrPops]
    dftemp1 = dftemp[dftemp.CountThis == 1]
    dftemp0 = dftemp[dftemp.CountThis == 0]
    
    sumOnes = ((dftemp1.apply(lambda row: row.str.contains(r"1\|1").sum()*2, axis = 1) +
    dftemp1.apply(lambda row: row.str.contains(r"1\|0").sum(), axis = 1)+
    dftemp1.apply(lambda row: row.str.contains(r"0\|1").sum(), axis = 1))/(len(nonAfrPops)*2-8)) > 0.05
    
    sumZeroes = ((dftemp0.apply(lambda row: row.str.contains(r"0\|0").sum()*2, axis = 1) +
    dftemp0.apply(lambda row: row.str.contains(r"1\|0").sum(), axis = 1)+
    dftemp0.apply(lambda row: row.str.contains(r"0\|1").sum(), axis = 1))/(len(nonAfrPops)*2-8)) > 0.05
    
    for ID, row in dftemp1[sumOnes].iterrows():
        finalDict[pop][ID] = row["ALT"] 
        
    for ID, row in dftemp0[sumZeroes].iterrows():
        finalDict[pop][ID] = row["REF"] 
        
#for key in finalDict:  #format 'PEL': {'rs537182016': 'A', 'rs572818783': 'T'}
#    print("{}\t{}".format(key, len(finalDict[key])))

derivedDict = defaultdict(int)   #defaultdicts are nice because they will add your key and then add to it
ancestralDict = defaultdict(int)
unknownDict = defaultdict(int)

for pop in finalDict:
    for locus in finalDict[pop]:
        if AADict[locus] == '.':
            unknownDict[pop] += 1
        elif AADict[locus] == finalDict[pop][locus]:
            ancestralDict[pop] += 1
        else:
            derivedDict[pop] += 1
            
match = re.search(r'ALL.(chr\d+)', sys.argv[1])
with open("{}.txt".format(match.group(1), "w") as results:
            
    for pop in derivedDict:
        print("{}\t{}".format(pop, derivedDict[pop]), file = results)

    for pop in ancestralDict:
        print("{}\t{}".format(pop, ancestralDict[pop]), file = results)

    for pop in unknownDict:
        print("{}\t{}".format(pop, unknownDict[pop]), file = results)









        
        
        
    

