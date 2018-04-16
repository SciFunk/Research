import gzip
import re
import sys
from collections import defaultdict
from multiprocessing import Pool
import pandas as pd

def snpAncestral(chromosome, header_line, pop_info):
    finalDict = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    
    afrPops = ["POS", "REF", "ALT"]
    for x in ["LWK", "MSL", "GWD", "ESN", "YRI"]:
        afrPops.extend(pop_info[x])
    #afrPops is a list of all samples NA19435 etc. in the African continent

    #Build initial dataframe
    df_iter = pd.read_table(chromosome, header = None, names = header_line, index_col = 2, engine = 'c',
                       compression = 'gzip', skiprows=251, comment = "#", chunksize = 40000, low_memory=False)

    for df in df_iter:
        df.drop(['#CHROM', 'QUAL', "FILTER", "FORMAT"], axis = 1, inplace=True) #removes useless columns

        #remove non-bi-allelic rows
        df['ALT'] = df['ALT'].apply(str)
        df['REF'] = df['REF'].apply(str)
        df.drop(df[df['ALT'].map(len) != 1].index, inplace = True)
        df.drop(df[df['REF'].map(len) != 1].index, inplace = True)

        #Select African pops
        dfAfr = df[afrPops]

        ##Calculate proportion of alleles for African pops different then reference for each SNP
        df = df.assign(sumOnes = ((dfAfr.apply(lambda row: row.str.contains(r"1\|1").sum()*2, axis = 1)+
                                   dfAfr.apply(lambda row: row.str.contains(r"1\|0").sum(), axis = 1)+
                                   dfAfr.apply(lambda row: row.str.contains(r"0\|1").sum(), axis = 1))/(len(afrPops)*2-6)))

        ##Calculate proportion of alleles for African pops same as reference for each SNP
        df = df.assign(sumZeroes = ((dfAfr.apply(lambda row: row.str.contains(r"0\|0").sum()*2, axis = 1)+
                      dfAfr.apply(lambda row: row.str.contains(r"1\|0").sum(), axis = 1)+
                      dfAfr.apply(lambda row: row.str.contains(r"0\|1").sum(), axis = 1))/(len(afrPops)*2-6)))

        df.drop(df.index[df.sumOnes > 0.01], inplace = True)
        ########################df.drop(df.index[df.sumZeroes > 0.01], inplace = True)
        #now df is only SNPS that occur <1% in Afr

        #creates a column that says whether to count 1's or 0's depending on which was <1% in Afr
        df = df.assign(CountThis = df['sumOnes'].apply(lambda row: 1 if row < 0.01 else 0))

        df.drop(["sumOnes", "sumZeroes"], axis = 1, inplace = True) #axis=1 means drop columns not rows


        #Extract Ancestral State
        df = df.assign(Ancestral = df['INFO'].str.extract('AA=(.)', expand=False).str.upper())
        df.drop("INFO", axis = 1, inplace = True)

        #get rows that are >5% in nonAfr pops
        for pop in pop_info:
            if pop in ["LWK", "MSL", "GWD", "ESN", "YRI"]:
                continue
            nonAfrPops = ["Ancestral", "POS", "REF", "ALT", "CountThis"]
            nonAfrPops.extend(pop_info[pop]) #nonAfrPops is a list of all samples HG19435 etc. not in African continent
            dftemp = df[nonAfrPops]
            dftemp1 = dftemp[dftemp.CountThis == 1]
            dftemp0 = dftemp[dftemp.CountThis == 0]

            sumOnes = ((dftemp1.apply(lambda row: row.str.contains(r"1\|1").sum()*2, axis = 1) +
            dftemp1.apply(lambda row: row.str.contains(r"1\|0").sum(), axis = 1)+
            dftemp1.apply(lambda row: row.str.contains(r"0\|1").sum(), axis = 1))/(len(nonAfrPops)*2-8)) > 0.05   
            sumZeroes = ((dftemp0.apply(lambda row: row.str.contains(r"0\|0").sum()*2, axis = 1) +
            dftemp0.apply(lambda row: row.str.contains(r"1\|0").sum(), axis = 1)+
            dftemp0.apply(lambda row: row.str.contains(r"0\|1").sum(), axis = 1))/(len(nonAfrPops)*2-8)) > 0.05  
             
            for sample in pop_info[pop]:
                for ID, row in dftemp1[sumOnes].iterrows():            
                    finalDict[sample][ID] = row["ALT"],row["POS"] 

                for ID, row in dftemp0[sumZeroes].iterrows():
                    finalDict[sample][ID] = row["REF"],row["POS"] 

        dfdictfinal = pd.DataFrame(finalDict)
        dfdictfinal = dfdictfinal.stack().reset_index()
        match = re.search(r'ALL.(chr\d+)', chromosome)
        dfdictfinal.to_csv(path_or_buf = r"{}.txt".format(match.group(1)), header=False, mode='w')

        chrome_dict = defaultdict(dict)
        chrome_dict[match.group(1)] = finalDict
        finalDict = chrome_dict[match.group(1)]
        finalDict = pd.DataFrame(finalDict)
        finalDict = finalDict.stack().reset_index()
        finalDict.to_csv(path_or_buf = r"chr{}.txt".format(match.group(1)), header=False, mode='w')
        return finalDict
    
def main():
    pop_info = defaultdict(list)  #if the key doesn't exist, it makes an empty list for key and then add key
    with open("pop_locations.txt", "r") as pop_infoFile:
        for line in pop_infoFile:
            if line.startswith("sample"):
                continue
            spline = line.strip().split()
            pop_info[spline[1]].append(spline[0])

    with open("1000.genome.header.txt", "r") as header:
        for line in header:
            header_line = line.strip().split()

    with Pool(processes = 22) as pool:
        jobs = []
        for chromo in sys.argv[1:]:
            jobs.append((chromo, header_line, pop_info))

	#apply snpAncestral fn to each file in jobs
        results = pool.starmap(snpAncestral, jobs)
        
if __name__ == "__main__":
    main()

