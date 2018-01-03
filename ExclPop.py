import sys
import gzip
from collections import Counter

def ExclPopCount(ExclPop):
    ExclPop_counts = pop_counts["LWK"] + pop_counts["MSL"] + pop_counts["GWD"] + pop_counts["ESN"] + pop_counts["YRI"] 
    + pop_counts["PEL"] + pop_counts["MXL"] + pop_counts["CLM"] + pop_counts["PUR"] + pop_counts["JPT"]
    + pop_counts["CDX"] + pop_counts["KHV"] + pop_counts["CHB"] + pop_counts["CHS"] + pop_counts["FIN"]
    + pop_counts["CEU"] + pop_counts["GBR"] + pop_counts["IBS"] + pop_counts["TSI"] + pop_counts["GIH"]
    + pop_counts["ITU"] + pop_counts["STU"] + pop_counts["PJL"] + pop_counts["BEB"] - pop_counts[ExclPop]
    ExclPop_total = pop_total["LWK"] + pop_total["MSL"] + pop_total["GWD"] + pop_total["ESN"] + pop_total["YRI"] 
    + pop_total["PEL"] + pop_total["MXL"] + pop_total["CLM"] + pop_total["PUR"] + pop_total["JPT"]
    + pop_total["CDX"] + pop_total["KHV"] + pop_total["CHB"] + pop_total["CHS"] + pop_total["FIN"]
    + pop_total["CEU"] + pop_total["GBR"] + pop_total["IBS"] + pop_total["TSI"] + pop_total["GIH"]
    + pop_total["ITU"] + pop_total["STU"] + pop_total["PJL"] + pop_total["BEB"] - pop_total[ExclPop]
    totalPercentage = (ExclPop_counts/ExclPop_total) * 100
    if totalPercentage < 0.5:
        for pop, count in pop_counts.items():
            percentage = (count / pop_total[ExclPop]) * 100
            if percentage > 5:
                final_pop_counts[ExclPop] += 1

pop_map = {}
#with open('pop_locations.txt') as data:
with open('C:\\Users\\SciFunk\\Downloads\\Working\\pop_locations.txt') as data:
    for line in data:
        cells = line.split()
        ID, pop = cells[:2]
        pop_map[ID] = pop
        
final_pop_counts = Counter()
#with gzip.open(sys.argv[1]) as data:
with gzip.open('C:\\Users\\SciFunk\\Downloads\\Working\\chr1small.vcf.gz', "rt") as data:
    for line in data:
        if line.startswith('"##') or line.startswith("##"):
            continue
        if line.startswith("#"): # is header
            header = line.strip().split('\t')[9:]
            header_length = len(header)  # set once for use in calculation below
            continue
        bits = line.strip().split('\t')[9:]
        bits = [sum(map(int, b.split('|'))) for b in bits]
        pop_counts = Counter()
        pop_total = Counter()
        for ID, snps in zip(header, bits):
            pop = pop_map[ID]    #gives, GBR, LWK, etc; pop_map is dict: {H1: GBR, H2: LWK, ...}
            pop_counts[pop] += snps #pop_counts is Counter that adds 0,1,2 based on SNPs
            pop_total[pop] += 2 #counts total population numbers
        ExclPopCount("LWK")
        print(final_pop_counts)
        ExclPopCount("MSL")
        ExclPopCount("GWD")
        ExclPopCount("ESN")
        ExclPopCount("YRI")
        ExclPopCount("PEL")
        ExclPopCount("MXL")
        ExclPopCount("CLM")
        ExclPopCount("PUR")
        ExclPopCount("JPT")
        ExclPopCount("CDX")
        ExclPopCount("KHV")
        ExclPopCount("CHB")
        ExclPopCount("CHS")
        ExclPopCount("FIN")
        ExclPopCount("CEU")
        ExclPopCount("GBR")
        ExclPopCount("IBS")
        ExclPopCount("TSI")
        ExclPopCount("GIH")
        ExclPopCount("ITU")
        ExclPopCount("STU")
        ExclPopCount("PJL")
        ExclPopCount("BEB")
        
print(final_pop_counts)

# with open("C:\\Users\\SciFunk\\Downloads\\Working\\1percentAfrican.txt", "a") as f:
#     for k, v in sorted(final_pop_counts.items()):
#         f.write('{}\t{}\n'.format(k, v))
