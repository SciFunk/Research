import sys
import gzip
from collections import Counter


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

        total_Afr = 0
        total_Afr += pop_counts["LWK"]
        total_Afr += pop_counts["ESN"]
        total_Afr += pop_counts["GWD"]
        total_Afr += pop_counts["MSL"]
        total_Afr += pop_counts["YRI"]
        if total_Afr > 997:   #this leaves <1% 0's
            for pop, count in pop_counts.items():
                percentage = (count / pop_total[pop]) * 100
                if 95 < percentage < 100:
                    final_pop_counts[pop] += 1

with open("C:\\Users\\SciFunk\\Downloads\\Working\\1percentAfrican0s.txt", "a") as f:
    for k, v in sorted(final_pop_counts.items()):
        f.write('{}\t{}\n'.format(k, v))
