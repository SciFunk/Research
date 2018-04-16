import glob
import os

out_filename = "allchrSampleswithPos.csv"

read_files = glob.glob("chr*.txt")
with open(out_filename, "w") as outfile:
    for filename in read_files:
        with open(filename) as infile:
            for line in infile:
                outfile.write('{},{}\n'.format(line.strip(), filename))
