* To create a file with all samples, rsID, chr, pos, allele (<1% in Afr pops, >5% in other pop): IDs_perSamplewithPos.py (submission file: IDs.qsub, runtime ~few hours)
  * creates a separate file for each chr, combine with combineFiles.py (adds a column with filename title)
  * combined file messy with extra characters, use cleanSamples.py to clean it up 
      * file has colums: chr, pos, rsID, sample, allele
      * *stored as Working/IDs/allchrSampleswithPos_CLEAN.csv*
