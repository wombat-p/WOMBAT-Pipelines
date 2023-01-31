#!/usr/bin/env python3

import pandas as pd

expd = pd.read_csv("exp_design.txt", sep='\t')
expd.columns = map(str.lower, expd.columns)  # convert column names to lower-case
sdrf_out = expd
# filling in some default values
sdrf_out.columns = ['comment[file uri]', 'characteristics[experimental samples]']
sdrf_out['comment[data file]'] = expd['comment[file uri]']
sdrf_out['comment[instrument]'] = "not available"
sdrf_out['comment[label]'] = 'AC=MS:1002038;NT=label free sample'
sdrf_out["comment[cleavage agent details]"] = 'NT=trypsin/P;AC=MS:1001313'
sample_nums = list(range(1, len(expd.index)+1))
sdrf_out.insert(0, 'source name', [m + str(n) for m,n in zip(["Sample"]*len(expd.index), sample_nums)])
sdrf_out['factor value[condition]'] = expd['characteristics[experimental samples]']

sdrf_out.to_csv("sdrf_local.tsv", sep="\t", index=False)
    

