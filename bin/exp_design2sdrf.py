#!/usr/bin/env python3

# Creating sdrf file template automatically from experimental design file
# We assume priority of biological replicates and fractions, and then fill 
# up technical replicates if not present in file. That means the non-existing
# will be filled without considering any technical ones.

import pandas as pd

expd = pd.read_csv("exp_design.txt", sep='\t')
expd.columns = map(str.lower, expd.columns)  # convert column names to lower-case
sdrf_out =  pd.DataFrame()
sdrf_out['comment[file uri]'] = expd.raw_file
sdrf_out['characteristics[experimental sample]'] = expd.exp_condition

# filling in some default values
sdrf_out['comment[data file]'] = expd['raw_file']
sdrf_out['comment[instrument]'] = "not available"
sdrf_out['comment[label]'] = 'AC=MS:1002038;NT=label free sample'
sdrf_out["comment[cleavage agent details]"] = 'NT=trypsin/P;AC=MS:1001313'

sdrf_out['factor value[condition]'] = sdrf_out['characteristics[experimental sample]']

# if fraction column not present, add it as comment[fraction identifier] filled with 1
if "fraction" not in expd.columns:
    sdrf_out["comment[fraction identifier]"] = "1"
else:
    sdrf_out["comment[fraction identifier]"] = expd["fraction"]

# if biorep column is not present, add it as characteristics[biological replicate]
if "biorep" not in expd.columns:
    # Take all rows with comment[fraction identifier] equal to one
    # and group them by the column characteristics[experimental sample]
    # and add a column with increasing numbers for each value in the column characteristics[experimental sample]
    expd.loc[sdrf_out["comment[fraction identifier]"] == "1", "biorep"] = expd.loc[sdrf_out["comment[fraction identifier]"] == "1", "exp_condition"].groupby(expd["exp_condition"]).cumcount() + 1
    # Take all rows with comment[fraction identifier] not equal to one
    # set their biorep column to the value of the biorep column of the row with the same exp_condition and comment[fraction identifier] equal to one
    expd.loc[sdrf_out["comment[fraction identifier]"] != "1", "biorep"] = expd.loc[sdrf_out["comment[fraction identifier]"] != "1"].merge(expd.loc[sdrf_out["comment[fraction identifier]"] == "1", ["exp_condition", "biorep"]], on="exp_condition", how="left")["biorep_y"]
sdrf_out["characteristics[biological replicate]"] = expd["biorep"]

# if techrep column not present, add it as comment[technical replicate] filled with 1
if "techrep" not in expd.columns:
    sdrf_out["comment[technical replicate]"] = "1"
    # Take all rows with comment[fraction identifier] equal to one
    # and group them by the column characteristics[experimental sample] and characteristics[biological replicate]
    # and add a column with increasing numbers for each value in the column characteristics[experimental sample] and characteristics[biological replicate]
    expd.loc[sdrf_out["comment[fraction identifier]"] == "1", "techrep"] = expd.loc[sdrf_out["comment[fraction identifier]"] == "1",["exp_condition", "biorep"]].groupby(['exp_condition', 'biorep']).cumcount() + 1
    # Take all rows with comment[fraction identifier] not equal to one
    # set their techrep column to the value of the techrep column of the row with the same exp_condition, biorep and comment[fraction identifier] equal to one
    expd.loc[sdrf_out["comment[fraction identifier]"] != "1", "techrep"] = expd.loc[sdrf_out["comment[fraction identifier]"] != "1"].merge(expd.loc[sdrf_out["comment[fraction identifier]"] == "1", ["exp_condition", "biorep", "techrep"]], on=["exp_condition", "biorep"], how="left")["techrep_y"]    
sdrf_out["comment[technical replicate]"] = expd["techrep"]

# Add a column "source name" to sdrf_out with increasing numbers for each value in the column characteristics[experimental sample] and characteristics[biological replicate]
sdrf_out["source name"] = expd[["exp_condition", "biorep"]].apply(lambda x: '_'.join(x.astype(str)), axis=1)

#sample_nums = list(range(1, len(expd.index)+1))
#sdrf_out.insert(0, 'source name', [m + str(n) for m,n in zip(["Sample"]*len(expd.index), sample_nums)])

# Sort columns of sdrf_out to start with source name, then all columns starting with characteristics, then all columns starting with comment, then all columns starting with factor value
sdrf_out = sdrf_out.reindex(sorted(sdrf_out.columns, key=lambda x: (x.startswith('factor value'),x.startswith('comment'), x.startswith('characteristics'), x.startswith('sourcename'))), axis=1)

sdrf_out.to_csv("sdrf_local.tsv", sep="\t", index=False)
    

