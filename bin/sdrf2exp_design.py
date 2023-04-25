#!/usr/bin/env python3
# This code is used to convert SDRF file into experiment design file.
# It takes the SDRF file and creates a new experiment design file called "exp_design.txt"
# with the "raw_file" and "exp_condition" columns. The "raw_file" column contains the 
# raw data file names, and the "exp_condition" column contains the experimental conditions
# for each raw data file. We further get a "biorep" column with the number of the biological 
# replicate for each raw data file.

import pandas as pd
import os

sdrf = pd.read_csv("sdrf_local.tsv", sep='\t')
sdrf = sdrf.astype(str)
sdrf.columns = map(str.lower, sdrf.columns)  # convert column names to lower-case
changing_columns = list()
# Identify the columns (characteristics) that change between samples.
changing_columns = []
for col in sdrf.columns:
    if "factor" in col or "characteristics" in col:
        if sdrf[col].value_counts().size > 1:
            changing_columns.append(col)


# In case of having factor columns, take only them
changing_columns = [x for x in changing_columns if "factor" in x  ]

sdrf_out = pd.DataFrame()
sdrf_out["raw_file"] = [os.path.basename(item) for item in sdrf["comment[file uri]"]]
if (len(changing_columns) > 0):
	sdrf_out["exp_conditions"] = sdrf[changing_columns].agg('_'.join, axis=1)
else: 
	sdrf_out["exp_conditions"] = "A"
        
# Error when there are not column characteristics[biological replicate] and comment[technical replicate] and comment[fraction identifier] in sdrf file
if "characteristics[biological replicate]" not in sdrf.columns:
    print("Error: There is not column characteristics[biological replicate] in sdrf file")
    exit()
if "comment[technical replicate]" not in sdrf.columns:
    print("Error: There is not column comment[technical replicate] in sdrf file")
    exit()
if "comment[fraction identifier]" not in sdrf.columns:
    print("Error: There is not column comment[fraction identifier] in sdrf file")
    exit()

# Add column biorep from the column "characteristics[biological replicate]"
sdrf_out["biorep"] = sdrf["characteristics[biological replicate]"]
# Add column fraction from the column "comment[fraction identifier]"
sdrf_out["fraction"] = sdrf["comment[technical replicate]"]
# Add column techrep from the column "comment[technical replicate]"
sdrf_out["techrep"] = sdrf["comment[technical replicate]"]

sdrf_out.to_csv("exp_design.txt", sep="\t", index=False)
    
