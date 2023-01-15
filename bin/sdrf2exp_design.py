#!/usr/bin/env python3

import pandas as pd
import os

# Read in the SDRF data from a tab-separated text file
sdrf_data = pd.read_csv("sdrf_local.tsv", sep='\t')

# Convert all columns to strings
sdrf_data = sdrf_data.astype(str)

# Convert the column names to lower-case
sdrf_data.columns = map(str.lower, sdrf_data.columns)

# Create a list of columns that have more than one unique value
changing_columns = []
for col in sdrf_data.columns:
    if "factor" in col or "characteristics" in col:
        if sdrf_data[col].value_counts().size > 1:
            changing_columns.append(col)

# Filter the list of columns to only include those with the "factor" prefix
changing_columns = [x for x in changing_columns if "factor" in x]

# Create an empty DataFrame for the output
output_data = pd.DataFrame()

# Add a 'raw_file' column with the base names of the file URIs
output_data["raw_file"] = [os.path.basename(item) for item in sdrf_data["comment[file uri]"]]

# Add an 'exp_conditions' column with the concatenated values from the changing columns
if (len(changing_columns) > 0):
    output_data["exp_conditions"] = sdrf_data[changing_columns].agg('_'.join, axis=1)
else: 
    output_data["exp_conditions"] = "A"

# Save the output data to a tab-separated text file
try:
    with open("exp_design.txt", "w") as f:
        output_data.to_csv(f, sep="\t", index=False)
except IOError as e:
    print("An error occurred while writing the output file:", e)
