#!/usr/bin/env python3

import pandas as pd

# Read in the experiment design data from a tab-separated text file
experiment_design = pd.read_csv("exp_design.txt", sep='\t')

# Convert the column names to lower-case
experiment_design.columns = map(str.lower, experiment_design.columns)

# Create a DataFrame for the SDRF output
sdrf_output = experiment_design

# Fill in some default values
sdrf_output.columns = ['comment[file uri]', 'characteristics[experimental samples]']
sdrf_output['comment[data file]'] = experiment_design['comment[file uri]']
sdrf_output['comment[instrument]'] = "not available"
sdrf_output['comment[label]'] = 'AC=MS:1002038;NT=label free sample'
sdrf_output["comment[cleavage agent details]"] = 'NT=trypsin/P;AC=MS:1001313'

# Generate a list of sample numbers
sample_nums = list(range(1, len(experiment_design.index)+1))

# Insert a 'source name' column into the DataFrame
sdrf_output.insert(0, 'source name', ["Sample{}".format(n) for n in sample_nums])

# Copy the 'characteristics[experimental samples]' column to a new 'factor value[condition]' column
sdrf_output['factor value[condition]'] = experiment_design['characteristics[experimental samples]']

# Save the SDRF output to a tab-separated text file
with open("sdrf_local.tsv", "w") as f:
    sdrf_output.to_csv(f, sep="\t", index=False)

