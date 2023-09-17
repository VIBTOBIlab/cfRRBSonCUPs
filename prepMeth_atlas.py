### this script prepares samples for deconvolution with methAtlas by renaming them, adding the right index and applying minimal coverage/count cut-off ###

## Import libraries
import glob
import os
import re
import numpy as np
import pandas as pd
from functools import reduce
import pybedtools
from pybedtools import BedTool
x = BedTool()
pybedtools.cleanup()
from scipy import stats

# Read arguments from command line
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='.csv file of all samples')
parser.add_argument('--metadata', '-m', help='.csv file of all sample metadata')
parser.add_argument('--output', '-o', help='output directory')
parser.add_argument('--names', '-n', help='list with sample names')
parser.add_argument('--refs', action='store_true', help='use this argument if you want to use your samples as references in methAtlas')
args = parser.parse_args()
input = args.input
metadata = args.metadata
names = args.names
refs = args.refs

# Output files
output = args.output
if refs == False:
    output = output + "methAtlas_test.csv"
if refs == True:
    output = output + "methAtlas_train.csv"

## determine minimal coverage or count
filter = 30

## apply coverage or count cut-off
samples_df = pd.read_csv(input, index_col='clusterID')
samples_df.sort_values(by='clusterID', inplace=True)
# calculate the amount of clusters in each sample
cluster_amount = dict.fromkeys(list(samples_df))
for sample in list(samples_df):
    cluster_amount[sample] = len(samples_df)-samples_df[sample].isna().sum()
    print(f"The amount of clusters in {sample}: {cluster_amount[sample]}")
# remove clusters with low coverage/count and calculate amount of remaining clusters
metadata_df = pd.read_csv(metadata, index_col='clusterID')
filter_df = metadata_df.filter(regex='depth')
filter_df.rename(columns = lambda x: x.rstrip(" depth"), inplace=True)
for sample in list(samples_df):
    samples_df[sample] = samples_df[sample].mask(filter_df[sample] < filter)
    remaining_clusters = len(samples_df)-samples_df[sample].isna().sum()
    print(f"The amount of clusters in {sample} with coverage/count lower than {filter}: {cluster_amount[sample] - remaining_clusters}")
    print(f"The amount of remaining clusters in {sample}: {remaining_clusters}")
samples_df.dropna(how='all', inplace=True)

## rename en drop redundant samples
names_df = pd.read_csv(names, sep=',',index_col=[0],usecols=['smapleID','name'])
samples_df.rename(columns=names_df['name'], inplace=True)
they_who_must_not_be_named = [name for name in samples_df.columns if name not in names_df['name'].values.tolist()]
print(f"there are {len(they_who_must_not_be_named)} samples who will be dropped because their name is not in the list")
samples_df.drop(columns = they_who_must_not_be_named, inplace=True)

## for references: group samples per tumor type and take median
if refs == True:
    print("the samples are processed as references")
### make sure all columns with same tumor type have same name
    samples_df.columns = samples_df.columns.str.split(' ').str[0]
### remove DMRs with to many Na values for a certain entity (here: less than half)
    NAs = samples_df.isnull().groupby(by=samples_df.columns, axis=1).sum()/samples_df.groupby(by=samples_df.columns, axis=1).size()
    NAs['max'] = NAs.max(axis=1)
    NAs = NAs[NAs['max'] < 0.5]
    NAs = NAs[['max']]
    samples_df = samples_df.merge(NAs, on='clusterID')
    samples_df.drop(columns=['max'],inplace=True)
### take median per group
    samples_df = samples_df.groupby(by=samples_df.columns, axis=1).median()

## add index for methAthlas
samples_df['CpGs'] = 'cg' + samples_df.index.astype(str)
samples_df = samples_df.set_index("CpGs")

## output
samples_df.to_csv(output)
