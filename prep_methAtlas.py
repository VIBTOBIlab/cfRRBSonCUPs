### this script prepares samples for deconvolution with methAtlas by renaming them, adding the right index and applying minimal coverage cut-off ###

## Import libraries
import glob
import os
import re
import numpy as np
import pandas as pd

## I/O files
sample_files = glob.glob("./reference/healthyWGBS_dennisLo/new-cluster_test/*samples-clustered.csv")
WBC = "./reference/cfRRBSreferences/new-clusters_test/2022-03-30_plasma-WBCref-samples_clustered_DMR200.csv"
refs = "./reference/cfRRBSreferences/new-clusters_test/2022-04-12_FFPE-ref-samples_clustered_allDMRs-8T.csv"
test = "./clustered_samples/2022-04-12_FFPE-CUPs-samples_clustered_allDMRs-8T.csv"

metadata_files = glob.glob("./reference/healthyWGBS_dennisLo/new-cluster_test/*metadata-clusters.csv")
WBC_metadata = "./reference/cfRRBSreferences/new-clusters_test/2022-03-30_plasma-WBCref-metadata_clusters_DMR200.csv"
ref_metadata = "./reference/cfRRBSreferences/new-clusters_test/2022-04-12_FFPE-ref-metadata_clusters_allDMRs-8T.csv"
test_metadata = "./clustered_samples/2022-04-12_FFPE-CUPs-metadata_clusters_allDMRs-8T.csv"

intersect_clusterfile = "./cfRRBS_450k_intersectClusters_without_sexChr.tsv"
names = "./sample_lists/sampleList_FFPE-CUPs.csv"
names = "./sample_lists/refList_age_v2.csv"
output = "./methAtlas/new-clusters_test/train_methAtlas_own/2022-04-01_cfRRBS_tumor+WBC_DMR200_lessthenhalfNAs_refs-methAtlas.csv"

## determine minimal coverage
coverage = 30

## apply coverage cut-off
samples_df = pd.read_csv(refs, index_col='clusterID')
samples_df.sort_values(by='clusterID', inplace=True)
# calculate the amount of clusters in each sample
cluster_amount = dict.fromkeys(list(samples_df))
for sample in list(samples_df):
    cluster_amount[sample] = len(samples_df)-samples_df[sample].isna().sum()
    print(f"The amount of clusters in {sample}: {cluster_amount[sample]}")
# remove clusters with low coverage and calculate amount of remaining clusters
metadata_df = pd.read_csv(ref_metadata, index_col='clusterID')
depth_df = metadata_df.filter(regex='depth')
depth_df.rename(columns = lambda x: x.split(' ')[0], inplace=True)
for sample in list(samples_df):
    samples_df[sample] = samples_df[sample].mask(depth_df[sample] < coverage)
    remaining_clusters = len(samples_df)-samples_df[sample].isna().sum()
    print(f"The amount of clusters in {sample} with coverage lower than {coverage}: {cluster_amount[sample] - remaining_clusters}")
    print(f"The amount of remaining clusters in {sample}: {remaining_clusters}")
samples_df.dropna(how='all', inplace=True)

## merge with array intersect-clusterfile
intersect_df = pd.read_csv('./2022-04-06_DMR200_v3clusters-merged.tsv', sep='\t', usecols=['clusterID'])
samples_df = intersect_df.merge(samples_df, on='clusterID', how='left')
# don't use this line when using samples as references:
samples_df = samples_df.fillna(samples_df.median())

## rename samples
names_df = pd.read_csv(names, sep=',',index_col=[0])
samples_df = samples_df.transpose()
samples_df = pd.merge(samples_df,names_df,how='left',left_index=True,right_index=True)
samples_df = samples_df.set_index('name')
samples_df = samples_df.transpose()
#samples_df['clusterID'] = samples_df['clusterID'].astype(int)
samples_df['clusterID'] = samples_df.index
samples_df.index = samples_df['clusterID']
samples_df.drop(columns=['clusterID'], inplace=True)
samples_df = samples_df[samples_df.columns.dropna()]

## using cfRRBS samples as references: group samples per tumor type and take median
### OR first group per tumor type
samples_df.columns = samples_df.columns.str.split(' ').str[0]
### remove DMRs with to many Na values (here: less than half)
NAs = samples_df.isnull().groupby(by=samples_df.columns, axis=1).sum()/samples_df.groupby(by=samples_df.columns, axis=1).size()
NAs['max'] = NAs.max(axis=1).astype(int)
NAs = NAs[NAs['max'] < 0.5]
NAs = NAs[['max']]
samples_df = samples_df.merge(NAs, on='clusterID')
samples_df.drop(columns=['max'],inplace=True)
### take median per group
samples_df = samples_df.groupby(by=samples_df.columns, axis=1).median()
###OR take median of all columns
samples_df = samples_df.median(axis=1).to_frame().rename(columns={0:'healthy'})


## add index for methAthlas
samples_df['CpGs'] = 'cg' + samples_df['clusterID'].astype(str)
samples_df = samples_df.set_index("CpGs")
samples_df.columns.name = None
samples_df.drop(columns=['clusterID'], inplace=True)

samples_df.to_csv(output)
