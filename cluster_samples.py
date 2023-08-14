### This script clusters samples with the cfRRBS cluster file

# import libraries
import glob
import os
import sys
import pandas as pd
import numpy as np
from functools import reduce
import subprocess
import pybedtools
from pybedtools import BedTool
pybedtools.cleanup()

# Parallelisation options
import multiprocessing
cpuCount = (multiprocessing.cpu_count())

# Read arguments from command line
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='path to samples files')
parser.add_argument('--clusterfile', '-c', help='path to cluster file')
parser.add_argument('--manifest', '-m', help='path to manifest file')
parser.add_argument('--output', '-o', help='output directory')
parser.add_argument('--NGS', action='store_true', help='use this argument if samples are Bismark coverage files')
parser.add_argument('--array', action='store_true', help='use this argument if samples are Illumina array files')
args = parser.parse_args()
NGS = args.NGS
array = args.array

# Input files
input = args.input
if array == True:
    samples = glob.glob(input + '*.txt')
if NGS == True:
    samples = glob.glob(input + '*.cov')
manifest = args.manifest
clusterfile = args.clusterfile

# Output files
output = args.output
sample_output = output + "samples_clustered.csv"
metadata_output = output + "metadata_clusters.csv"

# check if sample type is specified and import Illumina manifest file if necessary
if array == True:
    print("array enabled")
    # Load Infinium manifest file
    array450k = pd.read_csv(manifest, dtype={"CHR": str}, header = 7, usecols = (0,10,11,12), index_col="IlmnID", low_memory=False)
    array450k = array450k.dropna()
    array450k[['MAPINFO', 'Genome_Build']] = array450k[['MAPINFO', 'Genome_Build']].astype(int)
    # Extract locations with genome build GRCh37
    array450k = array450k[array450k['Genome_Build'] == 37]
    array450k = array450k.drop(['Genome_Build'], axis = 1)
    array450k[['CHR', 'MAPINFO']] = array450k[['CHR', 'MAPINFO']].astype(str)
    array450k.index.name = None
elif NGS == True:
    print("NGS enabled")
else:
    sys.exit("specify if the samples are NGS or array-based")

# prepare clusterfile to process with BEDtools
clusters_df = pd.read_csv(clusterfile, header=0, sep='\t', usecols=['chr','start','stop','clusterID'], low_memory=False)
clusters_BEDfile = BedTool.from_dataframe(clusters_df).sort()

# cluster samples and append to list
def cluster_samples(file):
    file_name = os.path.splitext(os.path.basename(file))[0]
    if array == True:
        name = file_name.split('.')[0]
        ID = file_name.split('.')[4] #temp change: ID = file_name.split('.')[0]
        label = ID.replace('TCGA',name) #temp change: label = ID.replace('GSM','melanocytes-')
        # convert TCGA sample to genome build GRCh37, and drop NAs
        sample_df = pd.read_csv(file, sep="\t", usecols=['Composite Element REF','Beta_value'], index_col=0, low_memory=False)#temp change: usecols=['ID_REF','VALUE']
        #sample_df.rename(columns={'VALUE':'Beta_value'},inplace=True)# temp change
        sample_df = pd.merge(array450k, sample_df, how="inner", left_index=True, right_index=True)
        sample_df.dropna(inplace=True)
        # add a stop and reorder the columns to use with BEDtools
        sample_df["MAPINFO_Stop"] = sample_df["MAPINFO"]
        sample_df = sample_df[["CHR", "MAPINFO", "MAPINFO_Stop", "Beta_value"]]
        sample_df.sort_values(by = ["CHR", "MAPINFO"], inplace=True)
        #don't forget if using start info of TCGA file directly: sample_df = sample_df[sample_df['Start'] != -1]
    elif NGS == True:
        label = file_name.split('_')[0]
        # cluster sample
        sample_df = pd.read_csv(file, sep="\t", header=None, names=['chr','start','stop','beta_value','meth','unmeth'], low_memory=False)
        # calculate beta value
        sample_df['beta_value'] = sample_df['beta_value']/100
    print(f"The number of CpGs in the {file_name} file is: {len(sample_df)}")
    # convert to bedfile, intersect with clusterfile, group by clusterID and count amount of CpGs per cluster
    sample_BEDfile = BedTool.from_dataframe(sample_df).sort()
    if array == True:
        sample_clusters_intersect = clusters_BEDfile.intersect(sample_BEDfile, sorted=True, wa=True, wb=True).groupby(g=[4], c=[8], o=['mean','count'])#alternative: ...loj=True).filter(lambda x: x[4] != ".")
        clusteredSample_df = sample_clusters_intersect.to_dataframe(disable_auto_names=True, header=None, names=['clusterID',label,label +' count'])
    elif NGS == True:
        sample_clusters_intersect = clusters_BEDfile.intersect(sample_BEDfile, sorted=True, wa=True, wb=True).groupby(g=[4], c=[8,8,9,10], o=['mean','count','sum','sum'])
        clusteredSample_df = sample_clusters_intersect.to_dataframe(disable_auto_names=True, header=None, names=['clusterID',label, label +' count','meth','unmeth'])
        # Get total depth (=methylated + unmethylated count)
        clusteredSample_df[label +' depth'] = clusteredSample_df['meth'].astype(float) + clusteredSample_df['unmeth'].astype(float)
        # only keep clusterID, beta value and total depth
        clusteredSample_df = clusteredSample_df.drop(columns=['meth','unmeth'])
    print(f"The number of clusters in {file_name} is: {len(clusteredSample_df)}")
    clusteredSample_df[label] = clusteredSample_df[label].astype(float)
    file_list.append(clusteredSample_df)

# create list with manager to be shared between processes
manager = multiprocessing.Manager()
file_list = manager.list()
# cluster samples in parallel
pool = multiprocessing.Pool(cpuCount)
pool.map(cluster_samples, samples)
# make dataframe from appended list and rename CpG ID
samples_df = reduce(lambda x,y: pd.merge(x,y,on='clusterID',how='outer'), file_list)
print("The total amount of clusters is: %s" % len(samples_df))
samples_df.sort_values(by='clusterID', inplace=True)
samples_df = samples_df.set_index('clusterID')

# create one file with beta values from all samples and another with amount of CpGs (for arrays) or amount of CpGs + coverage (for NGS) per cluster
metadata_list = list()
sample_list = list()
for col in list(samples_df):
    if 'depth' in col or 'count' in col:
        metadata_list.append(col)
    else:
        sample_list.append(col)
metadata_df = samples_df[metadata_list]
metadata_df.to_csv(metadata_output)
samples_df = samples_df[sample_list]
samples_df.to_csv(sample_output)
