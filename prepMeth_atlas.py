### this script prepares samples for deconvolution with methAtlas by renaming them, adding the right index and applying minimal coverage/count cut-off ###
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 20)

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

## I/O files
test = "./clustered_samples/2022-08-30_all-samples_clustered_DSS-DMRs-12G.csv"
test_metadata = "./clustered_samples/2022-08-30_all-metadata_clusters_DSS-DMRs-12G.csv"

names = "./sample_lists/final_refList.csv"
names = "./sample_lists/refList_age.csv"
names = "./sample_lists/refList_age_v2.csv"
names = "./sample_lists/refList_grouped_v4.csv"
names = "./sample_lists/sampleList_FFPE-CUPs.csv"
names = "./sample_lists/sampleList_plasma.csv"
names = "./sample_lists/sampleList_PE.csv"
names = "./sample_lists/sampleList_asc.csv"

## determine minimal coverage or count
# minimal coverage for NGS
filter = 30
# minimal count for array
filter = 3
# minimal stdev
stdev = 0.2

## apply coverage or count cut-off
samples_df = pd.read_csv(test, index_col='clusterID')
samples_df.sort_values(by='clusterID', inplace=True)
# calculate the amount of clusters in each sample
cluster_amount = dict.fromkeys(list(samples_df))
for sample in list(samples_df):
    cluster_amount[sample] = len(samples_df)-samples_df[sample].isna().sum()
    print(f"The amount of clusters in {sample}: {cluster_amount[sample]}")
# remove clusters with low coverage/count and calculate amount of remaining clusters
metadata_df = pd.read_csv(test_metadata, index_col='clusterID')
filter_df = metadata_df.filter(regex='depth') # if NGS: regex='depth', if array: regex='count'
filter_df.rename(columns = lambda x: x.rstrip(" depth"), inplace=True) # if TCGA or own samples: x.split(' ')[0]
for sample in list(samples_df):
    samples_df[sample] = samples_df[sample].mask(filter_df[sample] < filter)
    remaining_clusters = len(samples_df)-samples_df[sample].isna().sum()
    print(f"The amount of clusters in {sample} with coverage/count lower than {filter}: {cluster_amount[sample] - remaining_clusters}")
    print(f"The amount of remaining clusters in {sample}: {remaining_clusters}")
samples_df.dropna(how='all', inplace=True)
# remove clusters with high stdev
#filter_df = metadata_df.filter(regex='stdev')
#filter_df.rename(columns = lambda x: x.rstrip(" stdev"), inplace=True)
#for sample in list(samples_df):
#    samples_df[sample] = samples_df[sample].mask(filter_df[sample] > stdev)
#    remaining_clusters = len(samples_df)-samples_df[sample].isna().sum()
#    print(f"The amount of clusters in {sample} with stdev higher than {stdev}: {cluster_amount[sample] - remaining_clusters}")
#    print(f"The amount of remaining clusters in {sample}: {remaining_clusters}")
#samples_df.dropna(how='all', inplace=True) #how='all'
# check missing values
samples_df.isna().sum().sort_values()
samples_df

## merge with consistenly covered (min 30x) clusters
index = pd.read_csv('/Users/jilke/Scripts/RRBSonCUPs/vsc42927/cluster_files/cfRRBSclusters_v3_min30x.csv')
###### temp ########
index['CpGs'] = 'cg' + index['clusterID'].astype(str)
index = index.set_index("CpGs")
samples_df = index.merge(samples_df, left_index=True, right_index=True)
###################
samples_df = index.merge(samples_df, on='clusterID', how='left')
samples_df.index = samples_df['clusterID']
samples_df.drop(columns=['clusterID'], inplace=True)
samples_df.dropna(how='all', inplace=True)
samples_df = samples_df[samples_df.index < 203098]

## rename samples
names_df = pd.read_csv(names, sep=',',index_col=[0],usecols=['smapleID','name'])
names_df
samples_df = samples_df.transpose()
samples_df = pd.merge(samples_df,names_df,how='left',left_index=True,right_index=True)
samples_df = samples_df.set_index('name')
samples_df = samples_df.transpose()
samples_df['clusterID'] = samples_df.index
samples_df.index = samples_df['clusterID']
samples_df.drop(columns=['clusterID'], inplace=True)
samples_df = samples_df[samples_df.columns.dropna()]

## rename samples
names_df = pd.read_csv(names, sep=';',index_col=[0],usecols=['smapleID','diagnosis'])
names_df
samples_df.rename(columns=names_df['diagnosis'], inplace=True)
samples_df.drop(columns = list(samples_df.filter(regex='[0-9]')), inplace=True)

## filter refs for hierarchical deconvolution
samples_df = samples_df.filter(regex=r'(healthy|CHOL|GEJC)')
samples_df = samples_df.filter(regex=r'(healthy|LUSC|ESCC|CSCC)')
samples_df = samples_df.filter(regex=r'(healthy|BRCA)')

HIGI = samples_df[['CUP 2','CUP PAAD 1 -> CUP CHOL 1','CUP PAAD 2','CUP PAAD 3','CUP CHOL 2','CUP CARD 1','CUP CARD 2','CUP ESAD 1','CUP ESAD 2']]
HIGI
HIGIpl = samples_df[['cfCHOL 1','cfCHOL 2','cfPAAD 1','cfPAAD 2','cfPAAD 3','cfCUP 1','cfCUP 4']]
HIGIpl
samples_df = HIGI.merge(HIGIpl, on='clusterID')
SCCpl = samples_df[['cfCUP 8', 'cfSCC 1','cfSCC 2']]
SCCpl
SCCpe = samples_df[['cfpESSC 1', 'cfpLUSC 2', 'cfpLUSC 3']]
SCCpe
samples_df = SCCpe.merge(SCCpl, on='clusterID')
BRCA = samples_df[["D1819227","D1821215"]]
BRCA = samples_df[["DNA056033"]]
samples_df = BRCA
samples_df

## for references: group samples per tumor type and take median
### first remove rows with too many Na values
#samples_df = samples_df.loc[samples_df.isnull().mean(axis=1) < 0.5]
### OR first group per tumor type
samples_df.columns = samples_df.columns.str.split(' ').str[0]
samples_df.drop(columns=['ESADy'],inplace=True)
### remove DMRs with to many Na values (here: less than half)
NAs = samples_df.isnull().groupby(by=samples_df.columns, axis=1).sum()/samples_df.groupby(by=samples_df.columns, axis=1).size()
NAs['max'] = NAs.max(axis=1)
NAs = NAs[NAs['max'] < 0.5]
NAs = NAs[['max']]
len(NAs)
samples_df = samples_df.merge(NAs, on='clusterID')
samples_df.drop(columns=['max'],inplace=True)
### take median per group
samples_df = samples_df.groupby(by=samples_df.columns, axis=1).median()
### DMR selection: highest inter variance
inter_var = samples_df.var(axis=1).nlargest(2000)
samples_df = samples_df[samples_df.index.isin(inter_var.index)]
### take top 100 hypo-/hypermethylated
df_lst = []
for ref in list(samples_df):
    hypo = samples_df.nsmallest(400, ref)
    df_lst.append(hypo)
    hyper = samples_df.nlargest(400, ref)
    df_lst.append(hyper)
df = reduce(lambda x,y: pd.merge(x,y, right_index=True, left_index=True, how='outer'), df_lst)
samples_df = samples_df[samples_df.index.isin(df.index)]
### take max/min methylated
df_lst = []
for ref in list(samples_df):
    df = samples_df[[ref]]
    df = df[(df[ref] > 0.9) | (df[ref] < 0.03)]
    df_lst.append(df)
df = reduce(lambda x,y: pd.merge(x,y, right_index=True, left_index=True, how='outer'), df_lst)
samples_df = samples_df[samples_df.index.isin(df.index)]
###OR take median of all columns
samples_df = samples_df.median(axis=1).to_frame().rename(columns={0:'CONTR, WBC'})

## add index for methAthlas
samples_df['CpGs'] = 'cg' + samples_df.index.astype(str)
samples_df = samples_df.set_index("CpGs")

## output
samples_df.to_csv('./methAtlas/new-clusters_test/train_methAtlas_own/2022-09-15_cfRRBSwithmeso_12G_DSS-DMRs_lessthenhalfNAs_refs-methAtlas.csv')
samples_df.to_csv('./methAtlas/new-clusters_test/2022-09-12_all-samples_30x_DSS-DMRs-SCC_methAtlas.csv')

ESAD.to_csv('./methAtlas/new-clusters_test/2022-07-14_cfESAD_30x_DSS-DMRs-12G_samples-methAtlas.csv')
non_CUPs.to_csv('./methAtlas/new-clusters_test/2022-07-14_cfnon-CUPs_30x_DSS-DMRs-12G_samples-methAtlas.csv')
CUPs.to_csv('./methAtlas/new-clusters_test/2022-07-14_cfCUPs_30x_DSS-DMRs-12G_samples-methAtlas.csv')

samples_df.to_csv('./methAtlas/v3clusters/train_methAtlas_public/2022-06-16_TCGA+WBC_count3-30xindex-topHYPO-HYPER_refs-methAtlas.csv')
samples_df.to_csv('./methAtlas/v3clusters/2022-06-16_FFPE-CUPs_30x_cfRRBSclustersV3-methAtlas.csv')
samples_df.to_csv('./methAtlas/v3clusters/2022-05-10_plasma_10x_cfRRBSclustersV3-methAtlas.csv')

samples_df.to_csv('./methAtlas/old/train_methAtlas/2022-06-21_TCGA+WBC-sameOLD_count3-RRBS-450k_refs-methAtlas.csv')
samples_df.to_csv('./methAtlas/old/2022-06-16_FFPE-CUPs_30x-fillNA_RRBS-450k-methAtlas.csv')
samples_df.to_csv('./clustered_samples/2022-08-03_FFPE-refs_DSS-DMRs-17T_cov30x.csv')

samples_df

samples_df.to_csv('/Users/jilke/Scripts/cfRRBSvsEPIC-CNV/2022-06_methAtlas/v3clusters/train_methAtlas/2022-06-16_HB-count3-30xindex-minSTDEV-minHYPERmaxHYPO_refs-methAtlas.csv')
samples_df.to_csv('/Users/jilke/Scripts/cfRRBSvsEPIC-CNV/2022-06_methAtlas/v3clusters/2022-06-15_EPIC-FFPE_30xcov-noNA_samples-methAtlas.csv')

samples_df = pd.read_csv('./methAtlas/new-clusters_test/2022-06-23_PE_30x_DSS-DMRs-17T_methAtlas.csv',index_col='CpGs')
samples_df = samples_df.drop(columns=['BLCA'])
samples_df = samples_df[samples_df.index.isin(test_df.index)]
samples_df = samples_df[['ACC','BRCA','CEAD','CHOL','COAD','ESAD','KIRC','KIRP','LIHC','LUAD','PAAD','PRAD','SKCM','THCA','UCEC','healthy']]
ruby_subset = samples_df[samples_df.index.isin(size_df.index)]

# sort & select samples

lst = list(samples_df)
lst.sort()
samples_df = samples_df[lst]

ESAD = samples_df.filter(regex='ESAD')
ESAD.dropna(how='all', inplace=True)
ESAD

non_CUPs = samples_df[['cfPAAD 1','cfPAAD 2','cfLUAD 1','cfLUAD 2','cfSKCM 1','cfSKCM 2','cfPRAD 1','cfPRAD 2']]
non_CUPs.dropna(how='all', inplace=True)
non_CUPs

CUPs = samples_df.drop(columns=list(ESAD)+list(non_CUPs))
CUPs.dropna(how='all', inplace=True)
CUPs

names_df = pd.read_csv(names, sep=';',index_col=[0])
neg = names_df[names_df['diagnosis'].str.contains('?', regex=False) | names_df['diagnosis'].str.contains('-')]
pos = samples_df.drop(columns=list(neg['name']))
samples_df = samples_df[list(pos)+list(neg['name'])]

samples_df = samples_df[['DNA075044','DNA075069','DNA075049','DNA075047','DNA075050','DNA075054','DNA075040','DNA075045','DNA075037','DNA075070','DNA075055','DNA075063','DNA075072','DNA075071','DNA075058','DNA075064','DNA075068','DNA075056','DNA075057','DNA075035','DNA075041','DNA075042','DNA075066','DNA075034','DNA075043','DNA075067','DNA075065','DNA075052','DNA075062','DNA075059','DNA075038','DNA075060','DNA075033','DNA075046','DNA075048','DNA075051','DNA075053','DNA075036','DNA075061','DNA075039']]

## analyse methAtlas results
pd.set_option('display.max_columns', 300)
pd.set_option('display.max_rows', 300)
results = pd.read_csv('/Users/jilke/Downloads/test_methatlas_20190110_brain_MANUSCRIPT_filter30_all_GRCh37_deconv_output.csv', index_col=0)
results.columns = results.columns.str.split('_').str[0]
results = results[['DNA075036','DNA075066','DNA075071','DNA075044','DNA075043','DNA075040','DNA075061','DNA075050','DNA075034','DNA075059','DNA075042','DNA075052','DNA075049','DNA075045','DNA075048','DNA075039','DNA075055','DNA075064','DNA075067','DNA075046','DNA075069','DNA075070','DNA075063','DNA075054','DNA075037','DNA075072','DNA075058','DNA075068','DNA075056','DNA075057','DNA075035','DNA075041','DNA075065','DNA075062','DNA075038','DNA075060','DNA075033','DNA075051','DNA075053','DNA075047']]

results = pd.read_csv('/Users/jilke/Scripts/cfRRBSvsEPIC-CNV/2022-06_methAtlas/v3clusters/2022-06-14_EPIC-FFPE_30xcov_samples-methAtlas_deconv_output.csv', index_col=0)
results = pd.read_csv('./methAtlas/new-clusters_test/2022-07-14_asc_30x_DSS-DMRs-12G_methAtlas_deconv_output.csv', index_col=0)
# select results
results = results[list(neg['name'])]
results = results.drop(columns=list(neg['name']))
neg

# ranking
ranking = pd.DataFrame(columns=['1st diagnosis','1st score','2nd diagnosis','2nd score','3rd diagnosis','3rd score'])
for sample in list(results):
    df = results[[sample]].sort_values(by=sample, ascending=False)
    l = 0
    for d, s in zip(list(ranking.filter(regex='diagnosis')),list(ranking.filter(regex='score'))):
        for i in range(l,len(df)):
            if not 'healthy' in df.iloc[i].name:
                ranking.loc[sample,d] = df.iloc[i].name
                ranking.loc[sample,s] = df.iloc[i][sample]
            else:
                ranking.loc[sample,d] = 0
            if ranking.loc[sample,d] == df.iloc[i].name:
                l = i+1
                break

ranking
results
results.drop('healthy', inplace=True)
results.loc['total tumor'] = results.sum()
results = results.transpose()
results['name'] = results.index
results = results[['name','total tumor']]
resultsPE = results
results = pd.concat([results, resultsPE])

fa = pd.read_csv('/Users/jilke/Scripts/RRBSonCUPs/projects/2022-06_final-scripts_paper/FA-Tpct-corr.csv')
ranking['name'] = ranking.index
rankingPE = ranking
ranking = pd.concat([ranking, rankingPE])
ranking = ranking.merge(results, on='name', how='outer')

fa_class = fa.merge(ranking, on='name', how='outer')
fa_class.to_csv('/Users/jilke/Scripts/RRBSonCUPs/projects/2022-06_final-scripts_paper/FA-Tpct-prediction-corr.csv', index=False)

fa_class


# DMR selections that didn't work
### remove regions with low variance overall
var = samples_df.var(axis=1).sort_values()
var = var[var > 0.01]
samples_df = samples_df[samples_df.index.isin(var.index)]
### DMR selection: lowest intra variance
intra_var = samples_df.groupby(by=samples_df.columns, axis=1).var()
df_lst = []
for ref in list(intra_var):
    df = intra_var[[ref]]
    df = df.nsmallest(1000, ref)
    df_lst.append(df)
df = reduce(lambda x,y: pd.merge(x,y, right_index=True, left_index=True, how='outer'), df_lst)
samples_df = samples_df[samples_df.index.isin(df.index)]

HIGI = pd.read_csv('./methAtlas/new-clusters_test/2022-07-08_HIGI_30x_DSS-DMRs_methAtlas_deconv_output.csv', index_col=0)
HIGI_FFPE = HIGI.filter(regex = 'CUP')
HIGI_FFPE = HIGI_FFPE.drop(columns = list(HIGI_FFPE.filter(regex = 'cf')))
HIGI_FFPE.to_csv('./methAtlas/new-clusters_test/2022-07-08_HIGI-FFPE_30x_DSS-DMRs_methAtlas_deconv_output.csv')
HIGI_plasma = HIGI.filter(regex = 'cf')
HIGI_plasma.drop(columns = ["cfPAAD 1","cfPAAD 2","cfCUP 1"], inplace = True)
HIGI_plasma.to_csv('./methAtlas/new-clusters_test/2022-07-08_HIGI-plasma_30x_DSS-DMRs_methAtlas_deconv_output.csv')

HIGI_plasma
