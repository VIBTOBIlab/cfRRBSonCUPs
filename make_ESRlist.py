# import libraries
import sys
import glob
import os
import pandas as pd
pd.set_option('display.max_rows',200)
pd.set_option('display.max_columns',None)
from functools import reduce
import pybedtools
from pybedtools import BedTool
x = BedTool()
pybedtools.cleanup()
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# define functions
def check_chr_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False
    except TypeError:
        return False

def calculate_medians(df):
    df = df[df['trueDMR'] == True]
    # problem if > 0.5 of tumors are not in comparison anymore
    df['abs_median_pval'] = abs(df.filter(regex='pval')).median(axis=1)
    df['median_real_diff'] = df.filter(regex='real_diff').median(axis=1)
    df['abs_median_real_diff'] = abs(df['median_real_diff'])
    df.sort_values(by=['abs_median_pval'], ascending=False, inplace=True) #now: DSS, in case of DMRfinder: df.sort_values(by=['median_pval','abs_median_real_diff'], ascending=[True,False], inplace=True)
    return df

# Input files
DMRfolder = "./DSS/groups/"
nameList = "./refList_grouped.csv"

# list tumor types
names = pd.read_csv(nameList)
names['name'] = names['name'].map(lambda x: x.split(' ')[0])
tumor_types = names['name'].unique()
tumor_types = list(tumor_types)
## amount of comparisons
l = len(tumor_types)-1
tumor_types

# make a dict with the intersection of all comparisons per tumor type
df_dict = {}
for tumor in tumor_types:
    pybedtools.cleanup()
    files = glob.glob(DMRfolder + 'results*' + tumor + '*.csv')
    #if tumor == 'BRCA':
    #    files = [file for file in files if not 'BRCAl' in file] + [DMRfolder + 'resultsBRCA-BRCAl.csv']
    BED_list = []
    names = []
    print(f'for {tumor} these files are used: {files}')
    for file in files:
        ## make list with column names
        file_name = os.path.splitext(os.path.basename(file))[0]
        ## include if you don't want to include certain tumor types:
        #if any([tumor in file_name for tumor in removed_types]) == False:
        name = file_name.lstrip('results')
        names.append(name)
        ## make list with bedfiles
        DMR_df = pd.read_csv(file, sep='\t', dtype={'start':'int','end':'int'})
        DMR_BED = BedTool.from_dataframe(DMR_df).sort()
        BED_list.append(DMR_BED)
    ## take intersect of all files with bedtools
    DMR_intersect = x.multiinter(i=[b.fn for b in BED_list], header=True)
    df = pd.read_table(DMR_intersect.fn)
    ## rename columns
    df1 = df.drop(df.filter(regex='.tmp').columns, axis=1)
    df2 = df.filter(regex='.tmp')
    df2.columns=names
    df = df1.join(df2)
    df_dict[tumor] = df

# build DMR file with 200 DMRs per entity
DMRs_dict = {}
DMRcount_df = pd.DataFrame()
true_def_pct = pd.DataFrame(columns = tumor_types)
true_def_cnt = pd.DataFrame(columns = tumor_types)
for tumor in df_dict:
    tumor_dict = {}
    t = 0
    for i in range(l, 5, -1):
        df = df_dict[tumor]
# select the DMRs that all or most comparisons have in common per tumor type
        common_DMRs = df[df['num'] == i]
        if len(common_DMRs) > 0:
            BED = BedTool.from_dataframe(common_DMRs).sort()
# check if common DMRs within 1 tumor type vs others are always hypo- or hypermethylated
##first make list with all annotated DMRs per tumor vs other tumors
            files = glob.glob(DMRfolder + 'results*' + tumor + '*.csv')
            #if tumor == 'BRCA':
            #    files = files = [file for file in files if not 'BRCAl' in file] + [DMRfolder + 'resultsBRCA-BRCAl.csv']
            names = []
            df_list = []
            for file in files:
                ### make list with column names
                file_name = os.path.splitext(os.path.basename(file))[0]
                name = file_name.lstrip('results')
                names.append(name)
                #### for DMRfinder results
#                original_df = pd.read_csv(file, sep='\t')
                #### for DSS results (make input same as DMRfinder, this means incorrectly making areastat = pval)
                original_df = pd.read_csv(file, sep='\t', usecols=['chr','start','end','nCG','meanMethy1','meanMethy2','diff.Methy','areaStat'], dtype={'start':'int','end':'int'})
                T1 = name.split('-')[0]
                T2 = name.split('-')[1]
                orig_names = ['chr','start','end','CpG_'+T1+'-'+T2,T1+':mu_'+T1+'-'+T2,T2+':mu_'+T1+'-'+T2,T1+'->'+T2+':diff',T1+'->'+T2+':pval']
                ### make list with anotated dataframes
                original_BED = BedTool.from_dataframe(original_df).sort()
                intrsct_BED = original_BED.intersect(BED)
                intrsct_df = pd.read_table(intrsct_BED.fn, names=orig_names)
                ### make sure right difference is used (tumor first in equation)
                tmp = intrsct_df.drop(columns=[tumor+':mu_'+T1+'-'+T2])
                tmp['otherT'] = tmp.filter(regex=':mu')
                intrsct_df['real_diff_'+T1+'-'+T2] = intrsct_df[tumor+':mu_'+T1+'-'+T2] - tmp['otherT']
                df_list.append(intrsct_df)
# check per row if diff column always neg or always pos
        ### first remove DMRs on X chromosome or contigs
            lst_noX = []
            for df in df_list:
                mask = df['chr'].map(check_chr_float)
                df = df.loc[mask]
                df['chr'] = df['chr'].astype(int)
                lst_noX.append(df)
            df = reduce(lambda x,y: pd.merge(x,y,on=['chr','start','end'], how='outer'), lst_noX)
            # this line was added because sometimes chr, start, end is not start of dataframe anymore after reducing
            lst = [i for i in list(df) if i not in ['chr','start','end']]
            df = df[['chr','start','end']+lst]
            ### calculate mean and median count and check if difference is always > 0 or always < 0
            df['mean_count'] = df.filter(regex='CpG').mean(axis=1)
            df['median_count'] = df.filter(regex='CpG').median(axis=1)
            df.fillna(value=0, inplace=True)
            df['trueDMR'] = (df.filter(regex='real_diff') >= 0).all(axis=1) | (df.filter(regex='real_diff') <= 0).all(axis=1)
            ### report per tumor how many "true"(=all pos or all neg) or "false" DMRs there are
            df_false = df[df['trueDMR'] == False]
            if i == l:
                df_false.to_csv(DMRfolder+'falseDMRs/'+tumor+'.csv', sep='\t', header=True, index = False)
            df_true = df[df['trueDMR'] == True]
            print(f"{tumor} has {len(df)} defining autosomal DMRs against {i} tumor types with {len(df_true)} true DMRs and {len(df_false)} false")
            tumor_dict[i] = df
# check for each entity which other entity needs to be left out to get more defining DMRs ('true' DMRs)
            if len(df_true) > 0:
                t += 1
            print(t)
            if t >= 1:
                true = df[df['trueDMR'] == True]
                true = true.filter(regex='real_diff')
                true.columns = true.columns.str.lstrip('real_diff_')
                lst = true.columns.str.split('-')
                new_lst = []
                for j in lst:
                    j.remove(tumor)
                    newID = j[0]
                    new_lst.append(newID)
                true.columns = new_lst
                tot = len(true)
                # calculate percentage
                zero_df = pd.DataFrame(columns = true.columns, index = [tumor + '_' + str(i)])
                for col in true.columns:
                    count = (true[col] == 0).sum()
                    perc = count/tot*100
                    zero_df.loc[tumor + '_' + str(i),col] = round(perc,1)
                true_def_pct = pd.concat([true_def_pct,zero_df])
                # calculate count
                zero_df = pd.DataFrame(columns = true.columns, index = [tumor + '_' + str(i)])
                for col in true.columns:
                    zero_df.loc[tumor + '_' + str(i),col] = (true[col] == 0).sum()
                true_def_cnt = pd.concat([true_def_cnt,zero_df])
# make new DMR file using 200 DMRs from each tumor type
## select DMRs for each tumor type by first picking all DMRs that are in common in the comparisons with all other tumor types
## if not enough, select DMRs that all but 1 have in common
## if not enough, etc...
    i = 0
    check = False
    for key in tumor_dict:
        if i == 0:
            df = tumor_dict[key]
            df = calculate_medians(df)
            if len(df) >= 200:
                df = df.head(200)
        elif i > 0:
            dfi = tumor_dict[key]
            dfi = calculate_medians(dfi)
            # merge new DMRs with existing
            lst = list(dfi)
            BED = BedTool.from_dataframe(dfi).sort()
            BEDm = BED.merge(c=list(range(4,len(dfi.columns)+1)), o=(['sum']+['mean']*5)*int((len(dfi.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
            dfi = pd.read_table(BEDm.fn, names=list(dfi))
            dfh = dfi.head(200-len(df))
            df = pd.concat([df, dfh], ignore_index=True)
            if len(df) == 200:
                check = True
        if len(df) > 0:
            i += 1
            # merge new DMRs with existing
            lst = list(df)
            BED = BedTool.from_dataframe(df).sort()
            BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=(['sum']+['mean']*5)*int((len(df.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
            df = pd.read_table(BEDm.fn, names=list(df))
            DMRcount_df.loc[key, tumor] = len(df)
            # after merging, sometimes df that had 200 DMRs is again to short, so then next DMRs should be taken from same dfi dataframe
            if (check == True) & (len(df) < 200):
                while len(dfi) > 0:
                    dfi = dfi.iloc[len(dfh):]
                    dfh2 = dfi.head(200-len(df))
                    dfh = pd.concat([dfh, dfh2], ignore_index=True)
                    df = pd.concat([df, dfh2], ignore_index=True)
                    BED = BedTool.from_dataframe(df).sort()
                    BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=(['sum']+['mean']*5)*int((len(df.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
                    df = pd.read_table(BEDm.fn, names=list(df))
                    DMRcount_df.loc[key, tumor] = len(df)
                    if len(df) == 200:
                        break
        if len(df) == 200:
            break
    DMRs_dict[tumor] = df

## make intersect of all DMRs
BED_list = []
names = []
for tumor in DMRs_dict:
    ### make list with bedfiles
    DMR_df = DMRs_dict[tumor]
    DMR_BED = BedTool.from_dataframe(DMR_df).sort()
    BED_list.append(DMR_BED)
    ### make list with column names
    name = str(tumor)
    names.append(name)
### intersect
DMR_intersect = x.multiinter(i=[b.fn for b in BED_list], header=True)
df = pd.read_table(DMR_intersect.fn)
## rename columns
df1 = df.drop(df.filter(regex='.tmp').columns, axis=1)
df2 = df.filter(regex='.tmp')
df2.columns=names
df = df1.join(df2)
## merge adjacent DMRs
BED = BedTool.from_dataframe(df).sort()
BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=['sum','collapse']+['sum']*(len(df.columns)-5))
df = pd.read_table(BEDm.fn, names=list(df))
## rename cols to be used with cluster_samples.py script
df.rename(columns={'chrom':'chr'},inplace=True)
df.rename(columns={'end':'stop'},inplace=True)
df['clusterID'] = df.index
## check amount of unique/overlapping DMRs
df['distinct'] = df['list'].map(lambda x: set(x.split(",")))
df['num_distinct'] = df['distinct'].map(lambda x: len(x))
df.to_csv('./cluster_files/DSS-DMRs-12G_clusters.tsv', sep='\t', header=True, index=False)

# make specific DMR sets for hierarchical classification

## Input files
DMRfolder = "./DSS/subtypes/"
nameList = "./refList.csv"

## BRCA - BRCAl = BRCAc
### make a dict with the intersection of all comparisons per tumor type
BRCAc_dict = {}
for tumor in ['BRCA','BRCAl','healthy']:
    pybedtools.cleanup()
    files = glob.glob(DMRfolder + 'results*' + tumor + '*.csv')
    if tumor == 'BRCA':
        files = [file for file in files if not 'BRCAl' in file] + [DMRfolder + 'resultsBRCA-BRCAl.csv']
    BED_list = []
    names = []
    print(f'for {tumor} these files are used: {files}')
    for file in files:
        ## make list with column names
        file_name = os.path.splitext(os.path.basename(file))[0]
        ## include if you don't want to include certain tumor types:
        #if any([tumor in file_name for tumor in removed_types]) == False:
        name = file_name.lstrip('results')
        names.append(name)
        ## make list with bedfiles
        DMR_df = pd.read_csv(file, sep='\t', dtype={'start':'int','end':'int'})
        DMR_BED = BedTool.from_dataframe(DMR_df).sort()
        BED_list.append(DMR_BED)
    ## take intersect of all files with bedtools
    DMR_intersect = x.multiinter(i=[b.fn for b in BED_list], header=True)
    df = pd.read_table(DMR_intersect.fn)
    ## rename columns
    df1 = df.drop(df.filter(regex='.tmp').columns, axis=1)
    df2 = df.filter(regex='.tmp')
    df2.columns=names
    df = df1.join(df2)
    BRCAc_dict[tumor] = df

### build DMR file with 200 DMRs per entity
DMRs_BRCAc = {}
DMRcount_df = pd.DataFrame()
true_def_df = pd.DataFrame(columns = tumor_types)
for tumor in BRCAc_dict:
    tumor_dict = {}
    t = 0
    for i in range(l, 0, -1):
        df = BRCAc_dict[tumor]
# select the DMRs that all or most comparisons have in common per tumor type
        common_DMRs = df[df['num'] == i]
        if tumor == 'BRCA':
            common_DMRs = common_DMRs[(common_DMRs['BRCA-healthy'] == 1) & (common_DMRs['BRCA-BRCAl'] == 1)]
        if tumor == 'BRCAl':
            common_DMRs = common_DMRs[(common_DMRs['BRCAl-healthy'] == 1) & (common_DMRs['BRCA-BRCAl'] == 1)]
        if tumor == 'healthy':
            common_DMRs = common_DMRs[(common_DMRs['BRCA-healthy'] == 1) & (common_DMRs['BRCAl-healthy'] == 1)]
        if len(common_DMRs) > 0:
            BED = BedTool.from_dataframe(common_DMRs).sort()
# check if common DMRs within 1 tumor type vs others are always hypo- or hypermethylated
##first make list with all annotated DMRs per tumor vs other tumors
            files = glob.glob(DMRfolder + 'results*' + tumor + '*.csv')
            if tumor == 'BRCA':
                files = [file for file in files if not 'BRCAl' in file] + [DMRfolder + 'resultsBRCA-BRCAl.csv']
            names = []
            df_list = []
            for file in files:
                ### make list with column names
                file_name = os.path.splitext(os.path.basename(file))[0]
                name = file_name.lstrip('results')
                names.append(name)
                #### for DMRfinder results
#                original_df = pd.read_csv(file, sep='\t')
                #### for DSS results (make input same as DMRfinder, this means incorrectly making areastat = pval)
                original_df = pd.read_csv(file, sep='\t', usecols=['chr','start','end','nCG','meanMethy1','meanMethy2','diff.Methy','areaStat'], dtype={'start':'int','end':'int'})
                T1 = name.split('-')[0]
                T2 = name.split('-')[1]
                orig_names = ['chr','start','end','CpG_'+T1+'-'+T2,T1+':mu_'+T1+'-'+T2,T2+':mu_'+T1+'-'+T2,T1+'->'+T2+':diff',T1+'->'+T2+':pval']
                ### make list with anotated dataframes
                original_BED = BedTool.from_dataframe(original_df).sort()
                intrsct_BED = original_BED.intersect(BED)
                intrsct_df = pd.read_table(intrsct_BED.fn, names=orig_names)
                ### make sure right difference is used (tumor first in equation)
                tmp = intrsct_df.drop(columns=[tumor+':mu_'+T1+'-'+T2])
                tmp['otherT'] = tmp.filter(regex=':mu')
                intrsct_df['real_diff_'+T1+'-'+T2] = intrsct_df[tumor+':mu_'+T1+'-'+T2] - tmp['otherT']
                df_list.append(intrsct_df)
# check per row if diff column always neg or always pos
        ### first remove DMRs on X chromosome or contigs
            lst_noX = []
            for df in df_list:
                mask = df['chr'].map(check_chr_float)
                df = df.loc[mask]
                df['chr'] = df['chr'].astype(int)
                lst_noX.append(df)
            df = reduce(lambda x,y: pd.merge(x,y,on=['chr','start','end'], how='outer'), lst_noX)
            # this line was added because sometimes chr, start, end is not start of dataframe anymore after reducing
            lst = [i for i in list(df) if i not in ['chr','start','end']]
            df = df[['chr','start','end']+lst]
            ### calculate mean and median count and check if difference is always > 0 or always < 0
            df['mean_count'] = df.filter(regex='CpG').mean(axis=1)
            df['median_count'] = df.filter(regex='CpG').median(axis=1)
            df.fillna(value=0, inplace=True)
            df['trueDMR'] = (df.filter(regex='real_diff') >= 0).all(axis=1) | (df.filter(regex='real_diff') <= 0).all(axis=1)
            ### report per tumor how many "true"(=all pos or all neg) or "false" DMRs there are
            df_false = df[df['trueDMR'] == False]
            if i == l:
                df_false.to_csv(DMRfolder+'falseDMRs/'+tumor+'.csv', sep='\t', header=True, index = False)
            df_true = df[df['trueDMR'] == True]
            print(f"{tumor} has {len(df)} defining autosomal DMRs against {i} tumor types with {len(df_true)} true DMRs and {len(df_false)} false")
            tumor_dict[i] = df
# check for each entity which other entity needs to be left out to get more defining DMRs ('true' DMRs)
            if len(df_true) > 0:
                t += 1
            print(t)
            if t == 2:
                true = df[df['trueDMR'] == True]
                true = true.filter(regex='diff')
                true = true.filter(regex='>')
                true.columns = true.columns.str.rstrip(':diff')
                lst = true.columns.str.split('->')
                new_lst = []
                for i in lst:
                    i.remove(tumor)
                    newID = i[0]
                    new_lst.append(newID)
                true.columns = new_lst
                tot = len(true)
                zero_df = pd.DataFrame(columns = true.columns, index = [tumor])
                for col in true.columns:
                    column = true[col]
                    count = (column == 0).sum()
                    perc = count/tot*100
                    zero_df.loc[tumor,col] = round(perc,1)
                true_def_df = pd.concat([true_def_df,zero_df])
            if (len(df_true) > 300) & (t > 1):
                break
# make new DMR file using 200 DMRs from each tumor type
## select DMRs for each tumor type by first picking all DMRs that are in common in the comparisons with all other tumor types
## if not enough, select DMRs that all but 1 have in common
## if not enough, etc...
    i = 0
    check = False
    for key in tumor_dict:
        if i == 0:
            df = tumor_dict[key]
            df = calculate_medians(df)
            if len(df) >= 200:
                df = df.head(200)
                DMRcount_df.loc[key, tumor] = len(df)
            elif (len(df) < 200) & (len(df) > 0):
                DMRcount_df.loc[key, tumor] = len(df)
        elif i > 0:
            dfi = tumor_dict[key]
            dfi = calculate_medians(dfi)
            # 7/07/2022 new addition: add merging step
            lst = list(dfi)
            BED = BedTool.from_dataframe(dfi).sort()
            BEDm = BED.merge(c=list(range(4,len(dfi.columns)+1)), o=(['sum']+['mean']*5)*int((len(dfi.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
            dfi = pd.read_table(BEDm.fn, names=list(dfi))
            dfh = dfi.head(200-len(df))
            DMRcount_df.loc[key, tumor] = len(dfh)
            df = pd.concat([df, dfh], ignore_index=True)
            if len(df) == 200:
                check = True
        if len(df) > 0:
            i += 1
            # 7/07/2022 new addition: add merging step
            lst = list(df)
            BED = BedTool.from_dataframe(df).sort()
            BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=(['sum']+['mean']*5)*int((len(df.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
            df = pd.read_table(BEDm.fn, names=list(df))
            # after merging, sometimes df that had 200 DMRs is again to short, so then next DMRs should be taken from same dfi dataframe
            if (check == True) & (len(df) < 200):
                while len(dfi) > 0:
                    dfi = dfi.iloc[len(dfh):]
                    dfh2 = dfi.head(200-len(df))
                    dfh = pd.concat([dfh, dfh2], ignore_index=True)
                    DMRcount_df.loc[key, tumor] = len(dfh)
                    df = pd.concat([df, dfh2], ignore_index=True)
                    BED = BedTool.from_dataframe(df).sort()
                    BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=(['sum']+['mean']*5)*int((len(df.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
                    df = pd.read_table(BEDm.fn, names=list(df))
                    if len(df) == 200:
                        break
        if len(df) == 200:
            break
    DMRs_BRCAc[tumor] = df

## make intersect of all DMRs
BED_list = []
names = []
for tumor in DMRs_BRCAc:
    ### make list with bedfiles
    DMR_df = DMRs_BRCAc[tumor]
    DMR_BED = BedTool.from_dataframe(DMR_df).sort()
    BED_list.append(DMR_BED)
    ### make list with column names
    name = str(tumor)
    names.append(name)
### intersect
DMR_intersect = x.multiinter(i=[b.fn for b in BED_list], header=True)
df = pd.read_table(DMR_intersect.fn)
## rename columns
df1 = df.drop(df.filter(regex='.tmp').columns, axis=1)
df2 = df.filter(regex='.tmp')
df2.columns=names
df = df1.join(df2)
## merge adjacent DMRs
BED = BedTool.from_dataframe(df).sort()
BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=['sum','collapse']+['sum']*(len(df.columns)-5))
df = pd.read_table(BEDm.fn, names=list(df))
## rename cols to be used with cluster_samples.py script
df.rename(columns={'chrom':'chr'},inplace=True)
df.rename(columns={'end':'stop'},inplace=True)
df['clusterID'] = df.index
## check amount of unique/overlapping DMRs
df['distinct'] = df['list'].map(lambda x: set(x.split(",")))
df['num_distinct'] = df['distinct'].map(lambda x: len(x))
df.to_csv('./cluster_files/DSS-DMRs-BRCAc_clusters.tsv', sep='\t', header=True, index=False)

## PAAD/CHOL - GEJC = HIGI
### make a dict with the intersection of all comparisons per tumor type
HIGI_dict = {}
for tumor in ['healthy','GEJC','CHOL']:
    pybedtools.cleanup()
    files = glob.glob(DMRfolder + 'results*' + tumor + '*.csv')
    BED_list = []
    names = []
    print(f'for {tumor} these files are used: {files}')
    for file in files:
        ## make list with column names
        file_name = os.path.splitext(os.path.basename(file))[0]
        ## include if you don't want to include certain tumor types:
        #if any([tumor in file_name for tumor in removed_types]) == False:
        name = file_name.lstrip('results')
        names.append(name)
        ## make list with bedfiles
        DMR_df = pd.read_csv(file, sep='\t', dtype={'start':'int','end':'int'})
        DMR_BED = BedTool.from_dataframe(DMR_df).sort()
        BED_list.append(DMR_BED)
    ## take intersect of all files with bedtools
    DMR_intersect = x.multiinter(i=[b.fn for b in BED_list], header=True)
    df = pd.read_table(DMR_intersect.fn)
    ## rename columns
    df1 = df.drop(df.filter(regex='.tmp').columns, axis=1)
    df2 = df.filter(regex='.tmp')
    df2.columns=names
    df = df1.join(df2)
    HIGI_dict[tumor] = df

### build DMR file with 200 DMRs per entity
DMRs_HIGI = {}
DMRcount_df = pd.DataFrame()
true_def_df = pd.DataFrame(columns = tumor_types)
for tumor in HIGI_dict:
    tumor_dict = {}
    t = 0
    for i in range(l, 0, -1):
        df = HIGI_dict[tumor]
# select the DMRs that all or most comparisons have in common per tumor type
        common_DMRs = df[df['num'] == i]
        if tumor == 'CHOL':
            common_DMRs = common_DMRs[(common_DMRs['CHOL-healthy'] == 1) & (common_DMRs['CHOL-GEJC'] == 1)]
        if tumor == 'GEJC':
            common_DMRs = common_DMRs[(common_DMRs['GEJC-healthy'] == 1) & (common_DMRs['CHOL-GEJC'] == 1)]
        if tumor == 'healthy':
            common_DMRs = common_DMRs[(common_DMRs['CHOL-healthy'] == 1) & (common_DMRs['GEJC-healthy'] == 1)]
        if len(common_DMRs) > 0:
            BED = BedTool.from_dataframe(common_DMRs).sort()
# check if common DMRs within 1 tumor type vs others are always hypo- or hypermethylated
##first make list with all annotated DMRs per tumor vs other tumors
            files = glob.glob(DMRfolder + 'results*' + tumor + '*.csv')
            names = []
            df_list = []
            for file in files:
                ### make list with column names
                file_name = os.path.splitext(os.path.basename(file))[0]
                name = file_name.lstrip('results')
                names.append(name)
                #### for DMRfinder results
#                original_df = pd.read_csv(file, sep='\t')
                #### for DSS results (make input same as DMRfinder, this means incorrectly making areastat = pval)
                original_df = pd.read_csv(file, sep='\t', usecols=['chr','start','end','nCG','meanMethy1','meanMethy2','diff.Methy','areaStat'], dtype={'start':'int','end':'int'})
                T1 = name.split('-')[0]
                T2 = name.split('-')[1]
                orig_names = ['chr','start','end','CpG_'+T1+'-'+T2,T1+':mu_'+T1+'-'+T2,T2+':mu_'+T1+'-'+T2,T1+'->'+T2+':diff',T1+'->'+T2+':pval']
                ### make list with anotated dataframes
                original_BED = BedTool.from_dataframe(original_df).sort()
                intrsct_BED = original_BED.intersect(BED)
                intrsct_df = pd.read_table(intrsct_BED.fn, names=orig_names)
                ### make sure right difference is used (tumor first in equation)
                tmp = intrsct_df.drop(columns=[tumor+':mu_'+T1+'-'+T2])
                tmp['otherT'] = tmp.filter(regex=':mu')
                intrsct_df['real_diff_'+T1+'-'+T2] = intrsct_df[tumor+':mu_'+T1+'-'+T2] - tmp['otherT']
                df_list.append(intrsct_df)
# check per row if diff column always neg or always pos
        ### first remove DMRs on X chromosome or contigs
            lst_noX = []
            for df in df_list:
                mask = df['chr'].map(check_chr_float)
                df = df.loc[mask]
                df['chr'] = df['chr'].astype(int)
                lst_noX.append(df)
            df = reduce(lambda x,y: pd.merge(x,y,on=['chr','start','end'], how='outer'), lst_noX)
            # this line was added because sometimes chr, start, end is not start of dataframe anymore after reducing
            lst = [i for i in list(df) if i not in ['chr','start','end']]
            df = df[['chr','start','end']+lst]
            ### calculate mean and median count and check if difference is always > 0 or always < 0
            df['mean_count'] = df.filter(regex='CpG').mean(axis=1)
            df['median_count'] = df.filter(regex='CpG').median(axis=1)
            df.fillna(value=0, inplace=True)
            df['trueDMR'] = (df.filter(regex='real_diff') >= 0).all(axis=1) | (df.filter(regex='real_diff') <= 0).all(axis=1)
            ### report per tumor how many "true"(=all pos or all neg) or "false" DMRs there are
            df_false = df[df['trueDMR'] == False]
            if i == l:
                df_false.to_csv(DMRfolder+'falseDMRs/'+tumor+'.csv', sep='\t', header=True, index = False)
            df_true = df[df['trueDMR'] == True]
            print(f"{tumor} has {len(df)} defining autosomal DMRs against {i} tumor types with {len(df_true)} true DMRs and {len(df_false)} false")
            tumor_dict[i] = df
# check for each entity which other entity needs to be left out to get more defining DMRs ('true' DMRs)
            if len(df_true) > 0:
                t += 1
            print(t)
            if t == 2:
                true = df[df['trueDMR'] == True]
                true = true.filter(regex='diff')
                true = true.filter(regex='>')
                true.columns = true.columns.str.rstrip(':diff')
                lst = true.columns.str.split('->')
                new_lst = []
                for i in lst:
                    i.remove(tumor)
                    newID = i[0]
                    new_lst.append(newID)
                true.columns = new_lst
                tot = len(true)
                zero_df = pd.DataFrame(columns = true.columns, index = [tumor])
                for col in true.columns:
                    column = true[col]
                    count = (column == 0).sum()
                    perc = count/tot*100
                    zero_df.loc[tumor,col] = round(perc,1)
                true_def_df = pd.concat([true_def_df,zero_df])
            if (len(df_true) > 300) & (t > 1):
                break
# make new DMR file using 200 DMRs from each tumor type
## select DMRs for each tumor type by first picking all DMRs that are in common in the comparisons with all other tumor types
## if not enough, select DMRs that all but 1 have in common
## if not enough, etc...
    i = 0
    check = False
    for key in tumor_dict:
        if i == 0:
            df = tumor_dict[key]
            df = calculate_medians(df)
            if len(df) >= 200:
                df = df.head(200)
                DMRcount_df.loc[key, tumor] = len(df)
            elif (len(df) < 200) & (len(df) > 0):
                DMRcount_df.loc[key, tumor] = len(df)
        elif i > 0:
            dfi = tumor_dict[key]
            dfi = calculate_medians(dfi)
            # 7/07/2022 new addition: add merging step
            lst = list(dfi)
            BED = BedTool.from_dataframe(dfi).sort()
            BEDm = BED.merge(c=list(range(4,len(dfi.columns)+1)), o=(['sum']+['mean']*5)*int((len(dfi.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
            dfi = pd.read_table(BEDm.fn, names=list(dfi))
            dfh = dfi.head(200-len(df))
            DMRcount_df.loc[key, tumor] = len(dfh)
            df = pd.concat([df, dfh], ignore_index=True)
            if len(df) == 200:
                check = True
        if len(df) > 0:
            i += 1
            # 7/07/2022 new addition: add merging step
            lst = list(df)
            BED = BedTool.from_dataframe(df).sort()
            BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=(['sum']+['mean']*5)*int((len(df.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
            df = pd.read_table(BEDm.fn, names=list(df))
            # after merging, sometimes df that had 200 DMRs is again to short, so then next DMRs should be taken from same dfi dataframe
            if (check == True) & (len(df) < 200):
                while len(dfi) > 0:
                    dfi = dfi.iloc[len(dfh):]
                    dfh2 = dfi.head(200-len(df))
                    dfh = pd.concat([dfh, dfh2], ignore_index=True)
                    DMRcount_df.loc[key, tumor] = len(dfh)
                    df = pd.concat([df, dfh2], ignore_index=True)
                    BED = BedTool.from_dataframe(df).sort()
                    BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=(['sum']+['mean']*5)*int((len(df.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
                    df = pd.read_table(BEDm.fn, names=list(df))
                    if len(df) == 200:
                        break
        if len(df) == 200:
            break
    DMRs_HIGI[tumor] = df

## make intersect of all DMRs
BED_list = []
names = []
for tumor in DMRs_HIGI:
    ### make list with bedfiles
    DMR_df = DMRs_HIGI[tumor]
    DMR_BED = BedTool.from_dataframe(DMR_df).sort()
    BED_list.append(DMR_BED)
    ### make list with column names
    name = str(tumor)
    names.append(name)
### intersect
DMR_intersect = x.multiinter(i=[b.fn for b in BED_list], header=True)
df = pd.read_table(DMR_intersect.fn)
## rename columns
df1 = df.drop(df.filter(regex='.tmp').columns, axis=1)
df2 = df.filter(regex='.tmp')
df2.columns=names
df = df1.join(df2)
## merge adjacent DMRs
BED = BedTool.from_dataframe(df).sort()
BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=['sum','collapse']+['sum']*(len(df.columns)-5))
df = pd.read_table(BEDm.fn, names=list(df))
## rename cols to be used with cluster_samples.py script
df.rename(columns={'chrom':'chr'},inplace=True)
df.rename(columns={'end':'stop'},inplace=True)
df['clusterID'] = df.index
## check amount of unique/overlapping DMRs
df['distinct'] = df['list'].map(lambda x: set(x.split(",")))
df['num_distinct'] = df['distinct'].map(lambda x: len(x))
df.to_csv('./cluster_files/DSS-DMRs-HIGI_clusters.tsv', sep='\t', header=True, index=False)

## CSCC - ESCC - LUSC = SCC
SCC_dict = {}
for tumor in ['CSCC','ESCC','LUSC','healthy']:
    pybedtools.cleanup()
    files = glob.glob(DMRfolder + 'results*' + tumor + '*.csv')
    BED_list = []
    names = []
    print(f'for {tumor} these files are used: {files}')
    for file in files:
        ## make list with column names
        file_name = os.path.splitext(os.path.basename(file))[0]
        ## include if you don't want to include certain tumor types:
        #if any([tumor in file_name for tumor in removed_types]) == False:
        name = file_name.lstrip('results')
        names.append(name)
        ## make list with bedfiles
        DMR_df = pd.read_csv(file, sep='\t', dtype={'start':'int','end':'int'})
        DMR_BED = BedTool.from_dataframe(DMR_df).sort()
        BED_list.append(DMR_BED)
    ## take intersect of all files with bedtools
    DMR_intersect = x.multiinter(i=[b.fn for b in BED_list], header=True)
    df = pd.read_table(DMR_intersect.fn)
    ## rename columns
    df1 = df.drop(df.filter(regex='.tmp').columns, axis=1)
    df2 = df.filter(regex='.tmp')
    df2.columns=names
    df = df1.join(df2)
    SCC_dict[tumor] = df

# build DMR file with 200 DMRs per entity
DMRs_SCC = {}
DMRcount_df = pd.DataFrame()
true_def_df = pd.DataFrame(columns = tumor_types)
for tumor in SCC_dict:
    tumor_dict = {}
    t = 0
    for i in range(3, 0, -1):
        df = SCC_dict[tumor]
# select the DMRs that all or most comparisons have in common per tumor type
        common_DMRs = df[df['num'] == i]
        if len(common_DMRs) > 0:
            BED = BedTool.from_dataframe(common_DMRs).sort()
# check if common DMRs within 1 tumor type vs others are always hypo- or hypermethylated
##first make list with all annotated DMRs per tumor vs other tumors
            files = glob.glob(DMRfolder + 'results*' + tumor + '*.csv')
            names = []
            df_list = []
            for file in files:
                ### make list with column names
                file_name = os.path.splitext(os.path.basename(file))[0]
                name = file_name.lstrip('results')
                names.append(name)
                #### for DMRfinder results
#                original_df = pd.read_csv(file, sep='\t')
                #### for DSS results (make input same as DMRfinder, this means incorrectly making areastat = pval)
                original_df = pd.read_csv(file, sep='\t', usecols=['chr','start','end','nCG','meanMethy1','meanMethy2','diff.Methy','areaStat'], dtype={'start':'int','end':'int'})
                T1 = name.split('-')[0]
                T2 = name.split('-')[1]
                orig_names = ['chr','start','end','CpG_'+T1+'-'+T2,T1+':mu_'+T1+'-'+T2,T2+':mu_'+T1+'-'+T2,T1+'->'+T2+':diff',T1+'->'+T2+':pval']
                ### make list with anotated dataframes
                original_BED = BedTool.from_dataframe(original_df).sort()
                intrsct_BED = original_BED.intersect(BED)
                intrsct_df = pd.read_table(intrsct_BED.fn, names=orig_names)
                ### make sure right difference is used (tumor first in equation)
                tmp = intrsct_df.drop(columns=[tumor+':mu_'+T1+'-'+T2])
                tmp['otherT'] = tmp.filter(regex=':mu')
                intrsct_df['real_diff_'+T1+'-'+T2] = intrsct_df[tumor+':mu_'+T1+'-'+T2] - tmp['otherT']
                df_list.append(intrsct_df)
# check per row if diff column always neg or always pos
        ### first remove DMRs on X chromosome or contigs
            lst_noX = []
            for df in df_list:
                mask = df['chr'].map(check_chr_float)
                df = df.loc[mask]
                df['chr'] = df['chr'].astype(int)
                lst_noX.append(df)
            df = reduce(lambda x,y: pd.merge(x,y,on=['chr','start','end'], how='outer'), lst_noX)
            # this line was added because sometimes chr, start, end is not start of dataframe anymore after reducing
            lst = [i for i in list(df) if i not in ['chr','start','end']]
            df = df[['chr','start','end']+lst]
            ### calculate mean and median count and check if difference is always > 0 or always < 0
            df['mean_count'] = df.filter(regex='CpG').mean(axis=1)
            df['median_count'] = df.filter(regex='CpG').median(axis=1)
            df.fillna(value=0, inplace=True)
            df['trueDMR'] = (df.filter(regex='real_diff') >= 0).all(axis=1) | (df.filter(regex='real_diff') <= 0).all(axis=1)
            ### report per tumor how many "true"(=all pos or all neg) or "false" DMRs there are
            df_false = df[df['trueDMR'] == False]
            if i == l:
                df_false.to_csv(DMRfolder+'falseDMRs/'+tumor+'.csv', sep='\t', header=True, index = False)
            df_true = df[df['trueDMR'] == True]
            print(f"{tumor} has {len(df)} defining autosomal DMRs against {i} tumor types with {len(df_true)} true DMRs and {len(df_false)} false")
            tumor_dict[i] = df
# check for each entity which other entity needs to be left out to get more defining DMRs ('true' DMRs)
            if len(df_true) > 0:
                t += 1
            print(t)
            if t == 2:
                true = df[df['trueDMR'] == True]
                true = true.filter(regex='diff')
                true = true.filter(regex='>')
                true.columns = true.columns.str.rstrip(':diff')
                lst = true.columns.str.split('->')
                new_lst = []
                for i in lst:
                    i.remove(tumor)
                    newID = i[0]
                    new_lst.append(newID)
                true.columns = new_lst
                tot = len(true)
                zero_df = pd.DataFrame(columns = true.columns, index = [tumor])
                for col in true.columns:
                    column = true[col]
                    count = (column == 0).sum()
                    perc = count/tot*100
                    zero_df.loc[tumor,col] = round(perc,1)
                true_def_df = pd.concat([true_def_df,zero_df])
            if (len(df_true) > 300) & (t > 1):
                break
# make new DMR file using 200 DMRs from each tumor type
## select DMRs for each tumor type by first picking all DMRs that are in common in the comparisons with all other tumor types
## if not enough, select DMRs that all but 1 have in common
## if not enough, etc...
    i = 0
    check = False
    for key in tumor_dict:
        if i == 0:
            df = tumor_dict[key]
            df = calculate_medians(df)
            if len(df) >= 200:
                df = df.head(200)
                DMRcount_df.loc[key, tumor] = len(df)
            elif (len(df) < 200) & (len(df) > 0):
                DMRcount_df.loc[key, tumor] = len(df)
        elif i > 0:
            dfi = tumor_dict[key]
            dfi = calculate_medians(dfi)
            # 7/07/2022 new addition: add merging step
            lst = list(dfi)
            BED = BedTool.from_dataframe(dfi).sort()
            BEDm = BED.merge(c=list(range(4,len(dfi.columns)+1)), o=(['sum']+['mean']*5)*int((len(dfi.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
            dfi = pd.read_table(BEDm.fn, names=list(dfi))
            dfh = dfi.head(200-len(df))
            DMRcount_df.loc[key, tumor] = len(dfh)
            df = pd.concat([df, dfh], ignore_index=True)
            if len(df) == 200:
                check = True
        if len(df) > 0:
            i += 1
            # 7/07/2022 new addition: add merging step
            lst = list(df)
            BED = BedTool.from_dataframe(df).sort()
            BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=(['sum']+['mean']*5)*int((len(df.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
            df = pd.read_table(BEDm.fn, names=list(df))
            # after merging, sometimes df that had 200 DMRs is again to short, so then next DMRs should be taken from same dfi dataframe
            if (check == True) & (len(df) < 200):
                while len(dfi) > 0:
                    dfi = dfi.iloc[len(dfh):]
                    dfh2 = dfi.head(200-len(df))
                    dfh = pd.concat([dfh, dfh2], ignore_index=True)
                    DMRcount_df.loc[key, tumor] = len(dfh)
                    df = pd.concat([df, dfh2], ignore_index=True)
                    BED = BedTool.from_dataframe(df).sort()
                    BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=(['sum']+['mean']*5)*int((len(df.columns)-9)/6)+['sum']*2+['distinct']+['mean']*3)
                    df = pd.read_table(BEDm.fn, names=list(df))
                    if len(df) == 200:
                        break
        if len(df) == 200:
            break
    DMRs_SCC[tumor] = df

## make intersect of all DMRs
BED_list = []
names = []
for tumor in DMRs_SCC:
    ### make list with bedfiles
    DMR_df = DMRs_SCC[tumor]
    DMR_BED = BedTool.from_dataframe(DMR_df).sort()
    BED_list.append(DMR_BED)
    ### make list with column names
    name = str(tumor)
    names.append(name)
### intersect
DMR_intersect = x.multiinter(i=[b.fn for b in BED_list], header=True)
df = pd.read_table(DMR_intersect.fn)
## rename columns
df1 = df.drop(df.filter(regex='.tmp').columns, axis=1)
df2 = df.filter(regex='.tmp')
df2.columns=names
df = df1.join(df2)
## merge adjacent DMRs
BED = BedTool.from_dataframe(df).sort()
BEDm = BED.merge(c=list(range(4,len(df.columns)+1)), o=['sum','collapse']+['sum']*(len(df.columns)-5))
df = pd.read_table(BEDm.fn, names=list(df))
## rename cols to be used with cluster_samples.py script
df.rename(columns={'chrom':'chr'},inplace=True)
df.rename(columns={'end':'stop'},inplace=True)
df['clusterID'] = df.index
## check amount of unique/overlapping DMRs
df['distinct'] = df['list'].map(lambda x: set(x.split(",")))
df['num_distinct'] = df['distinct'].map(lambda x: len(x))
len(df.loc[df['num_distinct'] > 1])
#df = df.loc[df['num_distinct'] == 1]
df.to_csv('./cluster_files/DSS-DMRs-SCC_clusters.tsv', sep='\t', header=True, index=False)
