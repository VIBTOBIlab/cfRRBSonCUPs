# cfRRBSonCUPs
CUP diagnosis with FFPE references

1. DMR pipeline
- Use the DMR pipeline on https://github.ugent.be/DePreterLab/DMRpipeline to find DMRs between references.
- The input is a folder with .cov files (output from Bismark methylation calling) and a list with sample ID and group names.
- The output is a folder with .csv files. Each .csv file is a list of DMRs between 2 groups. Here, the "DSS/groups" folder contain all DMRs between the grouped entities (BRCAc, HIGI, LUAD, SCC, DLBL, MESO, OCVA, SCLC, PRAD, SKCM and healthy) and the "DSS/subtypes" folder contains the different subtypes (BRCA and BRCAl; CHOL and GEJC; CSCC, ESCC and LUSC)

2. Make ESR cluster files
- Use the make_ESR-list.py script to select the top 200 ESRs for each entity (=DMRs that can differentiate this entity from all - or most - other entities)
- The input is the output folders from the previous step (DSS/groups and DSS/subtypes) and a list with sample ID and group names (here: refList.csv and refList_grouped.csv)
- The output is a .tsv file which lists all ESRs of all entities. It has the following structure (from bedtools multiinter):

    * chromosome
    * 0-based start coordinate of the ESR.
    * 1-based end coordinate of the ESR.
    * The number of (tumor) entities for which this ESR is specific (usually only 1 ofcourse)
    * The list of (tumor) entity numbers (by order on the command line) for which this ESR is specific
    * Columns reflecting for each (tumor) entity if this ESR is specific (1) or not specific (0)

All output files are in the cluster_files folder.

3. Cluster samples
- Use the cluster_samples.py script to cluster the CpGs from the sample and reference files to the selected ESRs.
- usage: cluster_samples.py [-h] [-i INPUT] [-c CLUSTERFILE]
                          [-o OUTPUT] [--NGS]
- The input is a directory with all sample or reference .cov files
- The clusterfile is the output file from the previous step (in the cluster_files folder). A different clusterfile is used depending on the application (general or subclassification)
- The output is a directory of choice (here clustered_samples/). The script produces 2 files: samples_clustered.csv (containing all clustered samples) and metadata_clusters.csv (using the same sample names as header, but columns contain 'count', the amount of CpGs per cluster and 'depth', the amount of reads (coverage) per cluster)
- Here the --NGS flag should be used.

4. Prep for methAtlas

- Use the prep_methAtlas.py script to prep the references and the samples for deconvolution with methAtlas.
- The input are the clustered samples + metadata from the previous step
- The output is a dataframe in the format required for methatlas
- The references and the samples should be prepped seperately
