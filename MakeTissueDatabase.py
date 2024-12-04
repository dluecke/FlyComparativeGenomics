#!/usr/bin/python3

# MakeTissueDatabase.py a script to combine LIMS and cluster information
#  to track RNA samples by species, sex, tissue, rep, and status 
# Takes the LIMS output for Tissue, RNA, and Library as CSV (preremove special chars)
# Takes a FOFN of paths to RNAseq fastq files

# REQUIRES: pandas

# Written to run on Ceres cluster using environment /project/vpgru/software/python_envs/pandas/
# SHELL COMMANDS BEFORE RUNNING SCRIPT:
# $ ml python_3
# $ source /project/vpgru/software/python_envs/pandas/bin/activate

import pandas as pd
import numpy as np
import os

TISSUE_FILE = "AllTissue-export_elabinventory_20241017185222.csv"
RNA_FILE = "AllRNA-export_elabinventory_20241017185222.csv"
LIBRARY_FILE = "AllLibrary-export_elabinventory_20241017185222.csv"
FASTQ_FILE = "RNAfastqs-test.fofn"

TISSUE = pd.read_csv(TISSUE_FILE, dtype=str,
                     usecols=['SampleID', 'Sample Name', 
                              'Description', 'Notes', 'sampleMetaID', 
                              'Species', 'Tissue Type', 'Sex', 'Life Stage'])
TISSUE = TISSUE.add_suffix('_tis')

RNA = pd.read_csv(RNA_FILE, dtype=str,
                  usecols=['SampleID', 'Sample Name', 'Parent', 'Parent sampleID', 
                              'Description', 'Notes', 'sampleMetaID', 
                              'Species', 'Tissue', 'Sex'])
RNA = RNA.add_suffix('_rna')

LIBRARY = pd.read_csv(LIBRARY_FILE, dtype=str,
                      usecols=['SampleID', 'Sample Name', 
                              'Parent', 'Parent sampleID',
                              'Description', 'Notes', 
                              'sampleMetaID', 'Run ID'])
LIBRARY = LIBRARY.add_suffix('_lib')

CATALOG0 = pd.merge(TISSUE, RNA, 
                    left_on='SampleID_tis', 
                    right_on='Parent sampleID_rna', 
                    how='outer')

CATALOG = pd.merge(CATALOG0, LIBRARY, 
                   left_on='SampleID_rna', 
                   right_on='Parent sampleID_lib', 
                   how='outer')

CATALOG['Species'] = CATALOG[['Species_tis', 'Species_rna']].apply(
    lambda x: np.array2string(x.dropna().unique()[:1])[2:-2], axis=1)
CATALOG['Sex'] = CATALOG[['Sex_tis', 'Sex_rna']].apply(
    lambda x: np.array2string(x.dropna().unique()[:1])[2:-2], axis=1)
CATALOG['Tissue'] = CATALOG[['Tissue Type_tis', 'Tissue_rna']].apply(
    lambda x: np.array2string(x.dropna().unique()[:1])[2:-2], axis=1)

FASTQ = pd.read_csv("RNAfastqs-test.fofn", names=['Path'])
FASTQ['File'] = FASTQ['Path'].apply(os.path.basename)
FASTQ['Sample'] = FASTQ['File'].str.replace(
    "_S[0-9]+_", "_", regex=True).replace(
    "_R[1-2]_", "_", regex=True).replace(
    "_001", "", regex=True).replace(
    ".f[ast]*q.gz", "", regex=True).replace(
    "_", "", regex=True)

FASTQ_SAMPLES = pd.DataFrame({
    'Files': FASTQ.groupby('Sample')['File'].apply(list),
    'Paths': FASTQ.groupby('Sample')['Path'].apply(list) })

CATALOG['Seq File'] = CATALOG['Sample Name_lib'].str.replace(
    " |_|,", "", regex=True).apply(
        lambda x: FASTQ_SAMPLES.loc[ 
            FASTQ_SAMPLES.index.str.contains(str(x))]['Files'].tolist() )

CATALOG[['Species', 'Sex', 'Tissue', 'Sample Name_tis', 'Sample Name_rna', 'Sample Name_lib', 'Seq File']].to_csv("TissueCatalog01.tsv", sep="\t", index=False)
