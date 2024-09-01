#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os

import matplotlib
from matplotlib import rcParams


# sets backend of Matplotlib to "Agg"; a non-interactive backend for saving figures to files w/out displaying on screen
matplotlib.use('Agg')
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

############################################################
################### DIRECTORY STUFF ########################
############################################################

partial_path = '.../CellphoneDB_analysis/cpdb_outputs/'
txt_file = '/statistical_analysis_significant_means_.txt'

AD_Skin_GSE147424_Healthy = partial_path + 'AD_Skin_GSE147424/Healthy' + txt_file
AD_Skin_GSE147424_Lesional = partial_path + 'AD_Skin_GSE147424/Lesional' + txt_file
AD_Skin_GSE147424_Non_Lesional = partial_path + 'AD_Skin_GSE147424/Non_Lesional' + txt_file

AD_Skin_GSE153760_HC = partial_path + 'AD_Skin_GSE153760/HC' + txt_file
AD_Skin_GSE153760_AD = partial_path + 'AD_Skin_GSE153760/AD' + txt_file

COPD_Lung_GSE136831_Control = partial_path + 'COPD_Lung_GSE136831/Control' + txt_file
COPD_Lung_GSE136831_COPD = partial_path + 'COPD_Lung_GSE136831/COPD' + txt_file
COPD_Lung_GSE136831_IPF = partial_path + 'COPD_Lung_GSE136831/IPF' + txt_file

COPD_Lung_GSE171541_control = partial_path + 'COPD_Lung_GSE171541/control' + txt_file
COPD_Lung_GSE171541_copd = partial_path + 'COPD_Lung_GSE171541/copd' + txt_file

Human_Lung_Cell_Atlas_normal = partial_path + 'Human_Lung_Cell_Atlas/normal' + txt_file

ILD_Lung_GSE135893_Control = partial_path + 'ILD_Lung_GSE135893/Control' + txt_file
ILD_Lung_GSE135893_IPF = partial_path + 'ILD_Lung_GSE135893/IPF' + txt_file

IgAN_Cell_Reports_Zheng_Kidney_normal_control = partial_path + 'IgAN_Cell_Reports_Zheng_Kidney/normal_control' + txt_file
IgAN_Cell_Reports_Zheng_Kidney_IgAN = partial_path + 'IgAN_Cell_Reports_Zheng_Kidney/IgAN' + txt_file

PSO_Skin_GSE173706_Control = partial_path + 'PSO_Skin_GSE173706/Control' + txt_file
PSO_Skin_GSE173706_PSO_lesion = partial_path + 'PSO_Skin_GSE173706/PSO_lesion' + txt_file
PSO_Skin_GSE173706_PSO_non_lesion = partial_path + 'PSO_Skin_GSE173706/PSO_non_lesion' + txt_file

PSO_Skin_GSE220116_Control__pre_tx = partial_path + 'PSO_Skin_GSE220116/Control__pre_tx' + txt_file
PSO_Skin_GSE220116_PSO__pre_tx = partial_path + 'PSO_Skin_GSE220116/PSO__pre_tx' + txt_file

UC_Colon_GSE116222_Healthy = partial_path + 'UC_Colon_GSE116222/Healthy' + txt_file
UC_Colon_GSE116222_UC_Inflamed = partial_path + 'UC_Colon_GSE116222/UC_Inflamed' + txt_file
UC_Colon_GSE116222_UC_Non_Inflamed = partial_path + 'UC_Colon_GSE116222/UC_Non_Inflamed' + txt_file

UC_Colon_SCP259_Healthy = partial_path + 'UC_Colon_SCP259/Healthy' + txt_file
UC_Colon_SCP259_UC_Inflamed = partial_path + 'UC_Colon_SCP259/UC_Inflamed' + txt_file
UC_Colon_SCP259_UC_Non_Inflamed = partial_path + 'UC_Colon_SCP259/UC_Non_Inflamed' + txt_file

UC_Colon_GSE231993_Control = partial_path + 'UC_Colon_GSE231993/Control' + txt_file
UC_Colon_GSE231993_UC_inflamed = partial_path + 'UC_Colon_GSE231993/UC_inflamed' + txt_file
UC_Colon_GSE231993_UC_uninflamed = partial_path + 'UC_Colon_GSE231993/UC_uninflamed' + txt_file

Skin_Cell_Atlas_Healthy__non_lesion = partial_path + 'Skin_Cell_Atlas/Healthy__non_lesion' + txt_file
Skin_Cell_Atlas_Eczema__lesion = partial_path + 'Skin_Cell_Atlas/Eczema__lesion' + txt_file
Skin_Cell_Atlas_Eczema__non_lesion = partial_path + 'Skin_Cell_Atlas/Eczema__non_lesion' + txt_file
Skin_Cell_Atlas_Psoriasis__lesion = partial_path + 'Skin_Cell_Atlas/Psoriasis__lesion' + txt_file
Skin_Cell_Atlas_Psoriasis__non_lesion = partial_path + 'Skin_Cell_Atlas/Psoriasis__non_lesion' + txt_file

ILD_Lung_GSE1122960_Control = partial_path + 'ILD_Lung_GSE1122960/Control' + txt_file
ILD_Lung_GSE1122960_IPF = partial_path + 'ILD_Lung_GSE1122960/IPF' + txt_file

SLE_Phase2_scRNA_Control = partial_path + 'SLE_Phase2_scRNA/Control' + txt_file
SLE_Phase2_scRNA_SLE = partial_path + 'SLE_Phase2_scRNA/SLE' + txt_file

dirs = [AD_Skin_GSE147424_Healthy, AD_Skin_GSE147424_Lesional, AD_Skin_GSE147424_Non_Lesional,
AD_Skin_GSE153760_HC, AD_Skin_GSE153760_AD,
COPD_Lung_GSE136831_Control,COPD_Lung_GSE136831_COPD, COPD_Lung_GSE136831_IPF,
COPD_Lung_GSE171541_control, COPD_Lung_GSE171541_copd,
Human_Lung_Cell_Atlas_normal,
ILD_Lung_GSE135893_Control, ILD_Lung_GSE135893_IPF,
IgAN_Cell_Reports_Zheng_Kidney_normal_control, IgAN_Cell_Reports_Zheng_Kidney_IgAN,
PSO_Skin_GSE173706_Control, PSO_Skin_GSE173706_PSO_lesion, PSO_Skin_GSE173706_PSO_non_lesion,
PSO_Skin_GSE220116_Control__pre_tx, PSO_Skin_GSE220116_PSO__pre_tx,
UC_Colon_GSE116222_Healthy, UC_Colon_GSE116222_UC_Inflamed, UC_Colon_GSE116222_UC_Non_Inflamed,
UC_Colon_SCP259_Healthy, UC_Colon_SCP259_UC_Inflamed, UC_Colon_SCP259_UC_Non_Inflamed,
UC_Colon_GSE231993_Control, UC_Colon_GSE231993_UC_inflamed, UC_Colon_GSE231993_UC_uninflamed,
Skin_Cell_Atlas_Healthy__non_lesion, Skin_Cell_Atlas_Eczema__lesion, Skin_Cell_Atlas_Eczema__non_lesion, Skin_Cell_Atlas_Psoriasis__lesion, Skin_Cell_Atlas_Psoriasis__non_lesion,
ILD_Lung_GSE1122960_Control, ILD_Lung_GSE1122960_IPF,
SLE_Phase2_scRNA_Control, SLE_Phase2_scRNA_SLE] 

############################################################
########## DATA FRAME FUNCTION FOR EACH GROUP ##############
############################################################

# for loop
for directory in dirs: 
  
# tab delimiter, first column as index, first row is column name
  sig_means = pd.read_csv(directory, sep='\t', index_col=0, header=0)

# specify column, apply lambda function to each value in column, convert to a string and check if it starts w/ either "CC" or "CX" 
  gene_ligand = sig_means['gene_a'].apply(lambda x: str(x).startswith("CC") or str(x).startswith('CX'))
  gene_receptor = sig_means['gene_b'].apply(lambda x: str(x).startswith("CC") or str(x).startswith('CX'))

# select and retrieve (.index) indices where boolean condition is true in either and create a set 
  keeprs = set(gene_ligand[gene_ligand==True].index) | set(gene_receptor[gene_receptor==True].index)

# replace missing values with 0,  inplace = true means modify object directly  
  sig_means.fillna(0, inplace=True)

# if column name has "|" in it 
  interactome = set([i for i in sig_means.columns if '|' in i])

# sigmeans[interactome] = select columns in this set and make new df
# .apply(sum, 1) = apply this to each row and calc the sum
# assign the sums to the totals variable 
  totals = sig_means[interactome].apply(sum, 1)

# filter totals to keep only rows where sum is >0 
  totals_sig = totals[totals > 0].index.to_list()

# select columns specified by interactome & apply sum to each column 
  totals_inter = sig_means[interactome].apply(sum, 0)
  
# keep only columns where sum is >0, retrieve the columns and convert to a set list
  totals_inter_sig = set(totals_inter[totals_inter > 0].index.to_list())

# select rows specified 
  sig_only = sig_means.loc[totals_sig]

# loc selects rows from sig only based on keeprs and  sig_only 
  sig_only_liga = sig_only.loc[keeprs & set(sig_only.index)]

# retrieve values from interacting_pair column, convert to list 
  sig_only_liga.index = sig_only_liga['interacting_pair'].to_list()

# calculate sum of values down interactome columns, keep if sum >0
  totals_nonzero = sig_only_liga[interactome].apply(sum, 0)
  totals_nonzero_sig = totals_nonzero[totals_nonzero > 0].index.to_list()

# keep columns that don't have |
  cols_keep = [i for i in sig_only_liga.columns if '|' not in i] + totals_nonzero_sig
  sig_only_liga = sig_only_liga[cols_keep]

# keep specified column
  sig_plot = sig_only_liga[totals_nonzero_sig]

# loc selects rows from sig plot based on sig_plot 
  sig_plot = sig_plot.loc[sorted(sig_plot.index)]

# sort the columns
  sig_plot = sig_plot[sorted(sig_plot.columns)]

############################################################
####################### TABLE STUFF ########################
############################################################

# make table
  df = sig_plot
  df = df.rename_axis('ligand_receptor').reset_index()
  melted_df = pd.melt(df, id_vars=['ligand_receptor'], var_name='cell1_cell2', value_name='Value')
  melted_df_filtered = melted_df[melted_df['Value'] !=0]
  melted_df_filtered[['ligand', 'receptor']] = melted_df_filtered['ligand_receptor'].str.split('_', expand=True)
  melted_df_filtered[['cell1', 'cell2']] = melted_df_filtered['cell1_cell2'].str.split('|', expand=True)

# save table 
  table_name = directory + '_table_chemokine.csv'
  melted_df_filtered.to_csv(os.path.join(directory, table_name), index=False)
