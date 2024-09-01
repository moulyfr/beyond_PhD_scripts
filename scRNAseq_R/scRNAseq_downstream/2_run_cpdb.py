#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

def run_cpdb_func(cpdb_file_path, dataset_name, group_name, counts_data_var, partial_path):
    meta_file_path = f"{partial_path}cpdb_inputs/{dataset_name}_cellbridge.rds/{group_name}/meta.txt"
    counts_file_path = f"{partial_path}cpdb_inputs/{dataset_name}_cellbridge.rds/{group_name}"
    output_path = f"{partial_path}cpdb_outputs/{dataset_name}/{group_name}"
    
    cpdb_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path=cpdb_file_path,
        meta_file_path=meta_file_path,
        counts_file_path=counts_file_path,
        counts_data=counts_data,
        output_path=output_path,
        output_suffix=""
    )
    
    print(f"# # # # # # # # # # # # # # # # Hello, {dataset_name} - {group_name} is done!")

# define macro paths
cpdb_file_path = '.../cellphonedb.zip'
counts_data = 'hgnc_symbol'
partial_path = '.../CellphoneDB_analysis/'

# Run func for each dataset groups
#run_cpdb_func(cpdb_file_path, 'AD_Skin_GSE147424', 'Healthy', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'AD_Skin_GSE147424', 'Lesional', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'AD_Skin_GSE147424', 'Non_Lesional', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'AD_Skin_GSE153760', 'HC', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'AD_Skin_GSE153760', 'AD', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'COPD_Lung_GSE136831', 'Control', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'COPD_Lung_GSE136831', 'COPD', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'COPD_Lung_GSE136831', 'IPF', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'COPD_Lung_GSE171541', 'control', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'COPD_Lung_GSE171541', 'copd', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'Human_Lung_Cell_Atlas', 'normal', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'ILD_Lung_GSE135893', 'Control', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'ILD_Lung_GSE135893', 'IPF', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'IgAN_Cell_Reports_Zheng_Kidney', 'normal_control', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'IgAN_Cell_Reports_Zheng_Kidney', 'IgAN', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'PSO_Skin_GSE173706', 'Control', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'PSO_Skin_GSE173706', 'PSO_lesion', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'PSO_Skin_GSE173706', 'PSO_non_lesion', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'PSO_Skin_GSE220116', 'Control__pre_tx', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'PSO_Skin_GSE220116', 'PSO__pre_tx', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'UC_Colon_GSE116222', 'Healthy', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'UC_Colon_GSE116222', 'UC_Inflamed', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'UC_Colon_GSE116222', 'UC_Non_Inflamed', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'UC_Colon_SCP259', 'Healthy', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'UC_Colon_SCP259', 'UC_Inflamed', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'UC_Colon_SCP259', 'UC_Non_Inflamed', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'UC_Colon_GSE231993', 'Control', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'UC_Colon_GSE231993', 'UC_inflamed', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'UC_Colon_GSE231993', 'UC_uninflamed', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'Skin_Cell_Atlas', 'Healthy__non_lesion', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'Skin_Cell_Atlas', 'Eczema__lesion', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'Skin_Cell_Atlas', 'Eczema__non_lesion', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'Skin_Cell_Atlas', 'Psoriasis__lesion', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'Skin_Cell_Atlas', 'Psoriasis__non_lesion', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'ILD_Lung_GSE1122960', 'Control', counts_data, partial_path)
#run_cpdb_func(cpdb_file_path, 'ILD_Lung_GSE1122960', 'IPF', counts_data, partial_path)

#run_cpdb_func(cpdb_file_path, 'SLE_Phase2_scRNA', 'Control', counts_data, partial_path)
run_cpdb_func(cpdb_file_path, 'SLE_Phase2_scRNA', 'SLE', counts_data, partial_path)
