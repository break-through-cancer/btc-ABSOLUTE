#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd

# 1. Get parameters from cirro pipeline call
ds = PreprocessDataset.from_running()
ds.logger.info("List of starting params")
ds.logger.info(ds.params)

ds.logger.info('checking ds.files')
files = ds.files
ds.logger.info(files.head())
ds.logger.info(files.columns)

# 2. Add samplesheet parameter and set equal to ds.samplesheet
ds.logger.info("Checking samplesheet parameter")
ds.logger.info(ds.samplesheet)

suffix_map = {
    '.capseg.txt': 'seg_path',
    '.indel': 'indel_path',
    '.snp': 'snp_path',
    '.PP-modes.data.RData': 'rdata_path'
}

def map_path(filename):
    for suffix, path_type in suffix_map.items():
        if filename.endswith(suffix):
            return path_type
    return None

files['path_type'] = files['file'].apply(map_path)

result = (
    files.pivot(index='sample', columns='path_type', values='file')
         .rename_axis(columns=None)
         .reset_index()
)
samplesheet = pd.merge(ds.samplesheet, result, on='sample')

param_list = ["sample", "seg_path", "indel_path", "snp_path"]
samplesheet = samplesheet[param_list]

ds.logger.info("Print resulting samplesheet:")
ds.logger.info(samplesheet)

samplesheet.to_csv('samplesheet.csv', index=False)
ds.add_param("samplesheet", "samplesheet.csv")

ds.logger.info(ds.params)