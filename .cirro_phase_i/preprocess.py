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

ds.add_param('sample', ds.samplesheet.iloc[0]['sample'])

param_list = ["sample","seg_path","indel_path","snp_path"]

samplesheet = pd.DataFrame([{k: ds.params.get(k) for k in param_list}])
samplesheet.to_csv('samplesheet.csv', index=False)
ds.add_param("samplesheet", "samplesheet.csv")

keep_params = ['phase', 'samplesheet', 'outdir']

for key in list(ds.params.keys()):  # list() avoids modifying during iteration
    if key not in keep_params:
        ds.remove_param(key, force=True)

ds.logger.info(ds.params)