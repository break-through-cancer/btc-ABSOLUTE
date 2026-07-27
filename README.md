# btc-ABSOLUTE
Nextflow wrapper for ABSOLUTE v1.5

ABSOLUTE assesses tumor purity and ploidy, and runs in two phases:
- Phase I: generate all possible purity/ploidy solutions as well as Rdata object for force calling mutations at specific locations
- Phase II: after manually selecting a purity/ploidy solution, conduct force calling on the original Phase I input data.

To run:
```
nextflow run main.nf --samplesheet path/to/samplesheet.csv --outdir . --phase {1 or 2}
```

This pipeline is formatted for batch processing of samples, as specified in the samplesheet. Samplesheet file paths must be absolute or relative to the directory from which `main.nf` is run.

Samplesheet Phase I format:
```
sample,seg_path,indel_path,snp_path
sample_1,*.capseg.txt,*.indel,*.snp
sample_2,*.capseg.txt,*.indel,*.snp
```

Samplesheet Phase II format
```
sample,seg_path,indel_path,snp_path,purity,ploidy,rdata_path
sample_1,*.capseg.txt,*.indel,*.snp,0.00,0.00,phase_1_outdir/*.PP-modes.data.RData
sample_2,*.capseg.txt,*.indel,*.snp,0.00,0.00,phase_1_outdir*.PP-modes.data.RData
```

# ABSOLUTE Cirro Use Instructions 
To run ABSOLUTE, you need ICONICC outputs, specifically the processed_counts and seg.txt files. Moreover, in your ICONICC dataset there needs to be a samplesheet to link the files to the sample that is affiliated with it. 

The samplesheet for the ICONICC dataset needs to be structured as follows: 

```
sample,file
DFCI3-S2-L1,GBM1.DFCI3.S2.L1.seg.txt
DFCI3-S2-L1,GBM1.DFCI3.S2.L1_processed_counts.txt
```

The sample name should also match up with the sample name of the data used to generate the .maf file from Funcotator or a similar component. This consistency in sample naming will help downstream components to connect files coordinating to the same sample. 

The inputs to ABSOLUTE Preprocess are a maf file and ICONICC processed_counts and seg.txt files

The inputs to ABSOLUTE part 1 should be the output of ABSOLUTE Preprocess

After running ABSOLUTE part 1, update the sample's purity and ploidy by navigating to Cirro Projects -> Samples -> clicking on the right sample name and updating the purity and ploidy fields
<img width="1793" height="192" alt="Screenshot 2026-07-27 at 3 29 05 PM" src="https://github.com/user-attachments/assets/3cb95c9f-1042-4ba1-b49a-6aaa1ccc8a0a" />


The inputs to ABSOLUTE part 2 should be the output of both ABSOLUTE Preprocess and ABSOLUTE part 1, as well as making sure you have updated the purity/ploidy for the sample

When using the Cirro pipeline, you can choose datasets that contain the input files and, if the setup is correct, the samples that can be analyzed should appear.


