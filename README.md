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