# readQC

## Local run
```
nextflow run readQC/main.nf \
    --reads "path/to/reads/*{R,_}{1,2}.fastq.gz" --output_dir path/to/output/dir \
    --confindr_db_path path/to/confindr/db \
    -profile local \
    -resume
```

## MASSIVE run
```
nextflow run ~/js66_scratch/erkison/seroepi/scripts/nextflow/readQC/main.nf \
    --reads "path/to/reads/*{R,_}{1,2}.fastq.gz" --output_dir path/to/output/dir \
    --confindr_db_path path/to/confindr/db \
    -profile standard -resume
```
