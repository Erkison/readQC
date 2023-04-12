# readQC

## Local run
```
nextflow run ~/Ewomazino/Erkison/code/bioinformatics_scripts/massive/assemblyQC/main.nf \
    --reads "path/to/reads/*_{1,2}.fastq.gz" --output_dir path/to/output/dir \
    -profile local \
    -resume
```

## MASSIVE run
```
nextflow run ~/js66_scratch/erkison/seroepi/scripts/nextflow/assemblyQC/main.nf \
    --reads "path/to/reads/*_{1,2}.fastq.gz" --output_dir path/to/output/dir \
    -profile standard -resume
```
