# CombiPIPE
An analysis pipeline for the CombiGEM (done) and CombiSEAL (will add this next week) platforms
## Running the pipelines
The first three processes are the same in both pipelines.
1. sample extraction
2. barcode extraction
3. calculate lg Fold Change and -log10 P value between the initial and experimental group.
#### For CombiGEM (supporting pairwise screen only in this pipeline)
4. calculate gene-level genetic interaction
5. calculate Dunnett test p-values for each gene-level combinations
```
nextflow CombiGEM.nf --fastq <fastq> --sampinfo <sampleInfo.csv> --barcodes <barcode_list.csv> --pattern <1st 7 letters from the fastq file> --dimensions <number of gRNAs> --linker <barcode-connecting sequence>
```
#### For CombiSEAL (coming soon)
4. need to talk to people who were in that project.
5. epistasis calculation
```
nextflow CombiSEAL.nf --fastq <fastq> --sampinfo <sampleInfo.csv> --barcodes <barcode_list.csv> --pattern <1st 7 letters from the fastq file> --dimensions <number of modules> --linker <barcode-connecting sequence>
```
## CombiGEM publications
1. https://www.nature.com/articles/nbt.3326
2. https://www.pnas.org/content/113/9/2544.short
