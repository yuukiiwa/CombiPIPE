# CombiPIPE
An analysis pipeline for the CombiGEM and CombiSEAL platforms
### To run the pipeline
```
nextflow CombiPIPE.nf --fastq <fastq> --sampinfo <sampleInfo.csv> --barcodes <barcode_list.csv> --pattern <1st 7 letters from the fastq file> --dimensions <number of gRNAs> --linker <barcode-connecting sequence>
```
### Notes
For CombiSEAL libraries, please see [epistasisCalculator](https://github.com/AWHKU/epistasisCalculator).
### CombiGEM publications
1. https://www.nature.com/articles/nbt.3326
2. https://www.pnas.org/content/113/9/2544.short

