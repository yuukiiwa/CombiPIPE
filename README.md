# CombiPIPE
Analysis pipelines for combinatorial screens done using the CombiGEM-CRISPR and CombiSEAL (coming soon) platforms.
## Before running
Please install DescTools on RStudio. Then go to terminal (any command line window) and install regex, scipy, and numpy for python. If python3, simply use pip3 instead of pip.
```
pip install regex
pip install scipy
pip install numpy
```
Please install nextflow
```
curl -s https://get.nextflow.io | bash
export PATH=$PATH:</path/to/file>
```
Download JDK 8 from [here](https://www.oracle.com/java/technologies/javase-jdk8-downloads.html)
A dmg or exe file serves well. Open once downloaded, then follow the instructions on the installer.
```
export JAVA_HOME=$(/usr/libexec/java_home -v 1.8)
```
Generate a sample info. csv file containing sample index sequences, sample names, and conditions:
```
GATCAATGTTC,SA162,1
CGATCTGGCGAA,SA163,1
TCGTTCCTG,SA164,0
ATCAGAACAT,SA165,0
```
Generate a barcode csv file containing the barcodes and keys:
```
AAGCGAGT,1
CTCTAGGT,2
```
## Running the pipelines
The first three processes are the same in both pipelines.
1. extract samples
2. extract barcodes
3. calculate lg Fold Change and -log10 P value between the initial and experimental group.
#### For CombiGEM (This pipeline supports pairwise screens only)
4. calculate gene-level genetic interaction
5. calculate Dunnett test p-values for each gene-level combinations
```
nextflow CombiGEM.nf 
      --fastq <fastq> \
      --sampinfo <sampleInfo.csv> \
      --barcodes <barcode_list.csv> \
      --pattern <1st 7 letters from the fastq file> \ 
      --dimensions <number of gRNAs> \ 
      --linker <barcode-connecting sequence>
```
#### For CombiSEAL (coming soon)
4. (trying to figure this out)
5. epistasis calculation
```
nextflow CombiSEAL.nf      
      --fastq <fastq> \
      --sampinfo <sampleInfo.csv> \
      --barcodes <barcode_list.csv> \
      --pattern <1st 7 letters from the fastq file> \ 
      --dimensions <number of modules> \ 
      --linker <barcode-connecting sequence>
```
## Publications
#### CombiGEM-CRISPR
1. https://www.nature.com/articles/nbt.3326
2. https://www.pnas.org/content/113/9/2544.short
#### CombiSEAL
1. https://www.nature.com/articles/s41592-019-0473-0
