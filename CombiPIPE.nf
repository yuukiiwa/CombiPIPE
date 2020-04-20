def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow CombiPIPE.nf <blahblahblah.fastq> <sampleInfo.csv> <barcodes.txt> <the 1st 7 characters in the fastq file>
    Mandatory arguments:
     --fastq           A fastq file (raw NGS data)
     --sampinfo        A CSV file with sample name, index sequence, and group of each of the samples.
     --barcodes        A TXT file with all barcodes used (a barcode a line)
     --pattern	       The 1st 7 characters in the fastq file
     --dimensions      The number of modules in the laundary
     --linker          Linker sequence b/w the barcodes (AWp12 uses CAATTC)

    Other options:
     --help                         Generates the help page
    """.stripIndent()
}


// Show help message
if (params.help){
    helpMessage()
    exit 0
}
if (params.fastq) { ch_fastq = file(params.fastq, checkIfExists: true) } else { exit 1, "Please provide NGS data!" }
if (params.sampinfo) { ch_sampinfo = file(params.sampinfo, checkIfExists: true) } else { exit 1, "Please provide sample info (csv file)!" }
if (params.barcodes) { ch_barcodes = file(params.barcodes, checkIfExists: true) } else { exit 1, "Please provide all barcodes used!" }
if (params.pattern) { ch_pattern = params.pattern } else { exit 1, "Please provide the first 7 characters of the fastq file!" }
if (params.linker) { ch_linker = params.linker } else { exit 1, "Please provide the linker sequence! (AWp12 uses CAATTC)" }
if (params.dimensions) { ch_dimensions = params.dimensions } else { exit 1, "Please provide the number of dimensions (modules)!" }
ch_process1 = file(params.process1, checkIfExists: true)
ch_process2 = file(params.process2, checkIfExists: true)

process sortSamples {
  publishDir "./outputs_sortSamples", mode: 'copy',
        saveAs: { filename ->
                      if (!filename.endsWith(".version")) filename
                }
  input:
  file sampinfo from ch_sampinfo
  file fastq from ch_fastq
  val  pattern from ch_pattern
  file process1 from ch_process1
  
  output:
  path("*.fastq") into ch_process1_out
   
  script:
  """
  python $process1 $fastq $sampinfo $pattern
  """
}

process newBCanalyzer {
  publishDir "./outputs_newBCanalyzer", mode: 'copy',
        saveAs: { filename ->
                      if (!filename.endsWith(".version")) filename
                }
  input:
  file barcode_list from ch_barcodes
  val linker from ch_linker
  val dimensions from ch_dimensions  
  val pattern from ch_pattern
  file process2 from ch_process2 
  path fastq from ch_process1_out.flatten()

  output:
  file("*.csv") into ch_process2_out

  script:
  """  
  python $process2 -b $barcode_list -f $fastq -l $linker -n $dimensions -p $pattern
  """
}
