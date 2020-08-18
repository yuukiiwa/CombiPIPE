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
if (params.dummysgs) { ch_dummysgs = params.dummysgs }
ch_process1 = file(params.process1, checkIfExists: true)
ch_process2 = file(params.process2, checkIfExists: true)
ch_process3 = file(params.process3, checkIfExists: true)
ch_process4 = file(params.process4, checkIfExists: true)
ch_process5 = file(params.process5, checkIfExists: true)
ch_process6 = file(params.process6, checkIfExists: true)

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
  path "*.fastq"
  val "outputs_sortSamples" into ch_process1_out
   
  script:
  """
  python3 $process1 $fastq $sampinfo $pattern
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
  val fastq_dir from ch_process1_out

  output:
  file "*.csv"
  val "outputs_newBCanalyzer" into ch_process2_out

  script:
  """  
  python3 $process2 -b $barcode_list -f $PWD/$fastq_dir -l $linker -n $dimensions -p $pattern
  """
}

process log2FC_neglogPV {
  publishDir "./outputs_FCPV", mode: 'copy',
        saveAs: { filename ->
                      if (!filename.endsWith(".version")) filename
                }
  input:
  val BC_dir from ch_process2_out
  file sampinfo from ch_sampinfo
  val dimensions from ch_dimensions
  file process3 from ch_process3

  output:
  file "*.csv" into ch_process3_out

  script:
  """
  python3 $process3 $PWD/$BC_dir $sampinfo $dimensions
  """
}

process PairwiseGI_genDunnettInputs {
  publishDir "./outputs_pairwiseGI", mode: 'copy',
        saveAs: { filename ->
                      if (!filename.endsWith(".version")) filename
                }
  input:
  file fcpv from ch_process3_out
  file process4 from ch_process4
  val dimensions from ch_dimensions
  val dummysgs from ch_dummysgs

  output:
  file "genelevel*" into ch_genelvGI
  file "*.csv"
  val "outputs_pairwiseGI" into ch_process4_out 

  script:
  """
  python3 $process4 $fcpv $dimensions $dummysgs
  """
}

process batchDunnettTest {
  publishDir "./outputs_DunnettTest", mode: 'copy',
        saveAs: { filename ->
                      if (!filename.endsWith(".version")) filename
                }
  input:
  val dir from ch_process4_out
  file process5 from ch_process5

  output:
  file "*.csv"
  val "outputs_DunnettTest" into ch_process5_out

  script:
  """
  Rscript --vanilla $process5 $PWD/$dir
  """
}

process mergeDunnetts {
  publishDir "./output_mergedDunnettTests", mode: 'copy',
        saveAs: { filename ->
                      if (!filename.endsWith(".version")) filename
                }
  input:
  file gi from ch_genelvGI
  val dir from ch_process5_out
  file process6 from ch_process6

  output:
  file "*.csv"

  script:
  """
  python3 $process6 $PWD/$dir $gi
  """
}
