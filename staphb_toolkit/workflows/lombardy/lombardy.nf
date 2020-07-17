#!/usr/bin/env nextflow

//Description: Workflow for the hybrid assembly of Nanpore and Illumina reads, subsequent prophage identification/ annotation and plasmid/plasmid AMR predeiction
//Author: Rachael St. Jacques
//eMail: rachael.stjacques@dgs.virginia.gov

//starting parameters
params.reads = ""
params.metadata = ""
params.nanopore_reads =  ""
params.outdir = ""

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs(  "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { raw_reads }

Channel
    .fromPath(params.metadata)
    .set { metadata }

Channel
    .fromPath(params.nanopore_reads)
    .set { nanopore_reads }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(reads) into raw_reads_assemble, raw_reads_snp, raw_reads_qc

  script:
  if(params.name_split_on!=""){
    name = name.split(params.name_split_on)[0]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    mv ${reads[1]} ${name}_R2.fastq.gz
    """
  }else{
  """
  """
  }
}

//Step1a: FastQC
process fastqc {
  tag "$name"
  publishDir "${params.outdir}/logs/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from read_files_fastqc

  output:
  file("*_fastqc.{zip,html}") into fastqc_results

  script:
  """
  fastqc -q  ${reads}
  """
}

//Step1b: Trim with Trimmomatic
process trim {
  tag "$name"
  if(params.savetrimmedreads){
    publishDir "${params.outdir}/trimmed", mode: 'copy'
  }
  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}_trimmed{_1,_2}.fastq.gz") into trimmed_reads
  file("${name}.trim.stats.txt") into trimmomatic_stats

  script:
  """
  cpus=`grep -c ^processor /proc/cpuinfo`
  java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads \$cpus ${reads} -baseout ${name}.fastq.gz SLIDINGWINDOW:${params.windowsize}:${params.qualitytrimscore} MINLEN:${params.minlength} 2> ${name}.trim.stats.txt
  mv ${name}*1P.fastq.gz ${name}_trimmed_1.fastq.gz
  mv ${name}*2P.fastq.gz ${name}_trimmed_2.fastq.gz
  """
}

//Step2: Remove PhiX contamination
process cleanreads {
  tag "$name"
  publishDir "${params.outdir}/logs/cleanedreads", mode: 'copy',pattern:"*.stats.txt"

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads_cg
  file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads_snp
  file("${name}.phix.stats.txt") into phix_cleanning_stats
  file("${name}.adapters.stats.txt") into adapter_cleanning_stats

  script:
  """
  ram=`awk '/MemTotal/ { printf "%.0f \\n", \$2/1024/1024 - 1 }' /proc/meminfo`
  ram=`echo \$ram | awk '{\$1=\$1;print}'`
  repair.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh -Xmx\${ram}g in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh -Xmx\${ram}g in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_1.clean.fastq.gz out2=${name}_2.clean.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  """
}

//Step whatever who cares at this point: Choppin nanopore_reads
process porechop {
  tag "$name"
  publishDir "${params.outdir}/porechopped", mode: 'copy'

  input:
  set val(name), file(reads) from nanopore_reads

  output:
  tuple
}
