#!/usr/bin/env nextflow

//Description: Workflow for identifying and annotating prophage genomes in bactieral isolates
//Author: Rachael St. Jacques
//email: rachael.stjacques@dgs.virginia.gov

//starting parameters
params.reads = ""
params.outdir = ""

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.fastq.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
    .set { raw_reads }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file("*{R1,R2,_1,_2}.fastq.gz") into read_files_fastqc, read_files_trimming

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
  tuple name, file("${name}*{_1,_2}.fastq.gz") into trimmed_reads
  file("${name}.trim.stats.txt") into trimmomatic_stats

  script:
  """
  cpus=`grep -c ^processor /proc/cpuinfo`
  java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads \$cpus ${reads} -baseout ${name}.fastq.gz SLIDINGWINDOW:${params.windowsize}:${params.qualitytrimscore} MINLEN:${params.minlength} 2> ${name}.trim.stats.txt
  mv ${name}*1P.fastq.gz ${name}_1.fastq.gz
  mv ${name}*2P.fastq.gz ${name}_2.fastq.gz
  """
}
//Step2: Remove PhiX contamination
process cleanreads {
  tag "$name"
  publishDir "${params.outdir}/logs/cleanedreads", mode: 'copy',pattern:"*.stats.txt"

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads
  file("${name}.phix.stats.txt") into phix_cleanning_stats
  file("${name}.adapters.stats.txt") into adapter_cleanning_stats

  script:
  """
  repair.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_1.clean.fastq.gz out2=${name}_2.clean.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  """
}

//Step3: Assemble trimmed reads with Shovill
process shovill {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/logs/assemblies", mode: 'copy'

  input:
  set val(name), file(reads) from cleaned_reads

  output:
  tuple name, file("${name}.fasta") into assembled_genomes_quality, assembled_genomes

  shell:
  '''
  ram=`awk '/MemAvailable/ { printf "%.0f \\n", $2/1024/1024 }' /proc/meminfo`
  shovill --cpus 0 --ram $ram  --outdir . --R1 !{reads[0]} --R2 !{reads[1]} --force
  mv contigs.fa !{name}.fasta
  '''
}

//Step3a: Assembly Quality Report
process quast {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/logs/quast",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.quast.tsv") into quast_report

  script:
  """
  quast.py ${assembly} -o .
  mv report.txt ${name}.quast.tsv
  """
}

//Step4: Identify and annotate prophage genomes from assembled raw_reads
process prophage {
  publishDir "${params.outdir}/", mode: 'copy'
  tag "$name"

  input:
  set val(name), file(assembly) from assembled_genomes

  output:
  file("${name}") into prophage_out

  script:
  """
  mkdir assemblies
  mv *.fasta ./assemblies
  VIBRANT_run.py -i ./assemblies/*.fasta
  mv VIBRANT_${name} ${name}
  """
}
