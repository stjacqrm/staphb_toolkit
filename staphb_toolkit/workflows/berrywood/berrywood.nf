#!/usr/bin/env nextflow

//Description: dealing with the problem children
//Author: Rachael St. Jacques
//email: rachael.stjacques@dgs.virginia.gov

//starting parameters
params.reads = ""
params.outdir = ""
params.primers =""
params.report = ""

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
    .set { raw_reads }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(reads) into read_files_fastqc, read_files_trimming

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
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads
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

//Step3: Assemble cleaned reads with iVar
process ivar {
  tag "$name"

  publishDir "${params.outdir}/logs/alignments", mode: 'copy',pattern:"*.sorted.bam"
  //publishDir "${params.outdir}/SC2_reads", mode: 'copy',pattern:"*_SC2*.fastq.gz"
  publishDir "${params.outdir}/logs/assemblies", mode: 'copy', pattern: "*_consensus.fasta"


  input:
  set val(name), file(reads) from cleaned_reads

  output:
  tuple name, file("${name}_consensus.fasta") into assembled_genomes, annotation_genomes, converting_genomes
  tuple name, file("${name}.sorted.bam") into alignment_file
  tuple name, file("${name}_SC2*.fastq.gz") into sc2_reads

  shell:
"""
ln -s /reference/nCoV-2019.reference.fasta ./nCoV-2019.reference.fasta
minimap2 -K 20M -x sr -a ./nCoV-2019.reference.fasta !{reads[0]} !{reads[1]} | samtools view -u -h -F 4 - | samtools sort > SC2.bam
samtools index SC2.bam
samtools flagstat SC2.bam
samtools sort -n SC2.bam > SC2_sorted.bam
samtools fastq -f2 -F4 -1 ${name}_SC2_R1.fastq.gz -2 ${name}_SC2_R2.fastq.gz SC2_sorted.bam -s singletons.fastq.gz

ivar trim -i SC2.bam -b /reference/ARTIC-${params.primers}.bed -p ivar -e

samtools sort  ivar.bam > ${name}.sorted.bam
samtools index ${name}.sorted.bam
samtools flagstat ${name}.sorted.bam

samtools mpileup -f ./nCoV-2019.reference.fasta -d 1000000 -A -B -Q 0 ${name}.sorted.bam | ivar consensus -p ivar -m ${params.ivar_mindepth} -t ${params.ivar_minfreq} -n N
echo '>${name}' > ${name}_consensus.fasta

seqtk seq -U -l 50 ivar.fa | tail -n +2 >> ${name}_consensus.fasta
"""
}

//QC of read data
process samtools {
  tag "$name"

  input:
  set val(name), file(alignment) from alignment_file

  output:
  file "${name}_samtoolscoverage.tsv" into alignment_qc

  shell:
  """
  samtools coverage ${alignment} -o ${name}_samtoolscoverage.tsv
  """
}

//Collect and format report
process assembly_results{
  publishDir "${params.outdir}/logs/assemblies/", mode: 'copy', pattern: "*assembly_metrics.csv"

  echo true

  input:
  file(cg_pipeline_results) from alignment_qc.collect()
  file(assemblies) from assembled_genomes.collect()

  output:
  file "*assembly_metrics.csv"


  script:
"""
#!/usr/bin/env python3
import os, sys
import glob, csv
import xml.etree.ElementTree as ET
from datetime import datetime

today = datetime.today()
today = today.strftime("%m%d%y")

class result_values:
    def __init__(self,id):
        self.id = id
        self.aligned_bases = "NA"
        self.percent_cvg = "NA"
        self.mean_depth = "NA"
        self.mean_base_q = "NA"
        self.mean_map_q = "NA"
        self.status = "NA"


#get list of result files
samtools_results = glob.glob("*_samtoolscoverage.tsv")
assemblies = glob.glob("*consensus.fasta")

results = {}

# collect samtools results
for file in samtools_results:
    id = file.split("_samtoolscoverage.tsv")[0]
    result = result_values(id)
    status = []
    with open(file,'r') as tsv_file:
        tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))
        for line in tsv_reader:
            result.aligned_bases = line["covbases"]
            result.percent_cvg = line["coverage"]
            if float(line["coverage"]) < 98:
                status.append("coverage <98%")
            result.mean_depth = line["meandepth"]
            result.mean_base_q = line["meanbaseq"]
            if float(line["meanbaseq"]) < 30:
                status.append("meanbaseq < 30")
            result.mean_map_q = line["meanmapq"]
            if float(line["meanmapq"]) < 30:
                status.append("meanmapq < 30")
        if len(status) == 0:
            result.status = "PASS"
        else:
            result.status ="WARNING: " + '; '.join(status)

    results[id] = result



#create output file
with open(f"{today}_assembly_metrics.csv",'w') as csvout:
    writer = csv.writer(csvout,delimiter=',')
    writer.writerow(["sample","aligned_bases","percent_cvg", "mean_depth", "mean_base_q", "mean_map_q", "status"])
    for id in results:
        result = results[id]
        writer.writerow([result.id,result.aligned_bases,result.percent_cvg,result.mean_depth,result.mean_base_q,result.mean_map_q,result.status])
"""

}

//Step who knows: annotate the SC2 assembled_genomes
process annotate {
  publishDir "${params.outdir}/results/annotations/", mode: 'copy'
//  publishDir "${params.outdir}/convert_gff/", mode: 'copy'

  input:
  set val(name), file(assembly) from annotation_genomes

  output:
  set val(name), file("${name}_clean.gff3") into annotations, convert_gff

  script:
  """
  vigor4 -i ${assembly} -o ${name} -d sarscov2 -f 2
  grep -hr "gene" ${name}.gff3 > ${name}_clean.gff3 && grep -hr "exon" ${name}.gff3 >> ${name}_clean.gff3
  """

}

process convert_gff_gb {
  //errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/results/genbank_files", mode: 'copy' , pattern:"*.gb"

  input:
  set val(name), file(assemblies) from converting_genomes
  set val(name), file(gff_files) from convert_gff

  output:
  set val(name), file('*.gb') into genbank_files, genome_viz, genome_name
  //file("${name}.gbk") into genome_viz

  script:
  """
  gff_to_genbank.py ${gff_files} ${assemblies}
  for file in *.gb; do mv "\$file" "\${file/_*.gb/.gb}"; done
  """

}

process viz_genomes {
  //errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/results/genome_maps/sample_genomes", mode: 'copy', pattern:"*_genome_map.png"
  publishDir "${params.outdir}/results/genome_maps/mutations_genome_maps", mode: 'copy', pattern:"*mutations.png"
  publishDir "${params.outdir}/results/genome_maps/NC_045512_2_genome_map", mode: 'copy', pattern:"NC_045512_2_genome_map.png"
  publishDir "${params.outdir}/results/mutations/csv_files", mode: 'copy', pattern:"*.csv"
  publishDir "${params.outdir}/results/mutations/text_files/", mode: 'copy', pattern:"*.txt"

  input:
  set val(name), file(genbank) from genome_viz

  output:
  file("${name}_mutations.png") into mutations_viz
  file("NC_045512_2_genome_map.png") into master_map
  file("${name}_genome_map.png") into sample_genome
//  file("${name}.csv")
  file("${name}_diff_blocks_as_features_dict.txt")
  file("${name}_2.csv")
  file("${name}_3.csv")
//  file("${name}_diff_blocks_plain_dict.txt")


  script:
  """
  #!/usr/bin/env python
  import sys
  import argparse
  from geneblocks import CommonBlocks, load_record, DiffBlocks
  from dna_features_viewer import BiopythonTranslator
  from dna_features_viewer import GraphicRecord
  from dna_features_viewer import CircularGraphicRecord
  from dna_features_viewer import GraphicFeature
  import matplotlib.pyplot as plt
  import Bio
  import pandas as pd
  import csv
  from pylab import title, figure, xlabel, ylabel, xticks, bar, legend, axis, savefig

  # make the custom translator
  class CustomTranslator(BiopythonTranslator):

      label_fields = [
          "Name",
          "note",
          "gene",
          "product",
          "source",
          "locus_tag",
          "label",
      ]

      ignored_features_types = ("diff_equal","CDS","mRNA","source")
      default_box_color = None


      def compute_feature_color(self, feature):
          if feature.qualifiers.get("is_block", False):
              return BiopythonTranslator.compute_feature_color(self, feature)
          elif feature.type == "gene":
              return "green"
          elif feature.type == "exon":
              return "blue"
          elif feature.type == "note":
              return "pink"
          elif feature.type == "mRNA_with_minus_1_frameshift":
              return "gold"
          else:
              return "grey"

      @staticmethod
      def compute_feature_box_linewidth(f):
          return 1 if f.qualifiers.get("is_block", False) else 0

      @staticmethod
      def compute_feature_fontdict(f):
          return {"fontsize": 12 if f.qualifiers.get("is_block", False) else 9}

      def compute_feature_label(self, feature):
          if feature.type == 'CDS':
              return None
          elif feature.type == "mRNA":
              return None
          elif feature.type == "source":
              return None
          elif feature.type == "exon":
              return "exon"
          elif feature.type == "mRNA_with_minus_1_frameshift":
              return "frameshift"
          else:
              return BiopythonTranslator.compute_feature_label(self, feature)

  cr = CustomTranslator()

  seq_1 = cr.translate_record("/reference/NC_045512_2.gb")
  seq_1.ticks_resolution = 3000
  ax, _ = seq_1.plot(figure_width=30)
  plt.suptitle("Genome map of NC_045512_2")
  ax.figure.tight_layout()
  ax.figure.savefig("NC_045512_2_genome_map.png")

  seq_2 = cr.translate_record("${genbank}")
  seq_2.ticks_resolution = 3000
  ax, _ = seq_2.plot(figure_width=30)
  plt.suptitle("Genome map of ${name}")
  ax.figure.tight_layout()
  ax.figure.savefig("${name}_genome_map.png")

  seq_11 = load_record("/reference/NC_045512_2.gb")
  seq_22 = load_record("${genbank}")

  diff_blocks = DiffBlocks.from_sequences(seq_11, seq_22)
  dict = diff_blocks.diffs_as_features()
  f = open("${name}_diff_blocks_as_features_dict.txt", "a")
  f.write(str(dict))
  f.close()

  df = pd.DataFrame(dict)
  df.to_csv("${name}_1.csv")

  df.to_csv("${name}_2.csv")

  diff_features_df = pd.DataFrame(dict)
  diff_features_df.to_csv("${name}_3.csv")

  ax1, ax2 = diff_blocks.plot(figure_width=30)
  plt.suptitle("Changes in ${name} when compared to NC_045512_2.gb")
  ax1.figure.savefig("${name}_mutations.png")
  """
}

process build_report {
  //errorStrategy 'ignore'
  tag "$name"
//  publishDir "${params.outdir}/reports/complete_reports", mode: 'copy', pattern:"*report.pdf"
  publishDir "${params.outdir}/results/reports", mode: 'copy', pattern:"*report.md"
  publishDir "${params.outdir}/results/reports", mode: 'copy', pattern:"*report.pdf"

  input:
//  set val(genbank), file(genbank) from genome_viz
  set val(name), file(genbank) from genome_name
  file(mutation_map) from mutations_viz
  file(reference_gene_map) from master_map
  file(sample_genomes) from sample_genome
//  file(rmd) from report

  output:
  file("*report.md") //into mutations_viz
  file("*report.pdf")

  script:
  """
  date=\$(date '+%m%d%y')
  touch ${name}_report.md
  cat >> ${name}_report.md << EOF
  # ${name} Mutation Analysis

  <hr style="border:2px solid gray"> </hr>


  ### A comparison of a SARS-CoV-2 isolate to the Wuhan reference genome to determine areas of genomic abnormalities
  Genome maps of the NC_045512.2 reference genome and Query genome are generated to show where the predicetd genes would be encoded in the respective genomes. The reference genome map, NC_045512.2, will always be drawn in the same way.
  If there are mutations in the query genome, they will be relfected in the Query Genome Map.

  Each "query" genome is compared to the NC_045512.2 (Wuhan) reference genome to identify where mutations are in the genome. These are vizualized in Figure 3. Insertions are shown in green, deletions are shown in red,
  and "swaps" are shown in yellow. Mutations are shown in the figure as amino acids.

  This report is intended only for research purposes. No bioinformatics method is designed to replace a thourough, in-depth epidemiological investigation.


  ### NC_045512.2 Genome Map

  ![NC_045512.2 Reference Genome Map](${reference_gene_map})

  *This report is intended only for research purposes. No bioinformatics method is designed to replace a thourough, in-depth epidemiological investigation.*

  <div style="page-break-after: always; visibility: hidden">
  </div>

  ### Query Genome Map

  ![Query Genome Map](${sample_genomes})

  *This report is intended only for research purposes. No bioinformatics method is designed to replace a thourough, in-depth epidemiological investigation.*

  <div style="page-break-after: always; visibility: hidden">
  </div>

  ### Mutations in the Query genome compared to the NC_045512_2 genome

  ![Mutations in Query vs. Reference](${mutation_map})

  *This report is intended only for research purposes. No bioinformatics method is designed to replace a thourough, in-depth epidemiological investigation.*

  <div style="page-break-after: always; visibility: hidden">
  </div>


  ## Methods

  - The figures shown here were generated using sequence data processed with the [Berrywood](https://github.com/StaPH-B/staphb_toolkit/tree/dev) pipeline.

  - Berrywood uses [Vigor4](https://github.com/JCVenterInstitute/VIGOR4) to annotate the genomes of the SARS-CoV-2 samples provided.

  - Data from Berrywood is curated into a single pdf report using the [StaPH-B Berrywood Cluster Report Environment](https://hub.docker.com/r/staphb/cluster-report-env).

  - Questions and/or concerns can be directed to [Logan Fink](logan.fink@dgs.virginia.gov) and/or [Rachael St. Jacques](rachael.stjacques@dgs.virginia.gov).

  *This report is intended only for research purposes. No bioinformatics method is designed to replace a thourough, in-depth epidemiological investigation.*

  EOF
  pandoc ${name}_report.md -V geometry:"top=2cm, bottom=1.5cm, left=2cm, right=2cm" -o ${name}_report.pdf
  """
}
