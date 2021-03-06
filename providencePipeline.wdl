version 1.0

workflow providencePipeline {
input {
    File fastq1
    File fastq2
    File ref
    String samplePrefix
    String sampleLibrary
    String sampleFlowcell
    String sampleExt
    String sampleConstruct
  }

  parameter_meta {
    fastq1: "Read 1 fastq file, gzipped. Can be either targeted or whole transcriptome"
    fastq2: "Read 2 fastq file, gzipped. Can be either targeted or whole transcriptome"
    samplePrefix: "Prefix for output files"
    sampleLibrary: "Library name of the transcriptome"
    sampleFlowcell: "Flowcell name"
    sampleExt: "External name provided for this library"
    ref: "Reference to be used for alignment"
    sampleConstruct: "Library name of the construct used"
  }

  meta {
    author: "Yogi Sundaravadanam"
    email: "ysundaravadanam@oicr.on.ca"
    description: "A pipeline that will confirm the sequence of both the original construct as well as the transcribed product, identifying any variation from the expected sequence"
    dependencies: [
      {
        name: "bbmap/38.75",
        url: "https://sourceforge.net/projects/bbmap/"
      },
      {
        name: "bwa/0.7.12",
        url: "https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2"
      },
      {
        name: "samtools/1.9",
        url: "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2"
      },
      {
        name: "bcftools/1.9",
        url: "https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2"
      },
      {
        name: "vcftools/0.1.16",
        url: "https://vcftools.github.io/downloads.html"
      },
      {
        name: "seqtk/1.3",
        url: "https://github.com/lh3/seqtk/archive/v1.3.tar.gz"
      },
      {
        name: "blast/2.8.1",
        url: "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz"
      },
      {
        name: "rmarkdown/2.3",
        url: "https://cran.r-project.org/web/packages/rmarkdown/index.html"
      },
      {
        name: "emboss/6.6.0",
        url: "ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz"
      }
    ]
    
    output_meta: {
      bam: "Alignment file in BAM format",
      bai: "Alignment file index",
      consensusFasta: "Consenus of the sample in FASTA format",
      variantOnlyVcf: "Mutations in VCF format",
      vcf: "Information on all bases in VCF format",
      needleReport: "global alignment file",
      bl2seqReport: "local alignment file",
      alignmentStats: "Stats generated for the BAM file",
      readDistStats: "Read length distribution metrics",
			insertSizeStats: "Insert size metrics",
      bbMaplog: "Adapter stats",
      report: "PDF report generated by Rmarkdown",
      reportScript: "Rmarkdown script for manual generation of the report",
      orf: "ORF generated by EMBOSS tools"
      } 
  }


  call bbMap {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
      sample = samplePrefix
  }
  
  call bwa {
      input:
      fastq1 = bbMap.out1,  
      fastq2 = bbMap.out2,
      sample = samplePrefix,
			reference = ref
  }

  call variantCalling {
    input:
      sample = samplePrefix,
      bam = bwa.bwaMappedBam,
			reference = ref
  }
  
   call consensus2ReferenceSequence {
    input:
      consensusFasta = variantCalling.consensusFasta,
      sample = samplePrefix,
      reference = ref
  }

  call qcStats {
    input:
      sample = samplePrefix,
      bam = bwa.bwaMappedBam
  }

  call orfStats {
     input:
      consensusFasta = variantCalling.consensusFasta,
      sample = samplePrefix
  }

  call runReport {
    input:
      sample = samplePrefix,
      ext = sampleExt,
      library = sampleLibrary,
      flowcell = sampleFlowcell,
      reference = ref,
      construct = sampleConstruct,
      needleReport = consensus2ReferenceSequence.needleReport,
      consensusFasta = variantCalling.consensusFasta,
      vcf = variantCalling.vcfFile,
      mvcf = variantCalling.variantOnlyVcf,
      bbmapLog = bbMap.bbMapLog,
      samStats = qcStats.alignmentStats,
      readDistStats = qcStats.readDistStats,
      insertSizeStats = qcStats.insertSizeStats,
      orf = orfStats.orf
  }

  output {
    File bam = bwa.bwaMappedBam
    File bai = bwa.bwaMappedBai
    File vcf = variantCalling.vcfFile
    File consensusFasta = variantCalling.consensusFasta
    File variantOnlyVcf = variantCalling.variantOnlyVcf
    File bl2seqReport = consensus2ReferenceSequence.bl2seqReport
    File needleReport = consensus2ReferenceSequence.needleReport
    File alignmentStats = qcStats.alignmentStats
    File readDistStats = qcStats.readDistStats
    File insertSizeStats = qcStats.insertSizeStats
    File bbMaplog = bbMap.bbMapLog
    File report = runReport.rmarkdownReport
    File reportScript = runReport.scriptRmarkdown
    File orf = orfStats.orf
  }
}

task bbMap {
  input {
    String modules = "bbmap/38.75"
    File fastq1
    File fastq2
    String sample
    String reference = "$BBMAP_ROOT/share/bbmap/resources/adapters.fa"
    Int trimq = 25
    Int mem = 8
    Int timeout = 72
  }

  parameter_meta {
    fastq1: "Read 1 fastq file, gzipped"
    fastq2: "Read 2 fastq file, gzipped" 
    sample: "Sample name"
    trimq: "Quality trimming"
    mem: "Memory allocated to trimming task"
    timeout: "Timeout in hours for this task"
    modules: "Names and versions of modules needed for bbmap"
    reference: "Reference file for adaptoer"
  }

  command <<<
    set -euo pipefail

    #Remove adapters and quality trim

    bbmap bbduk in1=~{fastq1} in2=~{fastq2} \
    out1=~{sample}_qad_r1.fastq.gz out2=~{sample}_qad_r2.fastq.gz \
    ref=~{reference} \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=~{trimq} 2>~{sample}_bbmap.log
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File out1 = "~{sample}_qad_r1.fastq.gz"
    File out2 = "~{sample}_qad_r2.fastq.gz"
    File bbMapLog = "~{sample}_bbmap.log"
  }
}

task bwa {
  input {
    String modules = "bwa/0.7.12 samtools/1.9"
    File fastq1
    File fastq2
    String sample
		File reference
    Int mem = 12
    Int timeout = 72
    Int threads = 8
  }
  
  parameter_meta {
    fastq1: "Adapter trimmed Read 1 fastq file, gzipped"
    fastq2: "Adapter trimmed Read 2 fastq file, gzipped"
    sample: "Sample name"
    mem: "Memory allocated for alignment task"
    timeout: "Timeout in hours for this task"
    modules: "Names and versions of modules needed for bwa-mem"
    reference: "Reference file for alignment"
    threads: "Threads required for alignment"
  }

  String bwaMappedBam_ = "~{sample}.bwa.bam"
  String bwaMappedBai_ = "~{sample}.bwa.bam.bai"

  command <<<
    set -euo pipefail

    #index the reference
    ln -s ~{reference} .
    bwa index ~{basename(reference)}
    
    #Align fastq files
    bwa mem -M -t 8 ~{basename(reference)} \
    ~{fastq1}  ~{fastq2} | \
    samtools view -Sb | \
    samtools sort - -o ~{bwaMappedBam_}

    samtools index ~{bwaMappedBam_}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File bwaMappedBam = "~{bwaMappedBam_}"
    File bwaMappedBai = "~{bwaMappedBai_}"
  }
}

task variantCalling {
  input {
    String modules = "bcftools/1.9 samtools/1.9 vcftools/0.1.16 seqtk/1.3"
    File bam
    String sample
    File reference
    Int mem = 8
    Int timeout = 72
  }

  parameter_meta {
    bam: "Alignment file in BAM format"
    sample: "Sample name"
    mem: "Memory allocated for alignment task"
    timeout: "Timeout in hours for this task"
    modules: "Names and versions of modules needed for variant calling"
    reference: "Reference file for variant calling"
  }

  String vcfName = "~{sample}.vcf"
  String fastaName = "~{sample}.consensus.fasta"
  String variantOnlyVcf_ = "~{sample}.v.vcf"

  command <<<
  set -euo pipefail

  #Call consensus sequence
  samtools mpileup -aa -uf ~{reference} ~{bam} | \
     bcftools call --ploidy 1 -Mc | tee -a ~{vcfName} | \
    vcfutils.pl vcf2fq -d 10 | \
    seqtk seq -A - | sed '2~2s/[actg]/N/g' > ~{fastaName}

  bcftools mpileup -a "INFO/AD,FORMAT/DP,FORMAT/AD" \
    -d 8000 -f ~{reference} ~{bam} | \
    tee ~{sample}.m.vcf | bcftools call --ploidy 1 -m -v > ~{variantOnlyVcf_}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"

  }

  output {
    File vcfFile = "~{vcfName}"
    File consensusFasta = "~{fastaName}"
    File variantOnlyVcf = "~{variantOnlyVcf_}"
  }
}

task consensus2ReferenceSequence {
  input {
    String modules = "blast emboss/6.6.0"
    File reference
    File consensusFasta
    String sample
    Int mem = 8
    Int timeout = 72
  }
 
  parameter_meta {
    consensusFasta: "Consensus FASTA file for local alignment"
    sample: "Sample name"
    mem: "Memory allocated for alignment task"
    timeout: "Timeout in hours for this task"
    modules: "Names and versions of modules needed for BLAST"
    reference: "Reference file"
  }

  command <<<
    set -euo pipefail

    # Suppress error for negative controls or samples with very little reads
    if blastn -query ~{consensusFasta} -subject ~{reference} \
    -word_size 28 -reward 1 -penalty -2 -dust no > ~{sample}_bl2seq_report.txt 2>error.log; then
        echo 'blastn completed successfully' 1>&2
    elif grep -q -F 'BLAST engine error: Warning: Sequence contains no data' error.log; then
        # Copy the message to STDERR, and exit without an error status
        cat error.log 1>&2
    else
        echo 'Unexpected error' 1>&2
        cat error.log 1>&2
        exit 1
    fi

   needle -asequence ~{reference} -bsequence ~{consensusFasta} -gapopen 10 -gapextend .5 -sid1 Sbjct: -sid2 Query: -outfile ~{sample}_needle_report.txt 
   >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File bl2seqReport = "~{sample}_bl2seq_report.txt"
    File needleReport = "~{sample}_needle_report.txt"
  }
}

task qcStats {
  input {
    String modules = "bedtools/2.27 samtools/1.9 picard/2.21.2 rstats/3.6"
    String sample
    File bam
    Int mem = 8
    Int timeout = 72
  }
  
   parameter_meta {
    bam: "Alignment file in BAM format"
    sample: "Sample name"
    mem: "Memory allocated for alignment task"
    timeout: "Timeout in hours for this task"
    modules: "Names and versions of modules needed for QC"
  }
  
  command <<<
    set -euo pipefail

    bedtools genomecov -ibam ~{bam} > ~{sample}.genomecvghist.txt

    bedtools genomecov -d -ibam ~{bam} > ~{sample}.genome.cvgperbase.txt

    samtools stats ~{bam} > ~{sample}.samstats.txt
 
    java -jar $PICARD_ROOT/picard.jar CollectInsertSizeMetrics I=~{bam} O=~{sample}.insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5
    
    printf "Count\tRead_Length\n">read_dist.txt.tmp
    samtools view ~{bam}| cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c >>read_dist.txt.tmp
    awk 'BEGIN{ OFS="\t"}{ print $1, $2}' read_dist.txt.tmp >~{sample}.read_dist.txt
   >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File genomecvgHist = "~{sample}.genomecvghist.txt"
    File genomecvgPerBase = "~{sample}.genome.cvgperbase.txt"
    File alignmentStats = "~{sample}.samstats.txt"
    File readDistStats = "~{sample}.read_dist.txt"
    File insertSizeStats = "~{sample}.insert_size_metrics.txt"
  }
}

task orfStats {
  input {
    String modules = "emboss/6.6.0"
    String sample
    String consensusFasta
    Int mem = 8
    Int timeout = 72
  }
  parameter_meta {
    consensusFasta: "Consensus file in FASTA format"
    sample: "Sample name"
    mem: "Memory allocated for alignment task"
    timeout: "Timeout in hours for this task"
    modules: "Names and versions of modules needed for ORF generation"
  }

  command <<<
    set -euo pipefail
   
    getorf --minsize 150 ~{consensusFasta} -out ~{sample}.orf
   >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File orf = "~{sample}.orf"
  }
}

task runReport {
  input {
    String modules = "rmarkdown/0.1 pt-report-tools/2.1"
    String sample
    String ext
    String flowcell
    String library
    String reference
    File consensusFasta
    String needleReport
    String bbmapLog
    String samStats
    String vcf
    String mvcf
    String readDistStats
    String insertSizeStats
    String construct
    String orf
    String json = "~{sample}.report.json"
    Int mem = 8
    Int timeout = 72
  }

   parameter_meta {
    sample: "Sample name"
    ext: "External sample name"
    flowcell: "Flowcell the sample was run on"
    library: "Library name of the sample"
    reference: "File path of the reference used for alignment"
    consensusFasta: "File path of the consensus generated"
    needleReport: "File path of the needle report"
    bbmapLog: "File path of the bbmap log file"
    samStats: "File path of the samtools stats file"
    vcf: "File path of the VCF file"
    mvcf: "File path of the mutation VCF file"
    readDistStats: "File path of the read distribution"
    insertSizeStats: "File path of the insert size distrubution"
    construct: "Name of the construct library"
    orf: "File path of the ORFs"
    mem: "Memory allocated for alignment task"
    json: "Summary stats in JSON format"
    timeout: "Timeout in hours for this task"
    modules: "Names and versions of modules needed"
  }

  command <<<
    set -euo pipefail
    perl $PT_REPORT_TOOLS_ROOT/info_for_rmarkdown.pl ~{bbmapLog} ~{samStats} ~{vcf} ~{needleReport} ~{reference} ~{insertSizeStats} ~{orf} ~{mvcf} >~{json}
    cp $PT_REPORT_TOOLS_ROOT/rmarkdownProvidence.Rmd .
    cp $PT_REPORT_TOOLS_ROOT/OICR.png .

    cp ~{readDistStats} readdist.txt
    Rscript -e "rmarkdown::render('./rmarkdownProvidence.Rmd', params=list(ext='~{ext}',construct='~{construct}',sample='~{sample}',library='~{library}',flowcell='~{flowcell}', refpath='~{reference}',refname='~{basename(reference)}', needle='~{needleReport}', readdist='~{readDistStats}',json='~{json}'), output_file='~{sample}.pdf')"
    tar -cvhf ~{sample}.scriptRmarkdown.tar.gz *json *pdf OICR.png *.txt *.Rmd script
   >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File rmarkdownReport = "~{sample}.pdf"
    File scriptRmarkdown = "~{sample}.scriptRmarkdown.tar.gz"
  }
} 
