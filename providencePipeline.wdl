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
    fastq2: "Read 2 fastq file, gzipped. Can be either targeted or whole transcriptome."
    samplePrefix: "Prefix for output files"
    sampleLibary: "Library name of the transcriptome"
    sampleFlowcell: "Flowcell name"
    sampleExt: "External name provided for this library"
    ref: "Reference to be used for alignment"
    construct: "Library name of the construct used"
  }

  meta {
    author: "Yogi Sundaravadanam"
    email: "ysundaravadanam@oicr.on.ca"
    description: "."
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
      }
    ]
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
  
   call blast2ReferenceSequence {
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
      library = sampleLibrary,
      flowcell = sampleFlowcell,
      reference = ref,
      construct = sampleConstruct,
      bl2seqReport = blast2ReferenceSequence.bl2seqReport,
      consensusFasta = variantCalling.consensusFasta,
      vcf = variantCalling.vcfFile,
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
    File bl2seqReport = blast2ReferenceSequence.bl2seqReport
    File alignmentStats = qcStats.alignmentStats
    File readDistStats = qcStats.readDistStats
    File insertSizeStats = qcStats.insertSizeStats
    File bbMaplog = bbMap.bbMapLog
    File report = runReport.rmarkdownReport
    File reportScript = runReport.scriptRmarkdown
    File readDistStats = qcStats.readDistStats
    File insertSizeStats = qcStats.insertSizeStats
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

  String bwaMappedBam_ = "~{sample}.bwa.bam"
  String bwaMappedBai_ = "~{sample}.bwa.bam.bai"

  command <<<
    set -euo pipefail

    #index the reference
    bwa index ~{reference}
    
    #Align fastq files
    bwa mem -M -t 8 ~{reference} \
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

task blast2ReferenceSequence {
  input {
    String modules = "blast"
    File reference
    File consensusFasta
    String sample
    Int mem = 8
    Int timeout = 72
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
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File bl2seqReport = "~{sample}_bl2seq_report.txt"
  }
}

task qcStats {
  input {
    String modules = "bedtools samtools/1.9 picard/2.21.2"
    String sample
    File bam
    Int mem = 8
    Int timeout = 72
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
    String modules = "rmarkdown/0.1 pt-report-tools/1.0"
    String sample
    String flowcell
    String library
    String reference
    File consensusFasta
    String bl2seqReport
    String bbmapLog
    String samStats
    String vcf
    String readDistStats
    String insertSizeStats
    String construct
    String orf
    String json = "~{sample}.report.json"
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail
    perl $PT_REPORT_TOOLS_ROOT/info_for_rmarkdown.pl ~{bbmapLog} ~{samStats} ~{vcf} ~{bl2seqReport} ~{reference} ~{insertSizeStats} ~{orf} >~{json}
    cp $PT_REPORT_TOOLS_ROOT/rmarkdownProvidence.Rmd .
    cp ~{readDistStats} readdist.txt
    Rscript -e "rmarkdown::render('./rmarkdownProvidence.Rmd', params=list(ext='mPVT-S5000 B0011',construct='~{construct}',sample='~{sample}',library='~{library}',flowcell='~{flowcell}', refpath='~{reference}', blast='~{bl2seqReport}', readdist='~{readDistStats}',json='~{json}'), output_file='~{sample}.rmarkdown.pdf')"
    tar -cvhf ~{sample}.scriptRmarkdown.tar.gz *json *pdf *.txt *.Rmd script
   >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File rmarkdownReport = "~{sample}.rmarkdown.pdf"
    File scriptRmarkdown = "~{sample}.scriptRmarkdown.tar.gz"
  }
} 
