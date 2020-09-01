version 1.0

workflow providencePipeline {
input {
    File fastq1
    File fastq2
    String samplePrefix
		String refe
  }

  parameter_meta {
    fastq1: "Read 1 fastq file, gzipped. Can be either targeted or whole transcriptome"
    fastq2: "Read 2 fastq file, gzipped. Can be either targeted or whole transcriptome."
    samplePrefix: "Prefix for output files"
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
      sample = samplePrefix
  }

  call variantCalling {
    input:
      sample = samplePrefix,
      bam = bwa.bwaMappedBam
  }

  output {
    File bam = bwa.bwaMappedBam
    File bai = bwa.bwaMappedBai
    File vcf = variantCalling.vcfFile
    File consensusFasta = variantCalling.consensusFasta
    File variantOnlyVcf = variantCalling.variantOnlyVcf
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
    ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=~{trimq}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File out1 = "~{sample}_qad_r1.fastq.gz"
    File out2 = "~{sample}_qad_r2.fastq.gz"
  }
}

task bwa {
  input {
    String modules = "bwa/0.7.12 samtools/1.9"
    File fastq1
    File fastq2
    String sample
		String ref
    Int mem = 12
    Int timeout = 72
    Int threads = 8
  }

  String bwaMappedBam_ = "~{sample}.bwa.bam"
  String bwaMappedBai_ = "~{sample}.bwa.bam.bai"

  command <<<
    set -euo pipefail

    #Align fastq files
    bwa mem -M -t 8 $ref \
    ~{fastq1} -2 ~{fastq2} | \
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
    String reference = "$HG38_BOWTIE_INDEX_ROOT/hg38_random_index"
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
