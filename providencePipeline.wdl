version 1.0

workflow providencePipeline {
input {
    File fastq1
    File fastq2
    String samplePrefix
    File? primerBed
    File? panelBed
    Int? readCount
  }

  parameter_meta {
    fastq1: "Read 1 fastq file, gzipped. Can be either targeted or whole transcriptome"
    fastq2: "Read 2 fastq file, gzipped. Can be either targeted or whole transcriptome."
    samplePrefix: "Prefix for output files"
  }

  meta {
  	author: "Angie Mosquera"
  	email: "amosquera@oicr.on.ca"
  	description: "Classify samples as being either SARS-CoV-2 positive or negative, identify the strain of virus, and produce statistics about the mapping."
  	dependencies: [
    	{
    	    name: "bbmap/38.75",
    	    url: "https://sourceforge.net/projects/bbmap/"
    	},
        {
    	    name: "bowtie2/2.3.5.1",
    	    url: "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zi"
    	},
        {
    	    name: "samtools/1.9",
    	    url: "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2"
    	},
        {
    	    name: "kraken2/2.0.8",
    	    url: "https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads"
    	},
        {
    	    name: "ivar/1.0",
    	    url: "https://github.com/andersen-lab/ivar"
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
    	    name: "bedtools/2.27",
    	    url: "https://github.com/arq5x/bedtools2/releases/tag/v2.27.1"
    	},
        {
    	    name: "blast/2.8.1",
    	    url: "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz"
    	},
        {
    	    name: "spades/3.14.0",
    	    url: "http://cab.spbu.ru/files/release3.14.0/manual.html"
    	}
  	]
  }

  call bbMap {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
      sample = samplePrefix
  }

  call bowtie2HumanDepletion {
    input:
      fastq1 = bbMap.out1,
      fastq2 = bbMap.out2,
      sample = samplePrefix
  }

  call kraken2 {
    input:
      fastq1 = bowtie2HumanDepletion.out1,
      fastq2 = bowtie2HumanDepletion.out2,
      sample = samplePrefix
  }

  call bowtie2Sensitive {
    input:
      fastq1 = bowtie2HumanDepletion.out1,
      fastq2 = bowtie2HumanDepletion.out2,
      sample = samplePrefix
  }

  if (defined(primerBed)) {
    call articTrimming {
      input:
        bam = bowtie2Sensitive.bamFile,
        sample = samplePrefix,
        primerBed = select_first([primerBed])
    }
  }

  call variantCalling {
    input:
      sample = samplePrefix,
      bam = select_first([articTrimming.sortedBam, bowtie2Sensitive.bamFile])
  }

  call qcStats {
    input:
      sample = samplePrefix,
      bam = select_first([articTrimming.sortedBam, bowtie2Sensitive.bamFile]),
      hostMappedBam = bowtie2HumanDepletion.hostMappedBam,
      panelBed = panelBed
  }

  call blast2ReferenceSequence {
    input:
      consensusFasta = variantCalling.consensusFasta,
      sample = samplePrefix
  }

  if (!defined(readCount)) {
    call generateReadCount {
      input:
        fastq = fastq1
    }
  }

  call spadesGenomicAssembly {
    input:
      fastq1 = bowtie2HumanDepletion.out1,
      fastq2 = bowtie2HumanDepletion.out2,
      sample = samplePrefix,
      readCount = select_first([readCount, generateReadCount.readCount])
  }

  output {
    File hostRemovedR1Fastq = bowtie2HumanDepletion.out1
    File hostRemovedR2Fastq = bowtie2HumanDepletion.out2
    File hostMappedBam = bowtie2HumanDepletion.hostMappedBam
    File hostMappedBai = bowtie2HumanDepletion.hostMappedBai
    File taxonomicClassification = kraken2.out
    File bam = bowtie2Sensitive.bamFile
    File bai = bowtie2Sensitive.baiFile
    File? primertrimSortedBam = articTrimming.sortedBam
    File? primertrimSortedBai = articTrimming.sortedBai
    File vcf = variantCalling.vcfFile
    File consensusFasta = variantCalling.consensusFasta
    File variantOnlyVcf = variantCalling.variantOnlyVcf
    File bl2seqReport = blast2ReferenceSequence.bl2seqReport
    File? cvgHist = qcStats.cvgHist
    File genomecvgHist = qcStats.genomecvgHist
    File genomecvgPerBase = qcStats.genomecvgPerBase
    File hostMappedAlignmentStats = qcStats.hostMappedAlignmentStats
    File hostDepletedAlignmentStats = qcStats.hostDepletedAlignmentStats
    File spades = spadesGenomicAssembly.sampleSPAdes
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

task bowtie2HumanDepletion {
  input {
    String modules = "bowtie2/2.3.5.1 samtools/1.9 hg38-bowtie-index/2.3.5.1"
    File fastq1
    File fastq2
    String reference = "$HG38_BOWTIE_INDEX_ROOT/hg38_random_index"
    String sample
    Int mem = 12
    Int timeout = 72
    Int threads = 8
  }

  String hostMappedBam_ = "~{sample}.host.mapped.bam"
  String hostMappedBai_ = "~{sample}.host.mapped.bam.bai"

  command <<<
    set -euo pipefail

    #Align fastq files to hg38 & only keep unmapped

    bowtie2 -p ~{threads} -x ~{reference} \
    -1 ~{fastq1} -2 ~{fastq2} \
    --un-conc-gz ~{sample}_host_removed_r%.fastq.gz | \
    samtools view -b | \
    samtools sort - -o ~{hostMappedBam_}

    samtools index ~{hostMappedBam_}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File out1 = "~{sample}_host_removed_r1.fastq.gz"
    File out2 = "~{sample}_host_removed_r2.fastq.gz"
    File hostMappedBam = "~{hostMappedBam_}"
    File hostMappedBai = "~{hostMappedBai_}"
  }
}

task kraken2 {
  input {
    String modules = "kraken2/2.0.8 kraken2-database/1"
    File fastq1
    File fastq2
    String kraken2DB = "$KRAKEN2_DATABASE_ROOT/"
    String sample
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    kraken2 --paired ~{fastq1} ~{fastq2} \
    --db ~{kraken2DB} \
    --report ~{sample}.kreport2.txt
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File out = "~{sample}.kreport2.txt"
  }
}

task bowtie2Sensitive {
  input {
    String modules = "bowtie2/2.3.5.1 sars-covid-2-polymasked-bowtie-index/2.3.5.1 samtools/1.9"
    File fastq1
    File fastq2
    String sample
    String sarsCovidIndex = "$SARS_COVID_2_POLYMASKED_BOWTIE_INDEX_ROOT/MN908947.3.mask"
    Int mem = 8
    Int timeout = 72
    Int threads = 4
  }

  String bamFile_ = "~{sample}.bam"
  String baiFile_ = "~{sample}.bam.bai"

  command <<<
    set -euo pipefail

    #Align reads to reference

    bowtie2 --sensitive-local -p ~{threads} -x ~{sarsCovidIndex} \
    -1 ~{fastq1} -2 ~{fastq2} | samtools view -b | samtools sort - -o ~{bamFile_}

    samtools index ~{bamFile_}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File bamFile = "~{bamFile_}"
    File baiFile = "~{baiFile_}"
  }
}

task articTrimming {
  input {
    String modules = "ivar/1.0 bedtools"
    File bam
    File primerBed
    String sample
    Boolean? allowNoprimer
    Int mem = 8
    Int timeout = 72
  }

  String primertrim = "~{sample}.primertrim"
  String primertrimBam = "~{sample}.primertrim.bam"
  String sortedBam_ = "~{sample}.primertrim.sorted.bam"
  String sortedBai_ = "~{sample}.primertrim.sorted.bam.bai"

  command <<<
    set -euo pipefail

    ivar trim -i ~{bam} -b ~{primerBed} -p ~{primertrim} ~{true="-e" false="" allowNoprimer}

    samtools sort ~{primertrimBam} -o ~{sortedBam_}

    samtools index ~{sortedBam_}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File? sortedBam = "~{sortedBam_}"
    File? sortedBai = "~{sortedBai_}"
  }

  parameter_meta {
    bam: "Host depleted and aligned to mn908947 bam file."
    primerBed: "Bed file used to trim the primers off of the bam sequences."
    allowNoprimer: "Allow reads that don't have primer sequence? Ligation prep = false, nextera = true."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

task variantCalling {
  input {
    String modules = "bcftools/1.9 samtools/1.9 vcftools/0.1.16 seqtk/1.3 sars-covid-2-polymasked/mn908947.3"
    File bam
    String sample
    String sarsCovidRef = "$SARS_COVID_2_POLYMASKED_ROOT/MN908947.3.mask.fasta"
    Int mem = 8
    Int timeout = 72
  }

  String vcfName = "~{sample}.vcf"
  String fastaName = "~{sample}.consensus.fasta"
  String variantOnlyVcf_ = "~{sample}.v.vcf"

  command <<<
    set -euo pipefail

    #Call consensus sequence
    samtools mpileup -aa -uf ~{sarsCovidRef} ~{bam} | \
    bcftools call --ploidy 1 -Mc | tee -a ~{vcfName} | \
    vcfutils.pl vcf2fq -d 10 | \
    seqtk seq -A - | sed '2~2s/[actg]/N/g' > ~{fastaName}

    bcftools mpileup -a "INFO/AD,FORMAT/DP,FORMAT/AD" \
    -d 8000 -f ~{sarsCovidRef} ~{bam} | \
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

task qcStats {
  input {
    String modules = "bedtools samtools/1.9"
    String sample
    File? panelBed
    File bam
    File hostMappedBam
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    if [ -f "~{panelBed}" ]; then
      bedtools coverage -hist -a ~{panelBed} -b ~{bam} > ~{sample}.cvghist.txt
    fi

    bedtools genomecov -ibam ~{bam} > ~{sample}.genomecvghist.txt

    bedtools genomecov -d -ibam ~{bam} > ~{sample}.genome.cvgperbase.txt

    samtools stats ~{hostMappedBam} > ~{sample}.host.mapped.samstats.txt

    samtools stats ~{bam} > ~{sample}.samstats.txt
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File? cvgHist = "~{sample}.cvghist.txt"
    File genomecvgHist = "~{sample}.genomecvghist.txt"
    File genomecvgPerBase = "~{sample}.genome.cvgperbase.txt"
    File hostMappedAlignmentStats = "~{sample}.host.mapped.samstats.txt"
    File hostDepletedAlignmentStats = "~{sample}.samstats.txt"
  }
}

task blast2ReferenceSequence {
  input {
    String modules = "blast sars-covid-2-polymasked/mn908947.3"
    String reference = "$SARS_COVID_2_POLYMASKED_ROOT/MN908947.3.mask.fasta"
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

task generateReadCount {
  input {
    String modules = ""
    File fastq
    Int mem = 4
    Int timeout = 4
  }

  command <<<
    zcat ~{fastq} | sed -n '1~4p' | wc -l
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Int readCount = read_int(stdout())
  }
}

task spadesGenomicAssembly {
  input {
    String modules = "spades/3.14.0"
    Int readCount
    Int minReads = 100
    File fastq1
    File fastq2
    String sample
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    mkdir ~{sample}.SPAdes

    if [ "~{readCount}" -gt "~{minReads}" ]; then
      rnaspades.py --pe1-1 ~{fastq1} --pe1-2 ~{fastq2} -o ~{sample}.SPAdes
    else
      echo 'Not enough reads to run SPAdes genomic assembly.' 1>&2
    fi

    tar cf - ~{sample}.SPAdes | gzip --no-name > ~{sample}.SPAdes.tar.gz
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File sampleSPAdes = "~{sample}.SPAdes.tar.gz"
  }
}