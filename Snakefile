GENOME_URLS = {
  'hg19': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz',
  'hg38': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz',
  'mm10': 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz',
  'mm9': 'http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz',
  'sacCer3': 'http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz',
}

localrules: make_all

rule download:
  output:
    temp("{name}/{prefix}/fasta_unmasked/chromFa.tar.gz")
  resources:
    mem_mb=2048
  run:
    shell("mkdir -p {wildcards.name}/{wildcards.prefix}/fasta_unmasked")
    url = GENOME_URLS[wildcards.prefix]
    shell("wget -O {output} {url}")

rule untar:
  input:
    "{name}/{prefix}/fasta_unmasked/chromFa.tar.gz"
  output:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa"
  resources:
    mem_mb=2048
  shell:
    """
      cd {wildcards.name}/{wildcards.prefix}/fasta_unmasked
      tar xvzf chromFa.tar.gz
      cat `ls | egrep "chr[0-9MXY]+.fa" | sort -k 1.4,1.5n` > {wildcards.prefix}.fa
    """

rule twobit:
  input:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa"
  output:
    "{name}/{prefix}/2bit/{prefix}.2bit"
  conda:
    "env.yaml"
  resources:
    mem_mb=2048
  shell:
    """
      falist=$(ls {wildcards.name}/{wildcards.prefix}/fasta_unmasked/chr*.fa | egrep "chr[0-9MXY]+.fa" | sort -k 1.4,1.5n)
      falist=$(echo $falist | perl -ne 's/\\n/ /g; print $_')
      faToTwoBit $falist {output}
    """

rule fasta_index:
  input:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa"
  output:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa.fai"
  conda:
    "env.yaml"
  resources:
    mem_mb=5120
  shell:
    """
      samtools faidx {input}
    """

rule fasta_dict:
  input:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa"
  output:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa.dict"
  conda:
    "env.yaml"
  resources:
    mem_mb=5120
  shell:
    """
      picard CreateSequenceDictionary R={input} O={output}
    """

rule bowtie1_indexes:
  input:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa"
  output:
    ebwt1="{name}/{prefix}/unmasked_bowtie1_indexes/{prefix}.1.ebwt"
    ebwt2="{name}/{prefix}/unmasked_bowtie1_indexes/{prefix}.2.ebwt"
    ebwt3="{name}/{prefix}/unmasked_bowtie1_indexes/{prefix}.3.ebwt"
    ebwt4="{name}/{prefix}/unmasked_bowtie1_indexes/{prefix}.4.ebwt"
    ebwtrev1="{name}/{prefix}/unmasked_bowtie1_indexes/{prefix}.rev.1.ebwt"
    ebwtrev2="{name}/{prefix}/unmasked_bowtie1_indexes/{prefix}.rev.2.ebwt"
  conda:
    "env.yaml"
  resources:
    mem_mb=10240
  shell:
    """
      bowtie-build {input} {wildcards.name}/{wildcards.prefix}/unmasked_bowtie1_indexes/{wildcards.prefix}
    """

rule bowtie2_indexes:
  input:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa"
  output:
    bt21="{name}/{prefix}/unmasked_bowtie2_indexes/{prefix}.1.bt2"
    bt22="{name}/{prefix}/unmasked_bowtie2_indexes/{prefix}.2.bt2"
    bt23="{name}/{prefix}/unmasked_bowtie2_indexes/{prefix}.3.bt2"
    bt24="{name}/{prefix}/unmasked_bowtie2_indexes/{prefix}.4.bt2"
    bt2rev1="{name}/{prefix}/unmasked_bowtie2_indexes/{prefix}.rev.1.bt2"
    bt2rev2="{name}/{prefix}/unmasked_bowtie2_indexes/{prefix}.rev.2.bt2"
  conda:
    "env.yaml"
  resources:
    mem_mb=10240
  shell:
    """
      bowtie2-build {input} {wildcards.name}/{wildcards.prefix}/unmasked_bowtie2_indexes/{wildcards.prefix}
    """

rule bwa_indexes:
  input:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa"
  output:
    amb="{name}/{prefix}/unmasked_bwa_indexes/{prefix}.amb"
    ann="{name}/{prefix}/unmasked_bwa_indexes/{prefix}.ann"
    bwt="{name}/{prefix}/unmasked_bwa_indexes/{prefix}.bwt"
    pac="{name}/{prefix}/unmasked_bwa_indexes/{prefix}.pac"
    sa="{name}/{prefix}/unmasked_bwa_indexes/{prefix}.sa"
  conda:
    "env.yaml"
  resources:
    mem_mb=5120
  shell:
    """
      bwa index -a bwtsw {input} -p {wildcards.name}/{wildcards.prefix}/unmasked_bwa_indexes/{wildcards.prefix}
    """

rule make_all:
  input:
    bwa="{name}/{prefix}/unmasked_bwa_indexes/{prefix}.amb"
    bowtie2="{name}/{prefix}/unmasked_bowtie2_indexes/{prefix}.1.bt2"
    bowtie1="{name}/{prefix}/unmasked_bowtie1_indexes/{prefix}.1.ebwt"
    fadict="{name}/{prefix}/fasta_unmasked/{prefix}.fa.dict"
    fai="{name}/{prefix}/fasta_unmasked/{prefix}.fa.fai"
    twobit="{name}/{prefix}/2bit/{prefix}.2bit"
  output:
    "{name}/{prefix}/fasta_unmasked/done"
  shell:
    """
      touch {output}
    """
