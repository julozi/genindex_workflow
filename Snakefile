GENOME_URLS = {
    'hg19': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'
}

rule download:
  output:
    temp("{name}/{prefix}/fasta_unmasked/chromFa.tar.gz")
  run:
    shell("mkdir -p {wildcards.name}/{wildcards.prefix}/fasta_unmasked")
    url = GENOME_URLS[wildcards.prefix]
    shell("wget -O {output} {url}")

rule untar:
  input:
    "{name}/{prefix}/fasta_unmasked/chromFa.tar.gz"
  output:
    "{name}/{prefix}/fasta_unmasked/{prefix}.fa"
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
  shell:
    """
      falist=$(ls {wildcards.name}/{wildcards.prefix}/fasta_unmasked/chr*.fa | egrep "chr[0-9MXY]+.fa" | sort -k 1.4,1.5n)
      falist=$(echo $falist | perl -ne 's/\\n/ /g; print $_')
      faToTwoBit $falist {output}
    """
