
GENOME_URLS = {
    'hg19': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'
}

rule download:
    output:
        {name}/{prefix}/fasta_unmasked/chromFa.tar.gz
    run:
        import subprocess
        wget_cmd = ["wget", "-O", output, GENOME_URLS[wildcards.prefix]]
        subprocess.run(wget_cmd)

rule untar:
    input:
        {name}/{prefix}/fasta_unmasked/chromFa.tar.gz
    ouput:
        {name}/{prefix}/fasta_unmasked/{prefix}.fa
    shell:
        cd {wildcards.name}/{wildcards.prefix}/fasta_unmasked
        tar xvzf chromFa.tar.gz
        cat `ls | egrep "chr[0-9MXY]+.fa" | sort -k 1.4,1.5n` > {output}
        rm chr*
