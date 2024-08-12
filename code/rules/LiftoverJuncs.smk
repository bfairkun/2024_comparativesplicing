
rule Download_Chain_to_hg38:
    output:
        "ChainFiles/ToHg38/{GenomeName}.chain.gz"
    params:
        link = lambda wildcards: ChainFiles.loc[wildcards.GenomeName]['Genome_to_Hg38_Chain_Link']
    shell:
        """
        wget -O {output} {params.link}
        """

rule Download_Chain_from_hg38:
    output:
        "ChainFiles/FromHg38/{GenomeName}.chain.gz"
    params:
        link = lambda wildcards: ChainFiles.loc[wildcards.GenomeName]['Hg38_to_Genome_Chain_Link']
    shell:
        """
        wget -O {output} {params.link}
        """

rule Gather_Chains:
    input:
        expand("ChainFiles/FromHg38/{GenomeName}.chain.gz", GenomeName = ChainFiles.index),
        expand("ChainFiles/ToHg38/{GenomeName}.chain.gz", GenomeName = ChainFiles.index),
