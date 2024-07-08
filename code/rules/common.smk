import pandas as pd
import os


# try/except useful for running this script in isolation in interactive shell
# for debugging
try:
    samples = pd.read_csv(config["samples"],sep='\t', index_col=0)
    STAR_genomes = pd.read_csv(config["STAR_genomes"],sep='\t', index_col=0)
except (NameError, KeyError) as NameOrKeyError:
    samples = pd.read_csv("config/samples.tsv",sep='\t', index_col=0)
    STAR_genomes = pd.read_csv("config/STAR_Genome_List.tsv",sep='\t', index_col=0)

MazinGenomes = [g for g in STAR_genomes.index if g not in ["Lamprey_ensemblv_112"]]

# Add code for function definitions and other things that must be defined prior
# to rest of workflow (eg custom snakemake input functions)
def GetIndexingParams(wildcards):
    if STAR_genomes.loc[wildcards.GenomeName]['ChromLargerThan512Mbp'] == 'T':
        return '-C'
    else:
        return ''

def GetIndexSuffix(wildcards):
    if STAR_genomes.loc[wildcards.GenomeName]['ChromLargerThan512Mbp'] == 'T':
        return 'csi'
    else:
        return 'tbi'
