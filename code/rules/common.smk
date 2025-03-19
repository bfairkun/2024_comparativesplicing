import pandas as pd
import os
import filecmp
import pathlib


# try/except useful for running this script in isolation in interactive shell
# for debugging
try:
    samples = pd.read_csv(config["samples"],sep='\t', index_col=0)
    STAR_genomes = pd.read_csv(config["STAR_genomes"],sep='\t', index_col=0)
except (NameError, KeyError) as NameOrKeyError:
    samples = pd.read_csv("config/samples.tsv",sep='\t', index_col=0)
    STAR_genomes = pd.read_csv("config/STAR_Genome_List.tsv",sep='\t', index_col=0)

MazinGenomes = [g for g in STAR_genomes.index if g not in ["Lamprey_ensemblv_112"]]

SampleGenomes = samples['STARGenomeName'].unique()

CordosMoreiraGenomes = [g for g in SampleGenomes if g not in ['StarletSeaAnemone_RefSeq_GCF_932526225.1','SalpingoecaRosetta_ensemblv_59']]
CordosMoreiraGenomes_dict = dict(zip([Genome.split('_')[0] for Genome in CordosMoreiraGenomes], CordosMoreiraGenomes))

CordosoMoreira_df = pd.read_csv("config/Cordoso_Moreira_SampleList.tsv", sep='\t', index_col=0)
CordosoMoreira_df['formatted_column'] = CordosoMoreira_df.apply(
    lambda row: f"{row['fastq_aspera'].replace('fasp.sra.ebi.ac.uk:', '')} {os.getcwd()}/CordosoMoreira_Fastq/{row.name}.fastq.gz",
    axis=1
)
CordosoMoreira_df['formatted_column'].to_csv("config/batch.txt", index=False, header=False)

CordosoMoreira_df['JuncFileGroupForParralelization'] = CordosoMoreira_df.groupby('ID_Species').cumcount() % 10 + 1

Extra_Gtfs = pd.read_csv("config/CordosoGenomes_Extra_Gtfs.tsv", sep='\t', index_col=0)
GenomesWithRefSeq = Extra_Gtfs[Extra_Gtfs['RefSeqGtf'].notna()].index
GenomesWithEnsembl = Extra_Gtfs[Extra_Gtfs['EnsemblGtf'].notna()].index

ChainFiles = pd.read_csv("config/ChainFiles.tsv", sep='\t', index_col=0)

GTEx_juncFiles = pd.read_csv("config/GTEx_juncFileList.tsv", sep='\t', index_col=0)

Cordoso_contrasts = pd.read_csv("config/CordosoTimeSeriesContrasts.tsv", sep='\t', index_col=0)
Cordoso_TissueContrasts = pd.read_csv("config/CordosoTissueContrasts.tsv", sep='\t', index_col=0)
Cordoso_SamplesNoHumanEmbryo = pd.read_csv("../data/Cordoso_SampleList.WhichAreNonHumanEmbryo.tsv.gz", sep='\t', index_col=0)

maf_filelist = pd.read_csv("../data/multiz100way_maf_md5sum.txt", sep='\s+', names=["md5sum", "fn"])
maf_filelist["fn_base"] = maf_filelist['fn'].str.replace(".maf.gz$", "", regex=True)
maf_filelist['link'] = maf_filelist["fn"].apply(lambda fn: f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/{fn}")


def has_differences(dcmp):
    """
    https://stackoverflow.com/questions/4187564/recursively-compare-two-directories-to-ensure-they-have-the-same-files-and-subdi
    """
    try:
        differences = dcmp.left_only + dcmp.right_only + dcmp.diff_files
        if differences:
            return True
        return any([has_differences(subdcmp) for subdcmp in dcmp.subdirs.values()])
    except NotADirectoryError:
        return True

# Add code for function definitions and other things that must be defined prior
# to rest of workflow (eg custom snakemake input functions)

def CreateSymlinksOfDir1ContentsIntoDir2(Dir1, Dir2):
    """
    helper function to create symlinks for scripts... Imagine a snakemake rule
    defined in a module workflow with shell directive...
    shell: 'Rscript myRscriptWithRelativeFilepathRelativeToModuleWorkflow.R SomeMoreArgs'
    ...This rule will only work if
    myRscriptWithRelativeFilepathRelativeToModuleWorkflow.R exists in the
    workdir for the main Snakefile.
    """
    Dir1_sanitized = Dir1.rstrip("/") + "/"
    Dir2_sanitized = Dir2.rstrip("/") + "/"
    for filepath in pathlib.Path(Dir1_sanitized).glob('[!.]*'):
        module_script_file = os.path.abspath(filepath)
        new_script_link = Dir2_sanitized + os.path.basename(filepath)
        try:
            os.symlink(module_script_file, new_script_link)
            print(f'Making link: {new_script_link}->{module_script_file}', file=sys.stderr)
        except FileExistsError:
            if filecmp.cmp(module_script_file, new_script_link) or not has_differences(filecmp.dircmp(module_script_file, new_script_link)):
                # print(f'{new_script_link}->{module_script_file} already exists', file=sys.stderr)
                pass
            # elif os.readlink(new_script_link)==module_script_file:
            #     pass
            else:
                print(f'not making link, fix clashing file names: {new_script_link}->{module_script_file}', file=sys.stderr)

# Add code for function definitions and other things that must be defined prior
# to rest of workflow (eg custom snakemake input functions)
def GetIndexingParams(wildcards):
    if STAR_genomes.loc[wildcards.GenomeName]['ChromLargerThan512Mbp'] in ['T', 'True', 'TRUE', True]:
        return '--csi'
    else:
        return ''

def GetIndexSuffix(wildcards):
    if STAR_genomes.loc[wildcards.GenomeName]['ChromLargerThan512Mbp'] in ['T', 'True', 'TRUE', True]:
        return 'csi'
    else:
        return 'tbi'

def ExpandAllSamplesInFormatStringFromGenomeNameWildcard_NoHumanEmbryo(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=Cordoso_SamplesNoHumanEmbryo.loc[(Cordoso_SamplesNoHumanEmbryo['GenomeName']==wildcards.GenomeName) & (Cordoso_SamplesNoHumanEmbryo['NoHumanEmbryo'].isin(['T' 'True', 'TRUE', True])) ].index.unique())
    return InputFunctionToReturn

def ExpandAllSamplesInFormatStringFromGenomeNameWildcard_NoEmbryo(FormattedString):
    def InputFunctionToReturn(wildcards):
        return expand(FormattedString, sample=Cordoso_SamplesNoHumanEmbryo.loc[(Cordoso_SamplesNoHumanEmbryo['GenomeName']==wildcards.GenomeName) & (Cordoso_SamplesNoHumanEmbryo['NoEmbryo'].isin(['T' 'True', 'TRUE', True])) ].index.unique())
    return InputFunctionToReturn
