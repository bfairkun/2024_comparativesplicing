import argparse
import pyBigWig
import pandas as pd

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Append mean BigWig signal to BED intervals.")
    parser.add_argument('-b', '--bedfile', required=True, help='Input BED file')
    parser.add_argument('-w', '--bigwigfile', required=True, help='Input BigWig file')
    parser.add_argument('-o', '--outfile', required=True, help='Output BED file with mean signal appended')
    return parser.parse_args(args)

def calculate_mean_signal(row, bw):
    try:
        signal = bw.values(row[0], row[1], row[2], numpy=True)
        if len(signal) == 0:
            return 0.0
        return signal.mean()
    except RuntimeError as e:
        # Handle invalid interval bounds or other runtime errors
        print(f"Error processing interval {row[0]}:{row[1]}-{row[2]}: {e}")
        return 0.0

def process_bed_file(bed_file, bigwig_file, output_file):
    bw = pyBigWig.open(bigwig_file)
    chunk_size = 10000  # Adjust chunk size based on memory constraints and performance
    
    with pd.read_csv(bed_file, sep='\t', header=None, chunksize=chunk_size) as reader:
        with open(output_file, 'w') as out_file:
            for chunk in reader:
                chunk[3] = chunk.apply(lambda row: calculate_mean_signal(row, bw), axis=1)
                chunk.to_csv(out_file, sep='\t', header=False, index=False, mode='a')

    bw.close()

if __name__ == "__main__":
    if hasattr(sys, 'ps1'):
        # main("-i /project2/yangili1/bjf79/2024_comparativesplicing/code/GenomeFiles/Human_ensemblv75/Reference.gtf -o scratch/Human_ensemblv75.Reannotated.gtf -fa /project2/yangili1/bjf79/ReferenceGenomes/Human_ensemblv75/Reference.fa -v -translation_approach A".split(' '))
        main("-b conservation/human.AS_segments.bed -w conservation/PhyloP.hg19.bw -o scratch/test.bed".split(' '))
        # main("-i /project2/yangili1/bjf79/2024_comparativesplicing/code/GenomeFiles/Human_ensemblv75/Reference.gtf -o scratch/Human_ensemblv75.Reannotated.gtf -fa /project2/yangili1/bjf79/ReferenceGenomes/Human_ensemblv75/Reference.fa -v -translation_approach A".split(' '))
    main()


