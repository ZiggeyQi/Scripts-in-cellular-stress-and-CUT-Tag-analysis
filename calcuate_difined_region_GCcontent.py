from Bio import SeqIO
import pandas as pd
import argparse
import os

# loading the genome
def load_genome(genome_fasta):
    genome = {}
    for seq_record in SeqIO.parse(genome_fasta, "fasta"):
        genome[seq_record.id] = seq_record.seq
    return genome

# extract the defined sequence
def extract_sequence(genome, chr, site, strand, up, down):
    if strand == '+':
        start = max(0, site - up - 1)  # site was given as 1-based
        end = site + down - 1
    elif strand == '-':
        start = site - down -1
        end = site + up -1
    return genome[chr][start:end]

# calculate GC content 
def calculate_gc_content(seq):
    g = seq.count('G')
    c = seq.count('C')
    return (g + c) / len(seq) * 100

def main():
    # define the paramaters
    parser = argparse.ArgumentParser(description=f"This script will extract the sequence and its GC content \
                                     based on the defined up and down region refered to the given location")
    parser.add_argument('-g','--genome_file', type=str, required=True, help='The path of input genome file, fasta formated')
    parser.add_argument('-s','--site_information', type=str, required=True, help=f"The site information file, at least contain 3 column\
                        named chr, site, strand, Tab separated formated")
    parser.add_argument('-l','--location_up_down', type=str, required=True, help=f"Your interest region, up and down refered to \
                        the given site, provided as 1000,500")
    parser.add_argument('-o','--output_dir', type=str, required=True, help='The outputting directory')
    args = parser.parse_args()
    genome_file=args.genome_file
    site_file=args.site_information
    location=args.location_up_down
    output=args.output_dir

    # loading genome
    genome = load_genome(genome_file)
    
    # get the info file for processing
    gene_list = pd.read_csv(site_file, sep='\t')
    # get the defined up and down location
    site_up=abs(int(location.split(',')[0]))
    site_down=abs(int(location.split(',')[1]))

    results = []
    main_chr=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    for index, row in gene_list.iterrows():
        if row['chr'] in main_chr:
            chr = 'chr'+ row['chr']
        elif row['chr']=='MT':
            continue
        else:
            chr = row['chr']
        site = int(row['site'])
        strand = row['strand']
        defined_seq = extract_sequence(genome=genome, chr=chr, site=site, strand= strand, up=site_up, down=site_down)
        gc_content = calculate_gc_content(defined_seq)
        results.append({
            'gene': row['gene'],
            'transcript': row['transcript'],
            'entrezgene_id': row['entrezgene_id'],
            'chr': chr,
            'site': site,
            'strand': strand,
            'transcript_length': row['transcript_length'],
            'gene_gc_content': row['gene_gc_content'],
            'defined_seq': str(defined_seq),
            'promoter_gc_content': gc_content
        })
    
    # save results
    results_df = pd.DataFrame(results)
    results_file=os.path.join(output,f"hg38.p14_tss_{location}_gc_content.txt")
    results_df.to_csv(results_file, sep='\t', index=False)    

if __name__ == "__main__":
    main()
