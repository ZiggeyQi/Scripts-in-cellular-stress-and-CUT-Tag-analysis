import os
import argparse

#editting the gtf file to get the protein coding gene
def extract_gtf(out_dir, gtf_file, gene_type):
    types=gene_type.strip().split(',')
    typename='_'.join(t for t in types)
    output_gtf_file=os.path.join(out_dir,os.path.basename(gtf_file).replace('.gtf',f"_{typename}.gtf"))
    with open(gtf_file,'r') as in_gtf, open(output_gtf_file,'w') as op_gtf:
        for line in in_gtf:
            if line.startswith('##'): #preserve the header lines
                op_gtf.write(line)
                continue
            for index in range(len(types)):
                if line.split('gene_type "')[1].startswith(f'{types[index]}";'):
                    op_gtf.write(line)
                else:
                    continue

def main():
    # define the paramaters
    parser = argparse.ArgumentParser(description='get the interest gene type annotation information(gtf) from the gtf file')
    parser.add_argument('-g','--gtf', type=str, required=True, help='the path of GTF file')
    parser.add_argument('-t','--gene_type', type=str, required=True, help='the interest gene type to define the transcription orientaton of the reads, if more than one gene types provided, comma-separated, eg: protein_coding,lncRNA')
    parser.add_argument('-o','--output_dir', type=str, required=True, help='the dir of the outputting gtf file storage')
    args = parser.parse_args()
    gtf_file=args.gtf
    tpye=args.gene_type
    out_dir=args.output_dir

    extract_gtf(out_dir=out_dir,gtf_file=gtf_file,gene_type=tpye)
    
if __name__ == "__main__":
    main()
