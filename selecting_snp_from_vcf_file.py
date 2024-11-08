import argparse 
# get the T2C mutation snp info
def get_snp_positions(snp_file,mutation): #the input vcf file was outputed by VarScan2, fileformat=VCFv4.3
    snp_positions = {}
    T2C_snp_dictionary=snp_file.replace('.vcf','_T2C_snp_dictionary.txt')
    T2C_snp_position=snp_file.replace('.vcf','_T2C_snp_position.txt')
    with open(snp_file, 'r') as snp, open(T2C_snp_dictionary,'w') as tcsnp_dict, open(T2C_snp_position,'w') as tcsnp_pos:
        for line in snp:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                vaf=float(line.split('%:')[0].split(':')[-1]) #the Variant allele frequency value
                pvalue=float(line.split('%:')[1].split(':')[0])
                ref=parts[3] #reference base
                alt=parts[4] #altered base
                target_ref=mutation.split(',')[0] #the base to be mutated
                target_alt=mutation.split(',')[1] #the base will be mutated to
                if ref==target_ref and alt==target_alt and vaf>90 and pvalue<1e-5: #selectting T>C mutation and the calculated vaf should higher than 90% and pvalue lower than 1e-5
                    if parts[0] not in snp_positions:
                        snp_positions[parts[0]]=[int(parts[1])]
                    else:
                        snp_positions[parts[0]].append(int(parts[1]))
                else:
                    continue
            else:
                continue
        for key in snp_positions:
            snp_positions[key].sort()
        for key,snp in snp_positions.items():
            tcsnp_dict.write(f"{key}\t{snp}\n")
            for pos in snp:
                tcsnp_pos.write(f"{key}\t{pos}\n")

def main():
    # define the paramaters
    parser = argparse.ArgumentParser(description='selecting the single mutation snp from the VarScan2 outputed vcf file, fileformat=VCFv4.3.')
    parser.add_argument('-s','--snp', type=str, required=True, help='the path of vcf formated snp file, the target snp file will output to the input file path of the vcf file')
    parser.add_argument('-m','--mutation', type=str, required=True, help='the reference base and the base mutated to be, Uppercase and comma separated, eg: T,C')
    args = parser.parse_args()
    snp_file=args.snp
    mutation=args.mutation

    get_snp_positions(snp_file=snp_file,mutation=mutation)

if __name__ == "__main__":
    main()
