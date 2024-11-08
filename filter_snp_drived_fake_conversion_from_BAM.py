#filter the snp from BAM file to avoid preturbing the nascent reads identification
import os
import pysam
import argparse  

# get the raw nascent BAM file
def raw_nascent_reads(bam_file,threads):
    raw_nas_bam=bam_file.replace('.bam','_raw_nas.bam')
    with pysam.AlignmentFile(bam_file, "rb") as raw_bam, pysam.AlignmentFile(raw_nas_bam, "wb", template=raw_bam) as nas_bam:
        for read in raw_bam:
            if str(read.cigarstring) != 'None' and int(read.get_tag('Yf'))>0:
                nas_bam.write(read)
            else:
                continue
    
    samtools = '/gpfs/home/qihansong/biosoft/samtools/bin/samtools'
    os.system(f"{samtools} sort --threads {threads} {raw_nas_bam} -o {raw_nas_bam}") #sort the bam file

#filter the snp drived fake nascent reads to get the final nascent reads containing BAM file
def get_nascent_bam(raw_nascent_bam, snp_dict, nascent_bam_file, threads):
    #getting the T2C snp
    snp_positions={}
    with open(snp_dict,'r') as snp:
        snp_positions['non essential']=[]
        for line in snp:
            chromosome=line.strip().split('\t')[0]
            pos=line.strip().split('\t')[1].replace('[','').replace(']','').split(', ')
            snp_positions[chromosome]=pos
    #reading the bam file
    with pysam.AlignmentFile(raw_nascent_bam, "rb") as in_bam, pysam.AlignmentFile(nascent_bam_file, "wb", template=in_bam) as out_bam:
        snpcontain={}
        for read in in_bam:
            snpcount=0
            readid=read.query_name #get the read id
            chr='chr'+str(int(read.reference_id)+1) #get the chr info
            if int(read.reference_id)<22:
                chr=chr
            if chr=='chr23':
                chr='chrX'
            elif chr=='chr24':
                chr='chrY'
            elif chr=='chr25':
                chr='chrM'
            else:
                chr='non essential' #the T2C mutation of other chromosome will not be taken into account
            match_pos=[str(pos) for pos in read.get_reference_positions()] #get all the matched position
            Yf_tag=int(read.get_tag('Yf')) #get the original conversed base count
            chr_snp=[str(snp) for snp in snp_positions[chr]] #all the snp of the given chr
            for pos in match_pos: #the snp count of the read
                if pos in chr_snp:
                    snpcount+=1
                else:
                    continue
            new_Yf_tag=Yf_tag-snpcount
            if snpcount==0:
                out_bam.write(read)
            elif snpcount>0 and new_Yf_tag>0:
                snpcontain[readid]=[snpcount,Yf_tag,new_Yf_tag]
                read.set_tag('Yf',new_Yf_tag)
                out_bam.write(read)
            elif snpcount>0:
                snpcontain[readid]=[snpcount,Yf_tag,new_Yf_tag]
            else:
                continue
        samtools = '/gpfs/home/qihansong/biosoft/samtools/bin/samtools'
        os.system(f"{samtools} sort --threads {threads} {nascent_bam_file} -o {nascent_bam_file}") #sort the bam file
    return snpcontain

def main():
    # define the paramaters
    parser = argparse.ArgumentParser(description='take the snp into count for getting the real nascent reads from the hisat3N outputed BAM file, all the paramaters were required.')
    parser.add_argument('-s','--snp', type=str, required=True, help='the path of snp dictionary file, outputed by in house script selecting_snp_from_vcf_file.py')
    parser.add_argument('-p','--threads', type=int, required=True, help='the path of the output BAM file')
    parser.add_argument('-i','--input_bam', type=str, required=True, help='the path of input BAM file from the hisat3N alignment output')
    parser.add_argument('-o','--output_bam', type=str, required=True, help='the path of the output BAM file')
    args = parser.parse_args()
    snp_file=args.snp
    thread=args.threads
    bam_file=args.input_bam
    nascent_bam_file=args.output_bam

    raw_nascent_reads(bam_file=bam_file,threads=thread)
                
    raw_nas_bam=bam_file.replace('.bam','_raw_nas.bam')
    snp_contain_reads=get_nascent_bam(raw_nascent_bam=raw_nas_bam,snp_dict=snp_file,nascent_bam_file=nascent_bam_file,threads=thread)
    snp_reads=os.path.basename(nascent_bam_file).replace('_nascent.bam','_snp_contained_reads.txt')
    with open(snp_reads,'w') as snp_reads_op:
        snp_reads_op.write('read_id'+'\t'+'snp_count'+'\t'+'raw_conversion_count'+'\t'+'real_conversion_count'+'\n')
        for key,infos in snp_contain_reads.items():
            info='\t'.join([str(value) for value in infos])
            snp_reads_op.write(f"{key}\t{info}\n")

if __name__ == "__main__":
    main()
