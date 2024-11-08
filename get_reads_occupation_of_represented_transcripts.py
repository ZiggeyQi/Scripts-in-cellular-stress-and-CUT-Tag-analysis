import os
import pysam
import argparse

# get the represent transcript to further analysis
def get_represent_transcripts(raw_gtf,tss,tes,genetype,out_dir):
    final_gtf=os.path.join(out_dir,os.path.basename(raw_gtf)).replace('.gtf',f"_{genetype}_longest_transcript.gtf")
    seleted_gene_transcript=os.path.join(out_dir,os.path.basename(raw_gtf)).replace('.gtf',f"_{genetype}_longest_transcript.txt")
    geneid=[]
    transcript_info=dict()
    final_gene_id=[]
    final_transcript_id=[]
    with open(raw_gtf,'r') as rgtf, open(final_gtf,'w') as fgtf, open(seleted_gene_transcript,'w') as op_gt:
        lines = rgtf.readlines()
        # get the genes based on the defined length of the least(tss down plus tes up)
        tss_down=abs(int(tss.split(',')[1]))
        tes_up=abs(int(tes.split(',')[0]))
        print(f"the {genetype} genes, whose transcripts longer than {tss_down+tes_up}, will be selected")
        for line1 in lines:
            if not line1.startswith('##'):
                if line1.strip().split('\t')[2]=='gene' and line1.strip().split('gene_type "')[1].split('"')[0]==genetype:
                     info=line1.strip().split('\t')
                     id=line1.strip().split('gene_id "')[1].split('"')[0]
                     genelength=int(info[4])-int(info[3])
                     if genelength > tss_down+tes_up:
                        geneid.append(id)
                     else:
                        continue 
                else:
                    continue           
            else:
                fgtf.write(line1)
        print(f"raw selected genes were {len(geneid)}")
        ## get the represnted transcript of the longest transcript of each gene
        # get transcript info
        for line2 in lines:
            if not line2.startswith('##') and line2.strip().split('\t')[2]=='transcript':
                if any(gene_id in line2 for gene_id in geneid):
                    gene_ID=line2.strip().split('gene_id "')[1].split('"')[0]
                    trans_info=line2.strip().split('\t')
                    trans_id=line2.strip().split('transcript_id "')[1].split('"')[0]
                    trans_length=int(trans_info[4])-int(trans_info[3])
                    if trans_info[6]=='+':
                        tss_point=int(trans_info[3])
                    elif trans_info[6]=='-':
                        tss_point=int(trans_info[4])
                    else:
                        continue
                    if gene_ID not in transcript_info:
                        transcript_info[gene_ID]=[{trans_id:(tss_point,trans_length)}]
                    else:
                        transcript_info[gene_ID].append({trans_id:(tss_point,trans_length)})
                else:
                    continue
            else:
                continue
        print(f"selected genes were {len(transcript_info)}")
        # seleting the longest transcript
        op_gt.write('gene_id\ttranscript_id\tlenght\n')
        for key, list_of_dicts in transcript_info.items(): 
            sorted_dicts = sorted(list_of_dicts, key=lambda d: max(v[1] for v in d.values()), reverse=True)
            for key1,tup in sorted_dicts[0].items():
                if int(tup[1])> tss_down+tes_up:
                    final_gene_id.append(key)
                    final_transcript_id.append(key1)
                    op_gt.write(f"{key}\t{key1}\t{tup[1]}\n")
                else:
                    continue
        print(f"The final santified genes were {len(final_gene_id)} and the corresponding represented transcripts were {len(final_transcript_id)}")
        # outputing the final gtf file
        for line in lines:
            if not line.startswith('##') and line.strip().split('\t')[2]=='gene':
                gene=line.strip().split('gene_id "')[1].split('"')[0]
                if gene in final_gene_id:
                    fgtf.write(line)
            elif not line.startswith('##'):
                transcript=line.strip().split('transcript_id "')[1].split('"')[0]
                if transcript in final_transcript_id:
                    fgtf.write(line)
            else:
                continue



# counting the reads located in distinct regions and the mean depth of the regions, as well the elongation ratio
def get_count(bam,chr,start,end):
    count=0
    for read in bam.fetch(chr,start,end):
        count= count + 1
    return count

def get_mean_depth(bam,chr,start,end):
    coverage=bam.count_coverage(contig=chr,start=start,stop=end)
    coverage=str(coverage).replace("array('L', [",'').replace("])","").replace(")","").replace("(","").replace(", ","")
    mean_depth=sum(list(int(i) for i in coverage))/(end-start)
    return mean_depth

def get_info_matrix(bam_file,gtf_file,promoter_region,polyadenylation_region,output_dir):
    out_file=os.path.join(output_dir,os.path.basename(bam_file).replace('.bam','_occupation_information.txt'))
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(gtf_file,'r') as gtf, open(out_file,'w') as op_file:
        header='\t'.join(['geneID','transcriptID','gene_symbol','gene_count','promoter_count','genebody_count',
                          'polyadenylation_count','gene_depth','promoter_depth','genebody_depth','polyadenylation_depth',
                          'elongation_ratio','termination_ratio'])
        op_file.write(f"{header}\n")
        promoter=[abs(int(value1)) for value1 in promoter_region.split(',')]
        polyadenylation=[abs(int(value2)) for value2 in polyadenylation_region.split(',')]
        min_length=promoter[1]+polyadenylation[0]
        for line in gtf:
            if not line.startswith('##') and line.strip().split('\t')[2]=='transcript':
                gene_id=line.strip().split('gene_id "')[1].split('"')[0]
                gene_name=line.strip().split('gene_name "')[1].split('"')[0]
                trans_id=line.strip().split('transcript_id "')[1].split('"')[0]
                chr=line.strip().split('\t')[0]
                length=int(line.strip().split('\t')[4])-int(line.strip().split('\t')[3])
                if length > min_length:
                    if line.strip().split('\t')[6]=='+':
                        tss_site=int(line.strip().split('\t')[3])
                        tes_site=int(line.strip().split('\t')[4])
                        if tss_site-promoter[0] >0:
                            promoter_count=get_count(bam,chr,tss_site-promoter[0],tss_site+promoter[1])
                            genebody_count=get_count(bam,chr,tss_site+promoter[1],tes_site-polyadenylation[0])
                            polyadenylation_count=get_count(bam,chr,tes_site-polyadenylation[0],tes_site+polyadenylation[1])
                            gene_count=get_count(bam,chr,tss_site-promoter[0],tes_site+polyadenylation[1])
                            promoter_depth=get_mean_depth(bam,chr,tss_site-promoter[0],tss_site+promoter[1])
                            genebody_depth=get_mean_depth(bam,chr,tss_site+promoter[1],tes_site-polyadenylation[0])
                            polyadenylation_depth=get_mean_depth(bam,chr,tes_site-polyadenylation[0],tes_site+polyadenylation[1])
                            gene_depth=get_mean_depth(bam,chr,tss_site-promoter[0],tes_site+polyadenylation[1])
                            if genebody_depth==0 and promoter_depth==0:
                                elongation_ratio='NA'
                            else:
                                elongation_ratio=genebody_depth/(promoter_depth+genebody_depth)
                            if genebody_depth==0 and polyadenylation_depth==0:
                                termination_ratio='NA'
                            else:
                                termination_ratio=polyadenylation_depth/(genebody_depth+polyadenylation_depth)
                            op_file.write(f"{gene_id}\t{trans_id}\t{gene_name}\t{gene_count}\t{promoter_count}\t{genebody_count}\t\
                                        {polyadenylation_count}\t{gene_depth}\t{promoter_depth}\t{genebody_depth}\t\
                                            {polyadenylation_depth}\t{elongation_ratio}\t{termination_ratio}\n")
                        else:
                            promoter_count=get_count(bam,chr,1,tss_site+promoter[1])
                            genebody_count=get_count(bam,chr,tss_site+promoter[1],tes_site-polyadenylation[0])
                            polyadenylation_count=get_count(bam,chr,tes_site-polyadenylation[0],tes_site+polyadenylation[1])
                            gene_count=get_count(bam,chr,1,tes_site+polyadenylation[1])
                            promoter_depth=get_mean_depth(bam,chr,1,tss_site+promoter[1])
                            genebody_depth=get_mean_depth(bam,chr,tss_site+promoter[1],tes_site-polyadenylation[0])
                            polyadenylation_depth=get_mean_depth(bam,chr,tes_site-polyadenylation[0],tes_site+polyadenylation[1])
                            gene_depth=get_mean_depth(bam,chr,1,tes_site+polyadenylation[1])
                            if genebody_depth==0 and promoter_depth==0:
                                elongation_ratio='NA'
                            else:
                                elongation_ratio=genebody_depth/(promoter_depth+genebody_depth)
                            if genebody_depth==0 and polyadenylation_depth==0:
                                termination_ratio='NA'
                            else:
                                termination_ratio=polyadenylation_depth/(genebody_depth+polyadenylation_depth)
                            op_file.write(f"{gene_id}\t{trans_id}\t{gene_name}\t{gene_count}\t{promoter_count}\t{genebody_count}\t\
                                        {polyadenylation_count}\t{gene_depth}\t{promoter_depth}\t{genebody_depth}\t\
                                            {polyadenylation_depth}\t{elongation_ratio}\t{termination_ratio}\n")
                    elif line.strip().split('\t')[6]=='-':
                        tss_site=int(line.strip().split('\t')[4])
                        tes_site=int(line.strip().split('\t')[3])
                        if tes_site-polyadenylation[1] >0:
                            promoter_count=get_count(bam,chr,tss_site-promoter[1],tss_site+promoter[0])
                            genebody_count=get_count(bam,chr,tes_site+polyadenylation[0],tss_site-promoter[1])
                            polyadenylation_count=get_count(bam,chr,tes_site-polyadenylation[1],tes_site+polyadenylation[0])
                            gene_count=get_count(bam,chr,tes_site-polyadenylation[1],tss_site+promoter[0])
                            promoter_depth=get_mean_depth(bam,chr,tss_site-promoter[1],tss_site+promoter[0])
                            genebody_depth=get_mean_depth(bam,chr,tes_site+polyadenylation[0],tss_site-promoter[1])
                            polyadenylation_depth=get_mean_depth(bam,chr,tes_site-polyadenylation[1],tes_site+polyadenylation[0])
                            gene_depth=get_mean_depth(bam,chr,tes_site-polyadenylation[1],tss_site+promoter[0])
                            if genebody_depth==0 and promoter_depth==0:
                                elongation_ratio='NA'
                            else:
                                elongation_ratio=genebody_depth/(promoter_depth+genebody_depth)
                            if genebody_depth==0 and polyadenylation_depth==0:
                                termination_ratio='NA'
                            else:
                                termination_ratio=polyadenylation_depth/(genebody_depth+polyadenylation_depth)
                            op_file.write(f"{gene_id}\t{trans_id}\t{gene_name}\t{gene_count}\t{promoter_count}\t{genebody_count}\t\
                                        {polyadenylation_count}\t{gene_depth}\t{promoter_depth}\t{genebody_depth}\t\
                                            {polyadenylation_depth}\t{elongation_ratio}\t{termination_ratio}\n")
                        else:
                            promoter_count=get_count(bam,chr,tss_site-promoter[1],tss_site+promoter[0])
                            genebody_count=get_count(bam,chr,tes_site+polyadenylation[0],tss_site-promoter[1])
                            polyadenylation_count=get_count(bam,chr,1,tes_site+polyadenylation[0])
                            gene_count=get_count(bam,chr,1,tss_site+promoter[0])
                            promoter_depth=get_mean_depth(bam,chr,tss_site-promoter[1],tss_site+promoter[0])
                            genebody_depth=get_mean_depth(bam,chr,tes_site+polyadenylation[0],tss_site-promoter[1])
                            polyadenylation_depth=get_mean_depth(bam,chr,1,tes_site+polyadenylation[0])
                            gene_depth=get_mean_depth(bam,chr,1,tss_site+promoter[0])
                            if genebody_depth==0 and promoter_depth==0:
                                elongation_ratio='NA'
                            else:
                                elongation_ratio=genebody_depth/(promoter_depth+genebody_depth)
                            if genebody_depth==0 and polyadenylation_depth==0:
                                termination_ratio='NA'
                            else:
                                termination_ratio=polyadenylation_depth/(genebody_depth+polyadenylation_depth)
                            op_file.write(f"{gene_id}\t{trans_id}\t{gene_name}\t{gene_count}\t{promoter_count}\t{genebody_count}\t\
                                        {polyadenylation_count}\t{gene_depth}\t{promoter_depth}\t{genebody_depth}\t\
                                            {polyadenylation_depth}\t{elongation_ratio}\t{termination_ratio}\n")
                    else:
                        continue
                else:
                    continue
    print('Analysis finished -_- -_- -_-')                    


def main():
    # define the paramaters
    parser = argparse.ArgumentParser(description=f"This script will get the reads occupation details of the defined gene region \
                                     based on the aligned BAM and the annotation GTF file")
    parser.add_argument('-b','--bam_file', type=str, required=True, help='The path of input BAM file')
    parser.add_argument('-g','--gtf', type=str, required=True, help=f"Ensembl formated GTF file, we defaultly chosing the longest transcript\
                        to analysis, you can set -i parameter as TRUE to provide your interest transcript containing GTF file for analysis")
    parser.add_argument('-p','--promoter_region', type=str, default='-1000,500', help=f"Your interest promoter region, \
                        provided as -1000,500 the defined region was refered to TSS, default(-1000,500)")
    parser.add_argument('-t','--polyadenylation_region', type=str, default='-500,500', help=f"Your interest polyadenylation(termination) region, \
                        provided as -500,500 the defined region was refered to TES, default(-500,500)")
    parser.add_argument('-y','--gene_type', type=str, default='protein_coding', help=f"Your interest gene type, default(protein_coding)")
    parser.add_argument('-i','--provide_interest_gtf', type=str, default='FALSE', help=f"Set as TRUE to provide your interest transcript\
                         containing GTF file, only TRUE or FALSE, default(FALSE)")
    parser.add_argument('-o','--output_dir', type=str, required=True, help='The outputting directory')
    args = parser.parse_args()
    bam_file=args.bam_file
    gtf_file=args.gtf
    promoter=args.promoter_region
    termination=args.polyadenylation_region
    genetype=args.gene_type
    interest_gtf=args.provide_interest_gtf
    output=args.output_dir

    if interest_gtf=='FALSE':
        get_represent_transcripts(raw_gtf=gtf_file,tss=promoter,tes=termination,genetype=genetype,out_dir=output)
        final_gtf=os.path.join(output,os.path.basename(gtf_file)).replace('.gtf',f"_{genetype}_longest_transcript.gtf")
        get_info_matrix(bam_file=bam_file,gtf_file=final_gtf,promoter_region=promoter,polyadenylation_region=termination,output_dir=output)
    elif interest_gtf=='TRUE':
        get_info_matrix(bam_file=bam_file,gtf_file=gtf_file,promoter_region=promoter,polyadenylation_region=termination,output_dir=output)
    else:
        print('you are not properly provide the info for which gtf file will be used, check the --provide_interest_gtf parameter please')
    

if __name__ == "__main__":
    main()
