import pickle
import os
import random
import argparse
import pandas as pd
import utils
from pathlib import Path



def dicho(t,elt):
    s=0
    e=len(t)-1
    while(s<e):
        if(t[(s+e)//2]>elt):
            e=(s+e)//2-1
        elif(t[(s+e)//2]<elt):
            s=(s+e)//2+1
        else:
            return (s+e)//2
    return s



def intervals_overlap(a, b, c, d , proportion):
    # Calculate the overlap
    overlap_start = max(a, c)
    overlap_end = min(b, d)
    if (overlap_end-overlap_start) > proportion*min(b-a,d-c) :
        return True

    return False


#Look for a match between  segments
def look_for_match(segment,list_segments):
    for t in list_segments:
        if  intervals_overlap(segment[0],segment[1],t[0],t[1],0.8):
            return True
    return False

#Look for a no match between  segments
def look_for_close_match(segment,list_segments,min_dist):
    for t in list_segments:
        if  intervals_overlap(segment[0]-min_dist,segment[1]+min_dist,t[0],t[1],0):
            return True
    return False


#Compute the accuruacy of imputation on the region for all the SNPs
def compare_all_SNPs(t,seq_imputed,seq_original,out):
    i=dicho(seq_imputed[0],t[0])
    j=dicho(seq_original[0],t[0])

    while(i < len(seq_imputed[0])-1 and seq_imputed[0][i] < t[1]):
        while (j < len(seq_original[0])-1 and seq_original[0][j]<seq_imputed[0][i]):
            j+=1

        if (seq_original[0][j]==seq_imputed[0][i] ):     
            haplotype_original = seq_original[1][j]+seq_original[2][j]
            haplotype_imputed = seq_imputed[4][i][1]+seq_imputed[4][i][2]*2                   
            if (haplotype_original==0):
                mean_squared_error = (haplotype_original-haplotype_imputed)**2
                out.write(f'{1-seq_imputed[3][i]}\t{mean_squared_error}\n')
            elif (haplotype_original==2):
                mean_squared_error = (haplotype_original-haplotype_imputed)**2
                out.write(f'{seq_imputed[3][i]}\t{mean_squared_error}\n')
            elif (haplotype_original==1):
                mean_squared_error = (1-min(1,2-haplotype_imputed))**2
                out.write(f'{1-seq_imputed[3][i]}\t{mean_squared_error}\n')
                mean_squared_error = (1-min(1,haplotype_imputed))**2
                out.write(f'{seq_imputed[3][i]}\t{mean_squared_error}\n')
        i+=1
    return   


#Compare Imputation of RAF Archaic SNPs  vs Non Archaic SNPs, write the results in two files.
def compare_archaicSNPs(t,seq_imputed,seq_original,seq_outgroup,seq_archaic,fou_archaic,fou_non_archaic):
    i=dicho(seq_imputed[0],t[0])
    j=dicho(seq_original[0],t[0])
    k=dicho(seq_outgroup[0],t[0])
    l=dicho(seq_archaic[0],t[0])
    while(i < len(seq_imputed[0])-1 and seq_imputed[0][i] < t[1]):
        while (j < len(seq_original[0])-1 and seq_original[0][j]<seq_imputed[0][i]):
            j+=1
        while (k < len(seq_outgroup[0])-1 and seq_outgroup[0][k]<seq_imputed[0][i]):
            k+=1
        while (l < len(seq_archaic[0])-1 and seq_archaic[0][l]<seq_imputed[0][i]):
            l+=1
        if (seq_original[0][j]==seq_imputed[0][i] and seq_outgroup[0][k]==seq_imputed[0][i] and seq_archaic[0][l]==seq_imputed[0][i]):
            gen=-1
            if (seq_archaic[1][l] not in seq_outgroup[1][k] and seq_archaic[2][l] not in seq_outgroup[1][k]):
                gen=seq_archaic[1][l]
            if gen==-1:
                current_file=fou_non_archaic
            else:
                current_file=fou_archaic
            haplotype_original = seq_original[1][j]+seq_original[2][j]
            haplotype_imputed = seq_imputed[4][i][1]+seq_imputed[4][i][2]*2                   
            if (haplotype_original==0):
                mean_squared_error = (haplotype_original-haplotype_imputed)**2
                current_file.write(f'{1-seq_imputed[3][i]}\t{mean_squared_error}\n')
            elif (haplotype_original==2):
                mean_squared_error = (haplotype_original-haplotype_imputed)**2
                current_file.write(f'{seq_imputed[3][i]}\t{mean_squared_error}\n')
            elif (haplotype_original==1):
                mean_squared_error = (1-min(1,2-haplotype_imputed))**2
                current_file.write(f'{1-seq_imputed[3][i]}\t{mean_squared_error}\n')
                mean_squared_error = (1-min(1,haplotype_imputed))**2
                current_file.write(f'{seq_imputed[3][i]}\t{mean_squared_error}\n')
        i+=1
    return    


#Compare Imputation of RAF Archaic SNPs  vs Non Archaic SNPs, write the results in two files.
#We ask the SNPs to be private to the archaic genome we consider, absent from the other archaic genome.
def compare_private_archaicSNPs(t,seq_imputed,seq_original,seq_outgroup,seq_archaic_ref, seq_archaic_outgroup,fou_archaic,fou_non_archaic):
    i=dicho(seq_imputed[0],t[0])
    j=dicho(seq_original[0],t[0])
    k=dicho(seq_outgroup[0],t[0])
    l=dicho(seq_archaic_ref[0],t[0])
    m=dicho(seq_archaic_outgroup[0],t[0])
    while(i < len(seq_imputed[0])-1 and seq_imputed[0][i] < t[1]):
        while (j < len(seq_original[0])-1 and seq_original[0][j]<seq_imputed[0][i]):
            j+=1
        while (k < len(seq_outgroup[0])-1 and seq_outgroup[0][k]<seq_imputed[0][i]):
            k+=1
        while (l < len(seq_archaic_ref[0])-1 and seq_archaic_ref[0][l]<seq_imputed[0][i]):
            l+=1
        while (m < len(seq_archaic_outgroup[0])-1 and seq_archaic_outgroup[0][m]<seq_imputed[0][i]):
            m+=1
        if (seq_original[0][j]==seq_imputed[0][i] and seq_outgroup[0][k]==seq_imputed[0][i] and seq_archaic_ref[0][l]==seq_imputed[0][i]):
            gen=-1
            if (seq_archaic_ref[1][l] not in seq_outgroup[1][k] and seq_archaic_ref[2][l] not in seq_outgroup[1][k]):
                gen=seq_archaic_ref[1][l]
            if gen==-1 or (seq_archaic_outgroup[0][m]==seq_imputed[0][i] and (seq_archaic_outgroup[1][m]==gen or seq_archaic_outgroup[2][m]==gen) ):
                current_file=fou_non_archaic
            else:
                current_file=fou_archaic
            haplotype_original = seq_original[1][j]+seq_original[2][j]
            haplotype_imputed = seq_imputed[4][i][1]+seq_imputed[4][i][2]*2                   
            if (haplotype_original==0):
                mean_squared_error = (haplotype_original-haplotype_imputed)**2
                current_file.write(f'{1-seq_imputed[3][i]}\t{mean_squared_error}\n')
            elif (haplotype_original==2):
                mean_squared_error = (haplotype_original-haplotype_imputed)**2
                current_file.write(f'{seq_imputed[3][i]}\t{mean_squared_error}\n')
            elif (haplotype_original==1):
                mean_squared_error = (1-min(1,2-haplotype_imputed))**2
                current_file.write(f'{1-seq_imputed[3][i]}\t{mean_squared_error}\n')
                mean_squared_error = (1-min(1,haplotype_imputed))**2
                current_file.write(f'{seq_imputed[3][i]}\t{mean_squared_error}\n')
        i+=1
    return    




 # We look at all the SNPs in the chromosome and write their r² score depending on if the SNP is archaic or not.
def imputation_quality_archaic_SNPs():
    coverage = args.coverage
    chrom = args.chrom

    path = Path(f'../out/concordance/{coverage}/archaic')
    path.mkdir(parents=True, exist_ok=True)
    path = Path(f'./out/concordance/{coverage}/nonArchaic')
    path.mkdir(parents=True, exist_ok=True)

    outgroup = pickle.load(open(f"../out/pickledAfricans/chr{chrom}/outgroup","rb"))

    directory_archaic = f'../out/pickledArchaic/chr{chrom}'
    for filename in os.listdir(directory_archaic):
        with open(os.path.join(directory_archaic, filename),'rb') as fp:
            out_archaic =  open(f'../out/concordance/{coverage}/archaic/{coverage}_{filename}_{chrom}_archaic','w')
            out_nonarchaic =  open(f'../out/concordance/{coverage}/nonArchaic/{coverage}_{filename}_{chrom}_nonArchaic','w')
            out_archaic.write(f'RAF\tmean_squared_error\n')
            out_nonarchaic.write(f'RAF\tmean_squared_error\n')
            seq_archaic = pickle.load(fp)

            directory = f'../out/pickledData/{coverage}/chr{chrom}'
            
            for filename in os.listdir(directory):
                try:
                    with open(os.path.join(directory, filename),'rb') as fp:
                        seq_imputed = pickle.load(fp)
                        try:
                            with open(f'../out/pickledData/original/chr{chrom}/{filename}', 'rb') as fp:
                                seq_original = pickle.load(fp)
                                compare_archaicSNPs([int(seq_imputed[0][0]),int(seq_imputed[0][-1])],seq_imputed,seq_original,outgroup,\
                                                    seq_archaic,out_archaic,out_nonarchaic)
                        except IOError:
                            print(f'The reference non imputed individual {filename} could not be opened, make sure you created it using --read-vcf and --coverage original')
                except IOError:
                    print(f'error cannot open')
            out_archaic.close()
            out_nonarchaic.close()


 # We look at all the SNPs in the chromosome and write their imputation quality score depending on if they are in an archaic region or not.
def imputation_quality_archaic_segments():
    coverage = args.coverage
    chrom = args.chrom
    path = Path(f'../out/concordance/{coverage}/unshared_allSNPs')
    path.mkdir(parents=True, exist_ok=True)
    path = Path(f'../out/concordance/{coverage}/shared_allSNPs')
    path.mkdir(parents=True, exist_ok=True)
    path = Path(f'../out/concordance/{coverage}/nonintrogressed_allSNPs')
    path.mkdir(parents=True, exist_ok=True)

    out_archaicnew =  open(f'../out/concordance/{coverage}/unshared_allSNPs/{coverage}_{chrom}_archaicnew_allSNPs','w')
    out_archaicshared =  open(f'../out/concordance/{coverage}/shared_allSNPs/{coverage}_{chrom}_archaicshared_allSNPs','w')
    out_nonintrogressed =  open(f'../out/concordance/{coverage}/nonintrogressed_allSNPs/{coverage}_{chrom}_nonintrogressed_allSNPs','w')
    out_archaicnew.write(f'RAF\tmean_squared_error\n')
    out_archaicshared.write(f'RAF\tmean_squared_error\n')
    out_nonintrogressed.write(f'RAF\tmean_squared_error\n')

    directory = f'../out/pickledData/{coverage}/chr{chrom}'

    segments_imputed = {}
    segments_original = {}
    df = pd.read_csv(args.path_map, dtype={'Chrom': str})
    for _, row in df.iterrows():
        if row['Coverage'] == coverage and row['Chrom'] == chrom:
            if row['individual'] not in segments_imputed:
                segments_imputed[row['individual']] = []
            segments_imputed[row['individual']].append([int(row['start']),int(row['end'])])
        if row['Coverage'] == 'Original' and row['Chrom'] == chrom:
            if row['individual'] not in segments_original:
                segments_original[row['individual']] = []
            segments_original[row['individual']].append([int(row['start']),int(row['end'])])

    for filename in os.listdir(directory):
        try:
            with open(os.path.join(directory, filename),'rb') as fp:
                seq_imputed = pickle.load(fp)
                with open(f'../out/pickledData/original/chr{chrom}/{filename}', 'rb') as fp:
                    seq_original = pickle.load(fp)
 

                for segment in segments_imputed[filename]:
                    if look_for_match(segment,segments_original[filename]):
                        compare_all_SNPs(segment,seq_imputed,seq_original,out_archaicshared)
                    else:
                        compare_all_SNPs(segment,seq_imputed,seq_original,out_archaicnew)
                        
                    #We look at a random non introgressed region
                    offset = random.randint(-10_000_000, 10_000_000)
                    t=[segment[0]+offset,segment[1]+offset]

                    while t[0]<= seq_imputed[0][0] or t[1] >= seq_imputed[0][-1] or look_for_close_match(t,segments_imputed[filename],10_000) \
                        or  look_for_close_match(t,segments_original[filename],10_000):
                        offset = random.randint(-10_000_000, 10_000_000)
                        t = [segment[0]+offset,segment[1]+offset]

                    compare_all_SNPs(t,seq_imputed,seq_original,out_nonintrogressed)


        except IOError:
            print(f'error cannot open')
    out_archaicnew.close()
    out_archaicshared.close()
    out_nonintrogressed.close()


if __name__ == "__main__":
    '''
    If pickledArchaic and pickledOutgroup have not been created through the LAI script, run first:

        imputationAccuracy.py --read-vcf-archaic  --myfile *name_vcf_archaic*.vcf.gz 
        imputationAccuracy.py --read-outgroup  --directory *name of the directory containing the outgroup vcfs*

        Notes:
        --read-vcf-archaic assumes that the vcf given as input contains all chromosomes.
        If multiple archaic individuals are in the same vcf, all the archaic individuals will be separately used to define archaic SNPs
        It is possible to run the command multiple time if different archaic individuals are present in different vcfs

        The outgroup is used to define archaic SNPs
        The name of the outgroup vcf is hardcoded in the script and can be directly downloaded from:
        https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
        It is assumed that the directory given as input for --read-outgroup contains those files
        The name of the african samples to be extracted from these vcf are also directly written in utils.py
        It is possible to change these to use a different outgroup. 

        The output directory will be created at ../out

    Then for each coverage, run:

        imputationAccuracy.py --read-vcf  --file *name_vcf_sample*.vcf.gz --coverage *name of the dataset (coverage)*

        --read-vcf assumes that the vcf given as input contains all chromosomes.
        If multiple individuals are in the same vcf, the imputation accuracy will be computed for all individuals
        It is possible to run the command multiple time for the same coverage if different individuals are present in different vcfs
        <!> For the non imputed dataset (which will serve as a reference to compute imputation accuracy), the coverage should be exactly named "original"
        <!> For each coverage the same individual should have the same name, if not either rename it first in the vcf or the output files in ../out/pickledData/*name of the dataset (coverage)*

    Finally, run:

        imputationAccuracy.py --imputation-quality-archaic-snps --coverage *name of the dataset (coverage)* --chrom *chromosome*

        The output files will be in '../out/concordance/{coverage}/archaic' and '../out/concordance/{coverage}/nonArchaic'
        The first folder contains imputation accuracy for archaic SNPs, the second for non archaic SNPs

    '''
    parser = argparse.ArgumentParser(description="")

    parser.add_argument('--read-vcf-archaic', action='store_true', help='Enable reading archaic VCF file.')
    parser.add_argument('--read-vcf', action='store_true', help='Enable reading sample VCF file with RAF and GP.')
    parser.add_argument('--read-outgroup', action='store_true', help='Enable reading outgroup data.')
    parser.add_argument('--imputation-quality-archaic-segments', action='store_true', help='Compute the imputation quality depending if a SNP is in introgressed segment or not.')
    parser.add_argument('--imputation-quality-archaic-snps', action='store_true', help='Compute the imputation quality depending if a SNP is archaic or not.')

    parser.add_argument("--path_map", type=str, required=False,
                        help="path to the introgression map download from: FILL")

    parser.add_argument('--coverage', required=False, help='Name of the dataset (used to store output files)\
                        Important: if this is the non imputed reference dataset, call it "original"')

    parser.add_argument("--chrom", type=str, required=False,
                        choices=[str(i) for i in range(1, 23)],
                        help="runs the imputation quality analysis for the selected chromosome")
    parser.add_argument('--file', required=False, help='Path to input VCF file')
    parser.add_argument('--directory', required=False, help='Directory containing input VCF files')


    args = parser.parse_args()


    if args.file is None:
        pass
    else:
        file = args.file

    if args.coverage is None:
        pass
    else:
        coverage = args.coverage

    if args.directory is None:
        pass
    else:
        directory = args.directory


    if args.imputation_quality_archaic_segments:
        imputation_quality_archaic_segments()
    elif args.imputation_quality_archaic_snps:
        imputation_quality_archaic_SNPs()
    elif args.read_vcf_archaic:
        for chrom in range(1,23):
            res = utils.read_vcf(file,chrom)
            for i,ind in enumerate(res[1]):
                path = Path(f'../out/pickledArchaic/chr{chrom}')
                path.mkdir(parents=True, exist_ok=True)
                filehandler = open(f'../out/pickledArchaic/chr{chrom}/{ind}',"wb")
                pickle.dump(res[0][i],filehandler)
                filehandler.close()
    elif args.read_vcf:
        for chrom in range(1,23):
            if coverage == 'original':
                res = utils.read_vcf(file,chrom)
            else:
                res = utils.read_vcf_raf_gp(file,chrom)
            for i,ind in enumerate(res[1]):
                path = Path(f'../out/pickledData/{coverage}/chr{chrom}')
                path.mkdir(parents=True, exist_ok=True)
                filehandler = open(f'../out/pickledData/{coverage}/chr{chrom}/{ind}',"wb")
                pickle.dump(res[0][i],filehandler)
                filehandler.close()
    elif args.read_outgroup:
        for chrom in range(1,23):
            res = utils.read_vcf_sum(f'{directory}/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz',chrom)
            path = Path(f'../out/pickledAfricans/chr{chrom}')
            path.mkdir(parents=True, exist_ok=True)
            filehandler = open(f'../out/pickledAfricans/chr{chrom}/outgroup',"wb")
            pickle.dump(res,filehandler)
            filehandler.close()

