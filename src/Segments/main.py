import utils
import os
import pickle
from pathlib import Path
import argparse
import HMM

parser = argparse.ArgumentParser(description="Control options for processing genomic data.")


''' 
In order to generates the map, the script should be launched as such: 
main.py --read-vcf-archaic  --myfile *name_vcf_archaic*.vcf.gz 
main.py --read-vcf  --myfile *name_vcf_sample*.vcf.gz --name_dataset *name of the dataset*
main.py --read-outgroup 
main.py --run-hmm --name_dataset *name of the dataset* --chrom *chromosome number* 

Notes:
The name of the outgroup vcf is hardcoded in the script and can be directly downloaded from:
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
The name of the african samples to be extracted from these vcf are also directly written in utils.py
It is possible to change these to use a different outgroup.

This script does not create the output directories, before running it make sure they exist, in particular:

mkdir pickledArchaic
for i in {1..22}; do mkdir pickledArchaic/chr${i}; done
mkdir pickledAfricans
for i in {1..22}; do mkdir pickledAfricans/chr${i}; done
mkdir pickledData
mkdir pickledData/*name dataset*
for i in {1..22}; do mkdir pickledData/*name dataset*/chr${i}; done
mkdir maps
mkdir maps/*name dataset*
for i in {1..22}; do mkdir maps/*name dataset*/chr${i}; done



The output will be one map per individual and per chromosome, it is possible to merge these maps using:
main.py --merge-maps --name_dataset *name of the dataset*

'''
parser.add_argument('--read-vcf-archaic', action='store_true', help='Enable reading archaic VCF file.')
parser.add_argument('--read-vcf', action='store_true', help='Enable reading sample VCF file.')
parser.add_argument('--read-outgroup', action='store_true', help='Enable reading outgroup data.')
parser.add_argument('--run-hmm', action='store_true', help='Run the HMM.')
parser.add_argument('--merge-maps', action='store_true', help='Merge the maps.')


parser.add_argument('--myfile', required=False, help='Path to input VCF file')
parser.add_argument('--chrom', required=False, help='Path to chr')
parser.add_argument('--name_dataset', required=False, help='Path to dataset')



args = parser.parse_args()



if args.myfile is None:
    pass
else:
    myfile = args.myfile

if args.chrom is None:
    pass
else:
    chrom = int(args.chrom)

if args.name_dataset is None:
    pass
else:
    name_dataset = args.name_dataset


args = parser.parse_args()

read_vcf_archaic = args.read_vcf_archaic
read_vcf = args.read_vcf
read_outgroup = args.read_outgroup
run_HMM = args.run_hmm
merge_maps = args.merge_maps

# This script reads the vcf dataset and saves each individual's data into a separate pickle file.
if read_vcf_archaic:
    res = utils.read_vcf(myfile,chrom)
    for i,ind in enumerate(res[1]):
        filehandler = open(f'pickledArchaic/chr{chrom}/{ind}',"wb")
        pickle.dump(res[0][i],filehandler)
        filehandler.close()
elif read_vcf:
    for chrom in range(1,23):
        res = utils.read_vcf(myfile,chrom)
        for i,ind in enumerate(res[1]):
            filehandler = open(f'pickledData/{name_dataset}/chr{chrom}/{ind}',"wb")
            pickle.dump(res[0][i],filehandler)
            filehandler.close()
elif read_outgroup:
    for chrom in range(1,23):
        res = utils.read_vcf_sum(f'ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz',chrom)
        filehandler = open(f'pickledAfricans/chr{chrom}/outgroup',"wb")
        pickle.dump(res,filehandler)
        filehandler.close()

elif run_HMM:
    print(f"Processing chromosome {chrom}")
    outgroup = pickle.load(open(f"pickledAfricans/chr{chrom}/outgroup","rb"))
    neanderthal = pickle.load(open(f"pickledArchaic/chr{chrom}/Vindija33.19","rb"))
    denisova = pickle.load(open(f"pickledArchaic/chr{chrom}/Denisova","rb"))
    altai = pickle.load(open(f"pickledArchaic/chr{chrom}/AltaiNeanderthal","rb"))
    chagyrskaya = pickle.load(open(f"pickledArchaic/chr{chrom}/Chagyrskaya-Phalanx","rb"))
    directory_path = f'pickledData/{name_dataset}/chr{chrom}'
    for filename in os.listdir(directory_path):
        print(f"Processing file {filename} in chromosome {chrom}")
        file_path = os.path.join(directory_path, filename)
        if os.path.isfile(file_path):
            with open(file_path, 'rb') as f:
                sample = pickle.load(f)
                S = HMM.initS(0.04)
                A = HMM.initA(2000,1.2e-8,1000,0.04)
                B = HMM.initB(1.25e-8,1000,2400,30000)
                seq = HMM.create_obs(outgroup,sample,1000)
                res =  HMM.viterbi(seq, S, A,B)
                tractsHMM = HMM.get_HMM_tracts(res)
                if len(tractsHMM) > 1:
                    similarity_neanderthal = []
                    similarity_chagyrskaya = []
                    similarity_altai = []
                    similarity_denisova = []
                    snps_archaic = []
                    snps_neanderthal = []
                    snps_chagyrskaya = []
                    snps_altai = []
                    snps_denisova = []
                    for elt in tractsHMM[1]:
                        similarity_neanderthal.append(utils.get_similarity(elt, sample,neanderthal))
                        similarity_chagyrskaya.append(utils.get_similarity(elt, sample,chagyrskaya))
                        similarity_altai.append(utils.get_similarity(elt, sample,altai))
                        similarity_denisova.append(utils.get_similarity(elt, sample,denisova))
                        same,tot = utils.get_archaic_snps(elt, sample,neanderthal,outgroup)
                        snps_archaic.append(tot)
                        snps_neanderthal.append(same)
                        same,_ = utils.get_archaic_snps(elt, sample,chagyrskaya,outgroup)
                        snps_chagyrskaya.append(same)
                        same,_ = utils.get_archaic_snps(elt, sample,altai,outgroup)
                        snps_altai.append(same)
                        same,_ = utils.get_archaic_snps(elt, sample,denisova,outgroup)
                        snps_denisova.append(same)

                    map_path = Path(f'maps/{name_dataset}/chr{chrom}/{filename}.txt')
                    with map_path.open('w') as f_map:
                        f_map.write('chr\tstart\tend\tsimilarityVindijaNeanderthal\tsimilarityAltaiNeanderthal\tsimilarityChagyrskayaNeanderthal\tsimilarityDenisova\tindividual\tarchaicSnps\tsnpsVindijaNeanderthal\tsnpsAltaiNeanderthal\tsnpsChagyrskayaNeanderthal\tsnpsDenisova\n')
                        HMM.write_tracts(f_map, tractsHMM[1], chrom, filename, similarity_neanderthal, similarity_altai,similarity_chagyrskaya,\
                                         similarity_denisova,snps_archaic,snps_neanderthal, snps_altai,snps_chagyrskaya,\
                                         snps_denisova)
                        f_map.close()
                
                similarity_neanderthal = []
                similarity_chagyrskaya = []
                similarity_altai = []
                similarity_denisova = []
                snps_archaic = []
                snps_neanderthal = []
                snps_chagyrskaya = []
                snps_altai = []
                snps_denisova = []
                for elt in tractsHMM[0]:
                    similarity_neanderthal.append(utils.get_similarity(elt, sample,neanderthal))
                    similarity_chagyrskaya.append(utils.get_similarity(elt, sample,chagyrskaya))
                    similarity_altai.append(utils.get_similarity(elt, sample,altai))
                    similarity_denisova.append(utils.get_similarity(elt, sample,denisova))
                    same,tot = utils.get_archaic_snps(elt, sample,neanderthal,outgroup)
                    snps_archaic.append(tot)
                    snps_neanderthal.append(same)
                    same,_ = utils.get_archaic_snps(elt, sample,chagyrskaya,outgroup)
                    snps_chagyrskaya.append(same)
                    same,_ = utils.get_archaic_snps(elt, sample,altai,outgroup)
                    snps_altai.append(same)
                    same,_ = utils.get_archaic_snps(elt, sample,denisova,outgroup)
                    snps_denisova.append(same)
     

                map_path = Path(f'maps/{name_dataset}/chr{chrom}/{filename}_notArchaic.txt')
                with map_path.open('w') as f_map:
                    f_map.write('chr\tstart\tend\tsimilarityVindijaNeanderthal\tsimilarityAltaiNeanderthal\tsimilarityChagyrskayaNeanderthal\tsimilarityDenisova\tindividual\tarchaicSnps\tsnpsVindijaNeanderthal\tsnpsAltaiNeanderthal\tsnpsChagyrskayaNeanderthal\tsnpsDenisova\n')
                    HMM.write_tracts(f_map, tractsHMM[0], chrom, filename, similarity_neanderthal, similarity_altai,similarity_chagyrskaya,\
                                         similarity_denisova,snps_archaic,snps_neanderthal, snps_altai,snps_chagyrskaya,\
                                         snps_denisova)
                    f_map.close()
elif merge_maps:
    for chrom in range(1,23):
        directory_path = f'maps/{name_dataset}/chr{chrom}'
        for filename in os.listdir(directory_path):
            file_path = os.path.join(directory_path, filename)
            merged_map_path = Path(f'maps/{name_dataset}/merged/{filename}')
            if not merged_map_path.exists():
                with merged_map_path.open('w') as f_merged:
                    f_merged.write('chr\tstart\tend\tsimilarityVindijaNeanderthal\tsimilarityAltaiNeanderthal\tsimilarityChagyrskayaNeanderthal\tsimilarityDenisova\tindividual\tarchaicSnps\tsnpsVindijaNeanderthal\tsnpsAltaiNeanderthal\tsnpsChagyrskayaNeanderthal\tsnpsDenisova\n')
            with open(file_path,'r') as f_map:
                next(f_map)  # Skip header
                with merged_map_path.open('a') as f_merged:
                    for line in f_map:
                        f_merged.write(line)
