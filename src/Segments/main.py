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
main.py --read-vcf  --file *name_vcf_sample*.vcf.gz --name_dataset *name of the dataset*
main.py --read-outgroup  --directory *name of the directory containing the outgroup vcfs*
main.py --run-hmm --name-dataset *name of the dataset* --chrom *chromosome number* 

Notes:
--read-vcf-archaic assumes that the vcf given as input contains all chromosomes.
If multiple archaic individuals are in the same vcf, all the statistics (similarity, archaic snps) will be computed for all archaic individuals
It is possible to run the command multiple time if different archaic individuals are present in different vcfs

--read-vcf assumes that the vcf given as input contains all chromosomes.
If multiple individuals are in the same vcf, all the archaic segments and statistics will be computed for all individuals
It is possible to run the command multiple time if different individuals are present in different vcfs

The name of the outgroup vcf is hardcoded in the script and can be directly downloaded from:
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
It is assumed that the directory given as input for --read-outgroup contains those files
The name of the african samples to be extracted from these vcf are also directly written in utils.py
It is possible to change these to use a different outgroup.

The output directory will be created at ../out
The maps will be in ../out/maps/*name of the dataset*

The output will be one map per individual and per chromosome, it is possible to merge these maps using:
main.py --merge-maps --name-dataset *name of the dataset*
Then the maps will be in ../out/maps/*name of the dataset*/merged

The field --name-dataset is simply required to name output directories

Only .vcf.gz format is supported (not .bcf), it is assumed that the chromosome field is an integer (e.g. 1 and not chr1)
'''
parser.add_argument('--read-vcf-archaic', action='store_true', help='Enable reading archaic VCF file.')
parser.add_argument('--read-vcf', action='store_true', help='Enable reading sample VCF file.')
parser.add_argument('--read-outgroup', action='store_true', help='Enable reading outgroup data.')
parser.add_argument('--run-hmm', action='store_true', help='Run the HMM.')
parser.add_argument('--merge-maps', action='store_true', help='Merge the maps.')


parser.add_argument('--file', required=False, help='Path to input VCF file')
parser.add_argument('--chrom', required=False, help='Chromosome number')
parser.add_argument('--name-dataset', required=False, help='Name of the dataset (used to store output files)')
parser.add_argument('--directory', required=False, help='Directory containing input VCF files')


args = parser.parse_args()


if args.file is None:
    pass
else:
    file = args.file

if args.chrom is None:
    pass
else:
    chrom = int(args.chrom)

if args.name_dataset is None:
    pass
else:
    name_dataset = args.name_dataset

if args.directory is None:
    pass
else:
    directory = args.directory


read_vcf_archaic = args.read_vcf_archaic
read_vcf = args.read_vcf
read_outgroup = args.read_outgroup
run_HMM = args.run_hmm
merge_maps = args.merge_maps

# This script reads the vcf dataset and saves each individual's data into a separate pickle file.
if read_vcf_archaic:
    for chrom in range(1,23):
        res = utils.read_vcf(file,chrom)
        for i,ind in enumerate(res[1]):
            path = Path(f'../out/pickledArchaic/chr{chrom}')
            path.mkdir(parents=True, exist_ok=True)
            filehandler = open(f'../out/pickledArchaic/chr{chrom}/{ind}',"wb")
            pickle.dump(res[0][i],filehandler)
            filehandler.close()
elif read_vcf:
    for chrom in range(1,23):
        res = utils.read_vcf(file,chrom)
        for i,ind in enumerate(res[1]):
            path = Path(f'../out/pickledData/{name_dataset}/chr{chrom}')
            path.mkdir(parents=True, exist_ok=True)
            filehandler = open(f'../out/pickledData/{name_dataset}/chr{chrom}/{ind}',"wb")
            pickle.dump(res[0][i],filehandler)
            filehandler.close()
elif read_outgroup:
    for chrom in range(1,23):
        res = utils.read_vcf_sum(f'{directory}/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz',chrom)
        path = Path(f'../out/pickledAfricans/chr{chrom}')
        path.mkdir(parents=True, exist_ok=True)
        filehandler = open(f'../out/pickledAfricans/chr{chrom}/outgroup',"wb")
        pickle.dump(res,filehandler)
        filehandler.close()

elif run_HMM:
    print(f"Processing chromosome {chrom}")
    outgroup = pickle.load(open(f"../out/pickledAfricans/chr{chrom}/outgroup","rb"))
    archaic_individuals = {}
    archaic_individuals_names = []
    directory_path = f"../out/pickledArchaic/chr{chrom}"
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)
        f = open(file_path,"rb")
        archaic_individuals[filename] = pickle.load(f)
        archaic_individuals_names.append(filename)
    directory_path = f'../out/pickledData/{name_dataset}/chr{chrom}'
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

                
                similarity = {}
                snps = {}
                snps_archaic = []
                for _ , ind in enumerate(archaic_individuals_names):
                    similarity[ind] = []
                    snps[ind] = []
                for elt in tractsHMM[0]:
                    tot = 0
                    for _, ind in enumerate(archaic_individuals_names):
                        similarity[ind].append(utils.get_similarity(elt, sample,archaic_individuals[ind]))
                        same , tot = utils.get_archaic_snps(elt, sample,archaic_individuals[ind],outgroup)
                        snps[ind].append(same)
                    snps_archaic.append(tot)

                path = Path(f'../out/maps/{name_dataset}/chr{chrom}')
                path.mkdir(parents=True, exist_ok=True)
                map_path = Path(f'../out/maps/{name_dataset}/chr{chrom}/{filename}_notArchaic.txt')
                with map_path.open('w') as f_map:
                        header = 'chr\tstart\tend'
                        for _,ind in enumerate(archaic_individuals_names):
                            header += f'\tsimilarity{ind}'
                        header += '\tindividual\tarchaicSnps'
                        for _,ind in enumerate(archaic_individuals_names):
                            header += f'\tsnps{ind}'
                        header += '\n'
                        f_map.write(header)
                        HMM.write_tracts(f_map, tractsHMM[1], chrom, archaic_individuals_names, similarity,snps_archaic,snps)
                        f_map.close()

                if len(tractsHMM) > 1:
                    similarity = {}
                    snps = {}
                    snps_archaic = []
                    for _,ind in enumerate(archaic_individuals_names):
                        similarity[ind] = []
                        snps[ind] = []
                    for elt in tractsHMM[1]:
                        tot = 0
                        for _,ind in enumerate(archaic_individuals_names):
                            similarity[ind].append(utils.get_similarity(elt, sample,archaic_individuals[ind]))
                            same , tot = utils.get_archaic_snps(elt, sample,archaic_individuals[ind],outgroup)
                            snps[ind].append(same)
                        snps_archaic.append(tot)
                    
                    map_path = Path(f'../out/maps/{name_dataset}/chr{chrom}/{filename}.txt')
                    with map_path.open('w') as f_map:
                        header = 'chr\tstart\tend'
                        for _,ind in enumerate(archaic_individuals_names):
                            header += f'\tsimilarity{ind}'
                        header += '\tindividual\tarchaicSnps'
                        for _,ind in enumerate(archaic_individuals_names):
                            header += f'\tsnps{ind}'
                        header += '\n'
                        f_map.write(header)
                        HMM.write_tracts(f_map, tractsHMM[1], chrom, archaic_individuals_names, similarity,snps_archaic,snps)
                        f_map.close()
                
                
elif merge_maps:
    path = Path(f'../out/maps/{name_dataset}/merged/')
    path.mkdir(parents=True, exist_ok=True)
    for chrom in range(1,23):
        directory_path = f'../out/maps/{name_dataset}/chr{chrom}'
        for filename in os.listdir(directory_path):
            file_path = os.path.join(directory_path, filename)
            merged_map_path = Path(f'../out/maps/{name_dataset}/merged/{filename}')
            with open(file_path,'r') as f_map:
                header = next(f_map)
                if not merged_map_path.exists():
                    with merged_map_path.open('w') as f_merged:
                        f_merged.write(header)
                with merged_map_path.open('a') as f_merged:
                    for line in f_map:
                        f_merged.write(line)
