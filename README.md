# arcAncientImputed
## Data

The introgression map is in `data/maps/IntrogressionMaps.csv`, it contains the archaic segments detected for all 20 ancient samples (described in `data/other/InfoSamples.txt`), for each segment we provide:


| Column       | Description                                                                                           |
| ------------ | ----------------------------------------------------------------------------------------------------- |
| `Individual` | Identifier of the ancient individual in which the fragment is detected                                |
| `Coverage`   | Dataset type used to infer the fragment (e.g. `Original`, `ImputedOC`, `2X`, `1X`, `0.5X`, `0.0625X`) |
| `chr`        | Chromosome name (e.g. `"chr1"`)                                                                       |
| `Chrom`      | Chromosome number (integer)                                                                           |
| `start`      | Start genomic position of the fragment (bp, hg19)                                                     |
| `end`        | End genomic position of the fragment (bp, hg19)                                                       |
| `Length`     | Physical length of the fragment in base pairs                                                         |
| `Length.cM`  | Genetic length of the fragment in centiMorgans (computed using a recombination map)                   |
| `Segment_Type`     | Classification of the fragment:<br>• **Original**: inferred directly from original data<br>• **Archaic Shared**: overlaps archaic fragments detected in original data<br>• **Archaic New**: does not overlap archaic fragments detected in original data |
| `Filter`           | Post-processing filter status (e.g. `Length ≥ 40 kb`, `Similarity ≥ 0.99`)                                                                                                                                                                               |
| `StrictMask1000GP` | Whether the fragment lies within the strict 1000 Genomes accessibility mask                                                                                                                                                      
| `SimilarityOriginal_Neanderthal` | Similarity to Vindija Neanderthal genomes             |
| `SimilarityOriginal_Denisovan`   | Similarity to Denisovan genome                        |
| `SimilarityOriginal_NEAminusDEN` | Neanderthal minus Denisovan similarity                |
| `SimilarityOriginal_to_archaic`  | Maximum similarity to either Neanderthal or Denisovan |
| `DistanceOriginal_Neanderthal`   | Genetic distance to Vindija Neanderthal genomes       |
| `DistanceOriginal_Denisovan`     | Genetic distance to Denisovan genome                  |
| `DistanceOriginal_NEAminusDEN`   | Neanderthal minus Denisovan distance                  |
| `DistanceOriginal_to_archaic`    | Minimum distance to either archaic genome             |
| `SimilarityImp_Neanderthal` | Similarity to Vindija Neanderthal genomes             |
| `SimilarityImp_Denisovan`   | Similarity to Denisovan genome                        |
| `SimilarityImp_NEAminusDEN` | Neanderthal minus Denisovan similarity                |
| `SimilarityImp_to_archaic`  | Maximum similarity to either Neanderthal or Denisovan |
| `DistanceImp_Neanderthal`   | Genetic distance to Vindija Neanderthal genomes       |
| `DistanceImp_Denisovan`     | Genetic distance to Denisovan genome                  |
| `DistanceImp_NEAminusDEN`   | Neanderthal minus Denisovan distance                  |
| `DistanceImp_to_archaic`    | Minimum distance to either archaic genome             |
| `hmmPositionsImp`                             | Number of HMM-informative (African-divergent) SNPs |
| `PrivateSNPsImp_VindijaNeanderthal`           | Shared SNPs with Vindija Neanderthal               |
| `PrivateSNPsImputed_ChagyrskayaNeanderthal`   | Shared SNPs with Chagyrskaya Neanderthal           |
| `PrivateSNPsImputed_AltaiNeanderthal`         | Shared SNPs with Altai Neanderthal                 |
| `PrivateSNPsImp_Denisovan`                    | Shared SNPs with Denisovan genome                  |
| `SimilarityPrivateImp_VindijaNeanderthal`     | Ratio of shared SNPs to `hmmPositionsImp`          |
| `SimilarityPrivateImp_ChagyrskayaNeanderthal` | Ratio of shared SNPs to `hmmPositionsImp`          |
| `SimilarityPrivateImp_AltaiNeanderthal`       | Ratio of shared SNPs to `hmmPositionsImp`          |
| `SimilarityPrivateImp_Denisovan`              | Ratio of shared SNPs to `hmmPositionsImp`          |


## Scripts

### Computing the dstat, using GP

The script to compute the D-statistics using GP is available in `src/AlleleFrequencies`. To use it, we run:

```
./dstat.py -p1 <vcf containing the samples in P1>.vcf.gz -p2 <vcf containing the samples in P2>.vcf.gz -p3 <vcf containing the samples in P3>.vcf.gz -p4 <vcf containing the samples in P4>.vcf.gz -out <path for the output>
```

Only the population P2 is considered as imputed, for P1, P3 and P4, allele frequency is used. The script is based on a previous implementation of [RIPTA](https://github.com/David-Peede/RIPTA).

Test data with three ancient individuals and chromosome 2 is available here:
https://drive.google.com/drive/folders/1GQlXtcmGsybFq6E94M4BF3eesvqhPLdf?usp=sharing

### Local archaic ancestry inference

To identify the introgression segments in each individual and coverage we used an in-house version of hmmix. The script main.py can be found in the `src/Segment` folder in this repository.

As described in the script, we first convert all the vcf files in a pickle file format, we then run the hmm function to determine the introgressed segments and compute similarities to reference archaic genomes.

The script is not intended to be run as a standalone program; it is designed as part of a larger analysis pipeline where the pickle format allows to speedup the process of multiple experiments.

In order to generate the map, the script should be launched as such:
```
main.py --read-vcf-archaic  --myfile *name_vcf_archaic*.vcf.gz
main.py --read-vcf  --myfile *name_vcf_sample*.vcf.gz --name_dataset *name of the dataset*
main.py --read-outgroup
main.py --run-hmm --name_dataset *name of the dataset* --chrom *chromosome number*
```

Notes:

The name of the outgroup vcf is hardcoded in the script and can be directly downloaded from: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

The name of the african samples to be extracted from these vcf are also directly written in utils.py

It is possible to change these to use a different outgroup.

This script does not create the output directories, before running it make sure they exist, in particular: 

```
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
```

The output will be one map per individual and per chromosome, it is possible to merge these maps using:
```
main.py --merge-maps --name_dataset *name of the dataset*
```
