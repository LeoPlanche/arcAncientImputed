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

The script is not intended to be run as a standalone program, it is designed as part of a larger analysis pipeline where the pickle format allows to speedup the process of multiple experiments.

In order to generate the maps, the script should be launched as such: 
```
main.py --read-vcf-archaic  --file *name_vcf_archaic*.vcf.gz 
main.py --read-vcf  --file *name_vcf_sample*.vcf.gz --name-dataset *name of the dataset*
main.py --read-outgroup  --directory *name of the directory containing the outgroup vcfs*
main.py --run-hmm --name-dataset *name of the dataset* --chrom *chromosome number* 
```

#### Notes:

#### Archaic vcf  processing
```--read-vcf-archaic``` assumes that the vcf given as input contains all chromosomes.

If multiple archaic individuals are in the same vcf, all the statistics (similarity, archaic snps) will be computed for all archaic individuals.

It is possible to run the command multiple time if different archaic individuals are present in different vcfs.

#### Target data vcf processing
```--read-vcf``` assumes that the vcf given as input contains all chromosomes.

If multiple individuals are in the same vcf, all the archaic segments and statistics will be computed for all individuals.

It is possible to run the command multiple time if different individuals are present in different vcfs.

#### Outgroup vcf processing
The name of the outgroup vcfs is hardcoded in the script and the files (one per chromosome) can be directly downloaded from:

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

It is assumed that the directory given as input for ```--read-outgroup``` contains those files.

The name of the african samples to be extracted from these vcf are also directly written in utils.py

It is possible to change these to use a different outgroup.

#### Output
The output directory will be created at ```../out``` and the maps will be in ```../out/maps/*name of the dataset*```

The output will be one map per individual and per chromosome, it is possible to merge these maps using:

```
main.py --merge-maps --name-dataset *name of the dataset*
```

Then the maps will be in ```../out/maps/*name of the dataset*/merged```

The field ```--name-dataset``` is simply required to name output directories.

Only .vcf.gz format is supported (not .bcf), it is assumed that the chromosome field is an integer (e.g. ```1``` and not ```chr1```)

### Imputation Accuracy

To compute the imputation accuracy we use the script imputationAccuracy.py that can be found in the `src/Concordance` folder in this repository.

As described in the script, we first convert all the vcf files in a pickle file format, we then run the analysis to determine imputation accuracy.

If ```pickledArchaic``` and ```pickledOutgroup``` folders (and files) have not been created through the LAI script, first run:

```
imputationAccuracy.py --read-vcf-archaic  --file *name_vcf_archaic*.vcf.gz 
imputationAccuracy.py --read-outgroup  --directory *name of the directory containing the outgroup vcfs*
```

This script assumes that we have the same data for at least two coverages, the original one and the (imputed) downsampled.
For each coverage, run:

```
imputationAccuracy.py --read-vcf  --file *name_vcf_sample*.vcf.gz --coverage *name of the dataset (coverage)*
```

Finally, run:

```
imputationAccuracy.py --imputation-quality-archaic-snps --coverage *name of the dataset (coverage)* --chrom *chromosome*
```

    
#### Notes:

#### Archaic and outgroup vcf processing

  If the LAI has been run through the script described above, the preprocessing of the archaic and outgroup genomes should have already been done.

  ```--read-vcf-archaic``` assumes that the vcf given as input contains all chromosomes.
  
  If multiple archaic individuals are in the same vcf, all the archaic individuals will be separately used to define archaic SNPs
  
  It is possible to run the command multiple time if different archaic individuals are present in different vcfs

  The outgroup is used to define archaic SNPs
  
  The name of the outgroup vcf is hardcoded in the script and can be directly downloaded from:
  
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
  
  It is assumed that the directory given as input for ```--read-outgroup``` contains those files
  
  The name of the african samples to be extracted from these vcf are also directly written in utils.py
  
  It is possible to change these to use a different outgroup. 

#### Target data vcf processing

  This script assumes that we have the same data for at least two coverages, the original one and the (imputed) downsampled.

  ```--read-vcf``` assumes that the vcf given as input contains all chromosomes.
  
  If multiple individuals are in the same vcf, the imputation accuracy will be computed for all individuals
  
  It is possible to run the command multiple time for the same coverage if different individuals are present in different vcfs
  
  <!> For the non imputed dataset (which will serve as a reference to compute imputation accuracy), the coverage should be exactly named "original"
  
  <!> For each coverage the same individual should have the same name, if not either rename the individuals in the vcfs or in the output files in ```../out/pickledData/*name of the dataset (coverage)*```

#### Imputation accuracy 

  The output files will be in ```../out/concordance/{coverage}/archaic``` and ```../out/concordance/{coverage}/nonArchaic```
  
  The first folder contains imputation accuracy for archaic SNPs, the second for non archaic SNPs.
  
  The accuracy is given in a separate file for each chromosome. The files are tabs separared with two collumns, each row corresponds to an allele with as first value its allele frequency in the 1000 genome project and as second value its imputation accuracy when comparing the imputed data with the original data.


### Demo version of the D-statistic script

To test the new D-statistic approach, we prepared a demo version that can be run on a small dataset. The demo includes data from chromosome 2 for three individuals in our dataset: 2H10 (France, 5165 BCE), Mota (Ethiopia, 4470 BCE), and Ust’Ishim (Russia, 44,366 BCE).
The script dstat.py is written in Python 3 and requires only two standard libraries: argparse and gzip. It is a standalone script and does not require installation.
This demo version was executed from the terminal on Ubuntu 22.04.5 LTS using a laptop with 32 GB of RAM.

Download
A compressed folder containing all the files required to run the demo version can be downloaded from Zenodo (DOI: 10.5281/zenodo.16365479).

Running the demo
To run the demo version, use the following commands.

For 2H10
```
python3 dstat.py \
-p1 1000GP_YRI_ESN_MSL_ChimpImputedSNPs_Chr2.vcf.gz \
-p2 2H10_Raw_1X_Chr2.vcf.gz \
-p3 Vindija_ChimpImputedSNPs_Chr2.vcf.gz \
-p4 hg19_panTro6.archaicSNPs.annotated_Chr2.vcf.gz \
-out 2H10_Chr2_Vindija.txt
```

For Ust’Ishim
```
python3 dstat.py \
-p1 1000GP_YRI_ESN_MSL_ChimpImputedSNPs_Chr2.vcf.gz \
-p2 Ust_Ishim_Raw_1X_Chr2.vcf.gz \
-p3 Vindija_ChimpImputedSNPs_Chr2.vcf.gz \
-p4 hg19_panTro6.archaicSNPs.annotated_Chr2.vcf.gz \
-out Ust_Ishim_Chr2_Vindija.txt
```

For Mota
```
python3 dstat.py \
-p1 1000GP_YRI_ESN_MSL_ChimpImputedSNPs_Chr2.vcf.gz \
-p2 Mota_Raw_1X_Chr2.vcf.gz \
-p3 Vindija_ChimpImputedSNPs_Chr2.vcf.gz \
-p4 hg19_panTro6.archaicSNPs.annotated_Chr2.vcf.gz \
-out Mota_Chr2_Vindija.txt
```

Runtime
Each run takes approximately 8 minutes on the system described above.

Output
The output is a text file containing two columns:
The name of the individual in the VCF file
The calculated D-statistic value
Example results are provided in the Results folder.
