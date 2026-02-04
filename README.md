# arcAncientImputed
## Data

The introgression map is in `data/maps/IntrogressionMaps.csv`, it contains the archaic segments detected for all 20 ancient samples (described in `data/other/InfoSamples.txt`), for each segment we provide:

Column name  Explanation
Individual  Identifier of the ancient individual in which the fragment is detected.
Coverage  Dataset type used to infer the fragment (e.g. Original, ImputedOC, 2X, 1X, 0.5X, 0.0625X).
chr  Chromosome name (e.g. "chr1").
Chrom  Chromosome number as integer.
start  Start genomic position of the fragment (bp, hg19).
end  End genomic position of the fragment (bp, hg19).
Length  Physical length of the fragment in base pairs.
Length.cM.  Genetic length of the fragment in centiMorgans, computed using a recombination map.
Segment_Type  Classification of the fragment: Original: fragment inferred directly from the original data; Archaic Shared: archaic fragment overlapping introgressed regions detected with original data; Archaic New: archaic fragment not overlapping introgressed regions detected with original data.
Filter  Post-processing filter status (e.g. Length >=40Kb; Similarity to the closest archaic genome >=0.99).
StrictMask1000GP  Indicates whether the fragment lies within the strict 1000 Genomes accessibility mask.
SimilarityOriginal_Neanderthal  Similarity to Vindija Neanderthal genomes using original (high-coverage) data.
SimilarityOriginal_Denisovan  Similarity to Denisova genomes using original (high-coverage) data.
SimilarityOriginal_NEAminusDEN  Difference between Neanderthal and Denisovan similarity using original (high-coverage) data.
SimilarityOriginal_to_archaic  Maximum similarity to either Neanderthal or Denisovan using original (high-coverage) data.
DistanceOriginal_Neanderthal  Genetic distance to Vindija Neanderthal genomes using original (high-coverage) data.
DistanceOriginal_Denisovan  Genetic distance to Denisova genomes using original (high-coverage) data.
DistanceOriginal_NEAminusDEN  Difference between Neanderthal and Denisovan distance using original (high-coverage) data.
DistanceOriginal_to_archaic  Minimum distance to either archaic genome using original (high-coverage) data.
SimilarityImp_Neanderthal  Similarity to Vindija Neanderthal genomes using imputed data.
SimilarityImp_Denisovan  Similarity to Denisova genomes using imputed data.
SimilarityImp_NEAminusDEN  Difference between Neanderthal and Denisovan similarity using imputed data.
SimilarityImp_to_archaic  Maximum similarity to either Neanderthal or Denisovan using imputed data.
DistanceImp_Neanderthal  Genetic distance to Vindija Neanderthal genomes using imputed data.
DistanceImp_Denisovan  Genetic distance to Denisova genomes using imputed data.
DistanceImp_NEAminusDEN  Difference between Neanderthal and Denisovan distance using imputed data.
DistanceImp_to_archaic  Minimum distance to either archaic genome using imputed data.
hmmPositionsOriginal  Number of HMM-informative SNPs (African-divergent) present in the original (high-coverage) data.
PrivateSNPsOriginal_VindijaNeanderthal  Number of hmmPositionsOriginal shared with the Vindija Neanderthal genome using original (high-coverage) data.
PrivateSNPsOriginal_ChagyrskayaNeanderthal  Number of hmmPositionsOriginal shared with the Chagyrskaya Neanderthal genome using original (high-coverage) data.
PrivateSNPsOriginal_AltaiNeanderthal  Number of hmmPositionsOriginal shared with the Altai Neanderthal genome using original (high-coverage) data.
PrivateSNPsOriginal_Denisovan  Number of hmmPositionsOriginal shared with the Denisovan genome using original (high-coverage) data.
SimilarityPrivateOriginal_VindijaNeanderthal   Ratio of PrivateSNPsOriginal_VindijaNeanderthal to hmmPositionsOriginal
SimilarityPrivateOriginal_ChagyrskayaNeanderthal  Ratio of PrivateSNPsOriginal_ChagyrskayaNeanderthal to hmmPositionsOriginal
SimilarityPrivateOriginal_AltaiNeanderthal  Ratio of PrivateSNPsOriginal_AltaiNeanderthal to hmmPositionsOriginal
SimilarityPrivateOriginal_Denisovan  Ratio of PrivateSNPsOriginal_Denisovan to hmmPositionsOriginalhmm
PositionsImp  Number of HMM-informative SNPs (African-divergent) present in the imputed data.
PrivateSNPsImp_VindijaNeanderthal  Number of hmmPositions shared with the Vindija Neanderthal genome using imputed data.
PrivateSNPsImputed_ChagyrskayaNeanderthal  Number of hmmPositions shared with the Chagyrskaya Neanderthal genome using imputed data.
PrivateSNPsImputed_AltaiNeanderthal  Number of hmmPositions shared with the Altai Neanderthal genome using imputed data.
PrivateSNPsImp_Denisovan  Number of hmmPositions shared with the Denisovan genome using imputed data.
SimilarityPrivateImp_VindijaNeanderthal  Ratio of PrivateSNPsImp_VindijaNeanderthal to hmmPositions
SimilarityPrivateImp_ChagyrskayaNeanderthal  Ratio of PrivateSNPsImp_ChagyrskayaNeanderthal to hmmPositions
SimilarityPrivateImp_AltaiNeanderthal  Ratio of PrivateSNPsImp_AltaiNeanderthal to hmmPositions
SimilarityPrivateImp_Denisova  Ratio of PrivateSNPsImp_Denisovan to hmmPositions

## Computing the dstat, using GP

The script to compute the D-statistics using GP is available in `src/AlleleFrequencies`. To use it, we run:

```
./dstat.py -p1 <vcf containing the samples in P1>.vcf.gz -p2 <vcf containing the samples in P2>.vcf.gz -p3 <vcf containing the samples in P3>.vcf.gz -p4 <vcf containing the samples in P4>.vcf.gz -out <path for the output>
```

Only the population P2 is considered as imputed, for P1, P3 and P4, allele frequency is used. The script is based on a previous implementation of [RIPTA](https://github.com/David-Peede/RIPTA).

Test data with three ancient individuals and chromosome 2 is available here:
https://drive.google.com/drive/folders/1GQlXtcmGsybFq6E94M4BF3eesvqhPLdf?usp=sharing
