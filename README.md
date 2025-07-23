# arcAncientImputed
## Computing the dstat, using GP

The script to compute the D-statistics using GP is available in `src`. To use it, we run:

```
./dstat.py -p1 <vcf containing the samples in P1>.vcf.gz -p2 <vcf containing the samples in P2>.vcf.gz -p3 <vcf containing the samples in P3>.vcf.gz -p4 <vcf containing the samples in P4>.vcf.gz -out <path for the output>
```

Only the population P2 is considered as imputed, for P1, P3 and P4, allele frequency is used. The script is based on a previous implementation of [RIPTA](https://github.com/David-Peede/RIPTA).

Test data with three ancient individuals and chromosome 2 is available here:
https://drive.google.com/drive/folders/1GQlXtcmGsybFq6E94M4BF3eesvqhPLdf?usp=sharing
