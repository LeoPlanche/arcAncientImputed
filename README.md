# arcAncientImputed
## Computing the introgression maps

We compute the introgression maps using [DAIseg](https://github.com/LeoPlanche/DAIseg), which in the task of detecting a single archaic origin using a single outgroup is equivalent to [hmmix](https://github.com/LauritsSkov/Introgression-detection), on which the method is based.

DAIseg requires a demography and individual file, they are given in `src`. 

After running DAIseg (or hmmix), in order to compute the similarity to the archaic genomes, we run:

```
./main.py similarity -decode <decode folder obtained after running DAIseg> -p1 <vcf containing the ancient samples>.vcf.gz -p2 <vcf(s) containing the archaic genome>.vcf.gz -n2 <name of the archaic sample>
```

To compute the distance, using an outgroup:

```
./main.py similarity -decode <decode folder obtained after running DAIseg> -p1 <vcf containing the ancient samples>.vcf.gz -p2 <vcf(s) containing the archaic genome>.vcf.gz -n2 <name of the archaic sample> -p3 <vcf(s) containing the outgroup genome>.vcf.gz -n3 <name(s) of the outgroup sample(s)>.json
```

Then, to build the introgression map:

```
./main.py create_map -decode <decode folder obtained after running DAIseg and computing the similarities> -out <path for the introgression map>
```
