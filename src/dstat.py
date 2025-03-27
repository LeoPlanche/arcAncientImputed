import argparse
import gzip

# Intialize the command line options.
parser = argparse.ArgumentParser()

parser.add_argument(
    '-p1', '--p1_population', required=True,
    type=str, action='store',
    help='VCF for population P1, computed using allele frequency',
)
parser.add_argument(
    '-p2', '--p2_population', required=True,
    type=str, action='store',
    help='VCF for imputed samples (P2), each sample is computed separately, using GP',
)
parser.add_argument(
    '-p3', '--p3_population', required=True,
    type=str, action='store',
    help='VCF for population P3 (the source population), computed using allele frequency.',
)
parser.add_argument(
    '-p4', '--p4_population', required=True,
    type=str, action='store',
    help='VCF for population P4 (the outgroup population used for polarization), computed using allele frequency.',
)
parser.add_argument(
    '-out', '--out_path', required=True,
    type=str, action='store',
    help='Path for results.',
)

args = parser.parse_args()

if args.p1_population.endswith('.gz'):    
    f1 = gzip.open(args.p1_population, 'rt')
else:
    f1 = open(args.p1_population,'r')

if args.p2_population.endswith('.gz'):    
    f2 = gzip.open(args.p2_population, 'rt')
else:
    f2 = open(args.p2_population,'r')

if args.p3_population.endswith('.gz'):    
    f3 = gzip.open(args.p3_population, 'rt')
else:
    f3 = open(args.p3_population,'r')

if args.p4_population.endswith('.gz'):    
    f4 = gzip.open(args.p4_population, 'rt')
else:
    f4 = open(args.p4_population,'r')

names=[]
num=[]
dem=[]

line1=next(f1)
while line1[0] == '#':
    line1=next(f1)

line2=next(f2)
while line2[1] == '#':
    line2=next(f2)
sline2 = line2.split()
for elt in sline2[9:]:
    names.append(elt)
    num.append(0)
    dem.append(0)
line2=next(f2)
sline2 = line2.split()

line3=next(f3)
while line3[0] == '#':
    line3=next(f3)
sline3 = line3.split()
line4=next(f4)
while line4[0] == '#':
    line4=next(f4)
sline4 = line4.split()

for line1 in f1:
    sline1 = line1.split()
    chr=int(sline1[0])     
    pos = int(sline1[1]) 
    if chr >= int(sline2[0]) and chr >= int(sline3[0]) and chr >= int(sline4[0]):
        try:   
            while (int(sline2[1]) < pos and int(sline2[0]) ==  chr) or int(sline2[0]) <  chr :
                sline2 = next(f2).split() 
                
            while (int(sline3[1]) < pos and int(sline3[0]) ==  chr) or int(sline3[0]) <  chr :
                sline3 = next(f3).split() 

            while (int(sline4[1]) < pos and int(sline4[0]) ==  chr) or int(sline4[0]) <  chr :
                sline4 = next(f4).split() 
            
            # We are at the same position in all the files
            if pos==int(sline2[1]) and pos==int(sline3[1]) and pos == int(sline4[1]):
                
                #Define ancestral state using P4 
                if(sline4[9][0]=='0' and (sline4[9][2]=='0' or sline4[9][2]=='.')):
                    ancestral = sline4[3]
                    derived = sline4[4]
                elif(sline4[9][0]=='1' and (sline4[9][2]=='1' or sline4[9][2]=='.')): 
                    ancestral = sline4[4]
                    derived = sline4[3]
                else:
                    continue
            
                #Compute allele frequency for P1 
                n0=0
                n1=0
                p1=0
                for elt in sline1[9:]:
                    if(elt[0]=='0'):
                        n0+=1
                    elif(elt[0]=='1'): 
                        n1+=1
                    if(elt[2]=='0' ):
                        n0+=1
                    elif(elt[2]=='1'): 
                        n1+=1
                if n1+n0 > 0:
                    p1=n1/(n1+n0)
                else:
                    continue

                if ancestral==sline1[4] and derived==sline1[3]:
                    p1=1-p1
                elif ancestral!=sline1[3]:
                    continue

                #Compute allele frequency for P3 
                n0=0
                n1=0
                p3=0
                for elt in sline3[9:]:
                    if(elt[0]=='0'):
                        n0+=1
                    elif(elt[0]=='1'): 
                        n1+=1
                    if(elt[2]=='0' ):
                        n0+=1
                    elif(elt[2]=='1'): 
                        n1+=1
                if n1+n0 > 0:
                    p3=n1/(n1+n0)
                else:
                    continue

                if ancestral==sline3[4] and derived==sline3[3]:
                    p3=1-p3
                elif ancestral!=sline3[3]:
                    continue              
                
                #Compute allele frequency using GP, for all individuals in P3       
                for idx,elt in enumerate(sline2[9:]):
                    GP=(elt.split(sep=':')[2]).split(sep=',')               
                    p2=float(GP[2])+float(GP[1])/2
                    if ancestral==sline2[4] and derived==sline2[3]:
                        p2=1-p2
                    elif ancestral!=sline2[3]:
                        continue

                    num[idx] += (p1-p2)*p3
                    dem[idx] += (p1+p2-p1*p2)*p3                               
        except StopIteration:
            # End of file reached
            break
f1.close()
f2.close()
f3.close()
f4.close()

fou = open(args.out_path,'w')

for idx,elt in enumerate(names):
    fou.write(f'{elt}\t{num[idx] / dem[idx]} \n')
fou.close()