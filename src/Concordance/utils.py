import gzip 

outgroup = ["HG02922", "HG02923", "HG02938", "HG02941", "HG02943", "HG02944", "HG02946", "HG02947", "HG02952", "HG02953", "HG02968",\
 "HG02970", "HG02971", "HG02973", "HG02974", "HG02976", "HG02977", "HG02979", "HG02981", "HG03052", "HG03054", "HG03055", "HG03057", \
 "HG03058", "HG03060", "HG03061", "HG03063", "HG03064", "HG03066", "HG03069", "HG03072", "HG03073", "HG03074", "HG03077", "HG03078", \
 "HG03079", "HG03081", "HG03082", "HG03084", "HG03085", "HG03086", "HG03088", "HG03091", "HG03095", "HG03096", "HG03097", "HG03099", \
 "HG03100", "HG03103", "HG03105", "HG03108", "HG03109", "HG03111", "HG03112", "HG03114", "HG03115", "HG03117", "HG03118", "HG03120", \
 "HG03121", "HG03123", "HG03124", "HG03126", "HG03127", "HG03129", "HG03130", "HG03132", "HG03133", "HG03135", "HG03136", "HG03139", \
 "HG03157", "HG03159", "HG03160", "HG03162", "HG03163", "HG03166", "HG03168", "HG03169", "HG03172", "HG03175", "HG03189", "HG03190", \
 "HG03193", "HG03195", "HG03196", "HG03198", "HG03199", "HG03202", "HG03209", "HG03212", "HG03224", "HG03225", "HG03265", "HG03267", \
 "HG03268", "HG03270", "HG03271", "HG03279", "HG03280", "HG03291", "HG03294", "HG03295", "HG03297", "HG03298", "HG03300", "HG03301", \
 "HG03303", "HG03304", "HG03311", "HG03313", "HG03342", "HG03343", "HG03351", "HG03352", "HG03354", "HG03363", "HG03366", "HG03367", \
 "HG03369", "HG03370", "HG03372", "HG03376", "HG03378", "HG03380", "HG03382", "HG03385", "HG03388", "HG03391", "HG03394", "HG03397", \
 "HG03401", "HG03410", "HG03419", "HG03428", "HG03432", "HG03433", "HG03436", "HG03437", "HG03439", "HG03442", "HG03445", "HG03446", \
 "HG03449", "HG03451", "HG03452", "HG03455", "HG03457", "HG03458", "HG03460", "HG03461", "HG03464", "HG03469", "HG03470", "HG03472", \
 "HG03473", "HG03476", "HG03478", "HG03479", "HG03484", "HG03485", "HG03499", "HG03511", "HG03514", "HG03515", "HG03517", "HG03518", \
 "HG03520", "HG03521", "HG03547", "HG03548", "HG03556", "HG03557", "HG03558", "HG03559", "HG03563", "HG03565", "HG03567", "HG03571", \
 "HG03572", "HG03575", "HG03577", "HG03578", "HG03583", "NA18486", "NA18488", "NA18489", "NA18498", "NA18499", "NA18501", "NA18502", \
 "NA18504", "NA18505", "NA18507", "NA18508", "NA18510", "NA18511", "NA18516", "NA18517", "NA18519", "NA18520", "NA18522", "NA18523", \
 "NA18853", "NA18856", "NA18858", "NA18861", "NA18864", "NA18865", "NA18867", "NA18868", "NA18870", "NA18871", "NA18873", "NA18874", \
 "NA18876", "NA18877", "NA18878", "NA18879", "NA18881", "NA18907", "NA18908", "NA18909", "NA18910", "NA18912", "NA18915", "NA18916", \
 "NA18917", "NA18923", "NA18924", "NA18933", "NA18934", "NA19092", "NA19093", "NA19095", "NA19096", "NA19098", "NA19099", "NA19102", \
 "NA19107", "NA19108", "NA19113", "NA19114", "NA19116", "NA19117", "NA19118", "NA19119", "NA19121", "NA19129", "NA19130", "NA19131", \
 "NA19137", "NA19138", "NA19141", "NA19143", "NA19144", "NA19146", "NA19147", "NA19149", "NA19152", "NA19153", "NA19159", "NA19160", \
 "NA19171", "NA19172", "NA19175", "NA19184", "NA19185", "NA19189", "NA19190", "NA19197", "NA19198", "NA19200", "NA19201", "NA19204", \
 "NA19206", "NA19207", "NA19209", "NA19210", "NA19213", "NA19214", "NA19222", "NA19223", "NA19225", "NA19235", "NA19236", "NA19238", \
 "NA19239", "NA19247", "NA19248", "NA19256", "NA19257"]

def read_vcf(nameFile,chrom):
    f = gzip.open(nameFile, 'rt')
    res = []
    list_ind = []

    for line in f:       
        sline = line.split()
        if line.startswith('##'):
            continue
        elif line.startswith('#C'):    
            for i,elt in enumerate(sline):
                if i>8:
                    list_ind.append(elt)
                    var = []
                    var.append([])
                    var.append([])
                    var.append([])
                    res.append(var)
        elif read_chrom(sline[0]) == chrom:
            sline = line.split()
            for i,elt in enumerate(sline):
                if i>8:
                    if elt[0] != "." and elt[2] != ".":
                        res[i-9][0].append(int(sline[1]))
                        res[i-9][1].append(int(elt[0]))
                        res[i-9][2].append(int(elt[2]))
        elif read_chrom(sline[0]) > chrom:
            break

    return [res,list_ind]

def read_vsf_sum(nameFile,chr):
    f = gzip.open(nameFile, 'rt')
    var = []
    var.append([])
    var.append([])
    pos=[]
    if chr>9:
        c0=str(chr//10)
        c1=str(chr%10)
    else:
        c0=str(chr)
        c1="	"
        
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#C'):
            sline = line.split()
            for idx,elt in enumerate(sline):               
                if(elt in outgroup):
                    pos.append(idx)
        elif (line[0] == c0 and line[1] == c1):
            sline = line.split()
            var0=var1=False
            var[0].append(int(sline[1]))
            for p in pos:
                if(sline[p][0]=='0' or sline[p][2]=='0'):
                    var0=True
                if(sline[p][0]=='1' or sline[p][2]=='1'):   
                    var1=True
                if var0 and var1:
                    break

            if(var0 and var1):
                var[1].append([0,1])
            elif(var0):
                var[1].append([0])
            elif(var1):
                var[1].append([1])
            else:
                var[0].pop()
    f.close()    
    return var

#Same as read VCF, but also saves the MAF into the 4th position and the GP into the 5th position
def read_vcf_raf_gp(nameFile,chrom):
    
    f = open(nameFile, "r")

    if nameFile.endswith('.gz'):
        f = gzip.open(nameFile, 'rt')
    else:
        f = open(nameFile, 'rt')
   
    for line in f:
        sline=line.split()
        if line.startswith('##'):
            continue
        elif line.startswith('#C'):
            res = []
            lInd = []
            for elt in sline[9:]:
                print(elt,end=" ")
                var = []
                var.append([])
                var.append([])
                var.append([])
                var.append([])
                var.append([])
                res.append(var)
                name=elt
                if name[-1]=='\n':
                    name=name[:-1]
                lInd.append(name)
            if lInd==[]:
                return []

        elif (sline[0]==chrom):
            ssline = sline[7].split(';')
            RAF = -1
            for elt in ssline:
                if elt[:3]=='RAF':
                    RAF=float(elt[4:])
            for i,elt in enumerate(sline[9:]):
                if elt[0] != ".":
                    res[i][0].append(int(sline[1]))
                    res[i][1].append(int(elt[0]))
                    res[i][2].append(int(elt[2]))
                    res[i][3].append(RAF)
                    slineGP=(elt.split(':')[2]).split(',')
                    res[i][4].append([float(slineGP[0]),float(slineGP[1]),float(slineGP[2])])

    return [res,lInd]

def dicho_search(value, list):
    low = 0
    high = len(list) - 1
    while low <= high:
        mid = (low + high) // 2
        if list[mid] < value:
            low = mid + 1
        elif list[mid] > value:
            high = mid - 1
        else:
            return mid
    return low

def get_similarity(tract, sample, neanderthal):
    i = dicho_search(tract[0]*1000, sample[0])
    j = dicho_search(tract[0]*1000, neanderthal[0])
    same = tot = 0

    while i <len(sample[0]) and sample[0][i] < tract[1]*1000 and j < len(neanderthal[0]) and neanderthal[0][j] <tract[1]*1000:
        while i < len(sample[0]) and sample[0][i] < neanderthal[0][j]:
            i+=1

        if sample[0][i] == neanderthal[0][j]:
            genotype = [neanderthal[1][j], neanderthal[2][j]]
            if sample[1][i] in genotype or sample[2][i] in genotype:
                same += 1
            tot += 1
        j+=1
    return same / tot if tot > 0 else 'NaN'


def get_archaic_snps(tract, sample, neanderthal, outgroup):
    i = dicho_search(tract[0]*1000, sample[0])
    j = dicho_search(tract[0]*1000, neanderthal[0])
    k = dicho_search(tract[0]*1000, outgroup[0])
    same = tot = 0

    while j < len(neanderthal[0]) and neanderthal[0][j] <tract[1]*1000:
        while i < len(sample[0]) and sample[0][i] < neanderthal[0][j]:
            i+=1
        while k < len(outgroup[0]) and outgroup[0][k] < neanderthal[0][j]:
            k+=1

        if sample[0][i] == outgroup[0][k] and sample[0][i] == neanderthal[0][j]:
            archaic_allele = -1
            if sample[1][i] not in outgroup[1][k]:
                archaic_allele = sample[1][i]
            if sample[2][i] not in outgroup[1][k]:
                archaic_allele = sample[2][i]
            if archaic_allele != -1 :
                tot += 1
                if archaic_allele == neanderthal[1][j] or archaic_allele == neanderthal[2][j]:
                    same += 1
            
        j+=1
    return (same,tot)


def read_chrom(s):
    return int(s[3:]) if s.startswith('chr') else int(s)
