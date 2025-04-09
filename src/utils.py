import os
import gzip

def createMap(decodePath, outPath):
    fou = open(outPath,'w')
    for filename in os.listdir(decodePath):
        file_path = os.path.join(decodePath, filename)
        ind = (filename.split('.')[-1]).split('-')[0]
        f = open(file_path,'r')
        firstLine = next(f)
        if os.stat(outPath).st_size == 0:
            fou.write(firstLine[:-1])
            fou.write('\tIndividual\n')
        for line in f:
            sline = line.split()
            if sline[4] == 'Neanderthal':
                fou.write(line[:-1])
                fou.write(f'\t{ind}\n')
    fou.close()
    f.close()

def computeSimilarities(decodePath,mapped_individual,reference_individual,name_reference):
    for filename in os.listdir(decodePath):
        print(filename)
        file_path = os.path.join(decodePath, filename)
        name_mapped= filename.split('.')[0]
        for vcf_reference_individual in reference_individual:
            print(vcf_reference_individual)
            if mapped_individual.endswith('.gz'):    
                f1 = gzip.open(mapped_individual, 'rt')
            else:
                f1 = open(mapped_individual,'r')


            if vcf_reference_individual.endswith('.gz'):    
                f2 = gzip.open(vcf_reference_individual, 'rt')
            else:
                f2 = open(vcf_reference_individual,'r')


            map =  open(file_path, 'r') 
            processed_lines = []
            lineMap = next(map)
            slineMap = lineMap.split()

            if slineMap[-1] != f'similarity{name_reference}\n':
                processed_lines.append(f'{lineMap[:-1]}\tsimilarity{name_reference}\n')
            else:
                processed_lines.append(lineMap)

            line1=next(f1)
            while line1[1] == '#':
                line1= next(f1)
            line2=next(f2)
            while line2[1] == '#':
                line2= next(f2)

            sline1 = line1.split()
            try:
                index1 = sline1.index(name_mapped)
            except ValueError:
                print(f'Name of the mapped individual ({name_mapped}) is not in the vcf file.')
            sline2 = line2.split()
            try:
                index2 = sline2.index(name_reference)
            except ValueError:
                print(f'Name of the reference individual ({name_reference})is not in the vcf file.')

            match = 0
            total = 0
            lineMap = next(map)
            slineMap = lineMap.split()
            sline2 = next(f2).split()
            for line1 in f1:
                sline1 = line1.split()
                chr = int(sline1[0])
                try:
                    if int(sline1[1]) > int(slineMap[2]) and int(slineMap[0]) == chr and int(sline2[0]) ==  chr:
                        currLine = lineMap[:-1]
                        if total>0:
                            processed_lines.append(f'{currLine}\t{match/total:.4f}\n') 
                        else:
                            processed_lines.append(f'{currLine}\n') 
                        lineMap = next(map)
                        slineMap = lineMap.split()    
                        total=0
                        match=0
                    elif int(slineMap[0])<chr:
                        while int(slineMap[0])<chr:
                            lineMap = next(map)
                            slineMap = lineMap.split() 

                except StopIteration:
                    # End of file reached
                    break
                if chr >= int(sline2[0]):
                    try:
                        while (int(sline2[1]) < int(sline1[1]) and int(sline2[0]) ==  chr) or int(sline2[0]) <  chr:
                                sline2 = next(f2).split() 
                        if  int(sline1[1])==int(sline2[1]) and int(sline2[0]) ==  chr:
                            if sline1[index1][0] != '.' and sline2[index2][0] != '.': 
                                total += 1
                                gen=[]
                                if sline1[index1][0] == '0' or sline1[index1][2] =='0':
                                    gen.append('0')
                                if sline1[index1][0] == '1' or sline1[index1][2] =='1':
                                    gen.append('1')
                                if sline2[index2][0] in gen or sline2[index2][2] in gen:
                                    match += 1
                    except StopIteration:
                        # End of file reached
                        break
            processed_lines.append(lineMap)
            for lineMap in map:
                processed_lines.append(lineMap)
            with open(file_path, 'w') as outfile:
                outfile.writelines(processed_lines)
            f1.close()
            f2.close()