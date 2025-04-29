import os
import gzip
import json
import re

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
    filename = decodePath.split('/')[-1]
    name_mapped= filename.split('.')[0]
    print(filename)
    for vcf_reference_individual in reference_individual:
        print(vcf_reference_individual)
        if mapped_individual.endswith('.gz'):    
            f1 = gzip.open(mapped_individual, 'rt')
        else:
            f1 = open(mapped_individual,'r')


        if vcf_reference_individual.endswith('.gz'):    
            file_reference = gzip.open(vcf_reference_individual, 'rt')
        else:
            file_reference = open(vcf_reference_individual,'r')


        map =  open(decodePath, 'r') 
        processed_lines = []
        lineMap = next(map)
        slineMap = lineMap.split()

        if slineMap[-1] != f'similarity{name_reference}':
            processed_lines.append(f'{lineMap[:-1]}\tsimilarity{name_reference}\n')
        else:
            processed_lines.append(lineMap)

        line1=next(f1)
        while line1[1] == '#':
            line1= next(f1)
        line2=next(file_reference)
        while line2[1] == '#':
            line2= next(file_reference)

        sline1 = line1.split()
        try:
            index1 = sline1.index(name_mapped)
        except ValueError:
            print(f'Name of the mapped individual ({name_mapped}) is not in the vcf file.')
            break
        sline2 = line2.split()
        try:
            index2 = sline2.index(name_reference)
        except ValueError:
            print(f'Name of the reference individual ({name_reference})is not in the vcf file.')

        match = 0
        total = 0
        lineMap = next(map)
        slineMap = lineMap.split()
        sline2 = next(file_reference).split()
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
                        processed_lines.append(lineMap)
                        lineMap = next(map)
                        slineMap = lineMap.split() 

            except StopIteration:
                # End of file reached
                break
            if chr >= int(sline2[0]):
                try:
                    while (int(sline2[1]) < int(sline1[1]) and int(sline2[0]) ==  chr) or int(sline2[0]) <  chr:
                            sline2 = next(file_reference).split() 
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
                    processed_lines.append(lineMap)
                    for lineMap in map:
                        processed_lines.append(lineMap)
                    break

        with open(decodePath, 'w') as outfile:
            outfile.writelines(processed_lines)
        f1.close()
        file_reference.close()
        map.close()

def getChrom(s):
    return int(s[3:]) if s.startswith("chr") else int(s)


def computeDistances(decodePath,mapped_individual,reference_individual,name_reference, outgroup_individuals, names_outgroup):
    if names_outgroup.endswith('.json'):    
        with open(names_outgroup, 'r') as file:
            list_names_outgroup = json.load(file)
    else:
        print(f'File {names_outgroup} is not a json file.')
        return

    filename = decodePath.split('/')[-1]
    name_mapped= filename.split('.')[0]
    print(filename)

    reference_individual = sorted(reference_individual, key=extract_chromosome_from_name)
    outgroup_individuals = sorted(outgroup_individuals, key=extract_chromosome_from_name)

    for (vcf_reference_individual, vcf_outgroup) in zip(reference_individual, outgroup_individuals):
        print(vcf_reference_individual)
        print(vcf_outgroup)
        if mapped_individual.endswith('.gz'):    
            file_mapped = gzip.open(mapped_individual, 'rt')
        else:
            file_mapped = open(mapped_individual,'r')


        if vcf_reference_individual.endswith('.gz'):    
            file_reference = gzip.open(vcf_reference_individual, 'rt')
        else:
            file_reference = open(vcf_reference_individual,'r')

        if vcf_outgroup.endswith('.gz'):    
            file_outgroup = gzip.open(vcf_outgroup, 'rt')
        else:
            file_outgroup = open(vcf_outgroup,'r')


        map =  open(decodePath, 'r') 
        processed_lines = []
        lineMap = next(map)
        slineMap = lineMap.split()

        if slineMap[-1] != f'distance{name_reference}':
            processed_lines.append(f'{lineMap[:-1]}\tdistance{name_reference}\n')
        else:
            processed_lines.append(lineMap)

        line1=next(file_mapped)
        while line1[1] == '#':
            line1= next(file_mapped)

        line_reference=next(file_reference)
        while line_reference[1] == '#':
            line_reference= next(file_reference)

        line_outgroup=next(file_outgroup)
        while line_outgroup[1] == '#':
            line_outgroup= next(file_outgroup)


        sline1 = line1.split()
        try:
            index1 = sline1.index(name_mapped)
        except ValueError:
            print(f'Name of the mapped individual ({name_mapped}) is not in the vcf file.')
        sline_reference = line_reference.split()
        try:
            index_reference = sline_reference.index(name_reference)
        except ValueError:
            print(f'Name of the reference individual ({name_reference})is not in the vcf file.')
        sline_outgroup = line_outgroup.split()
        index_outgroup = []
        for elt in list_names_outgroup:
            try:
                index = sline_outgroup.index(elt)
                index_outgroup.append(index)
            except ValueError:
                print(f'Name of the outgroup individual ({elt}) is not in the vcf file.')

        match = 0
        total = 0
        match_outgroup = [0 for _ in range(len(index_outgroup))]
        total_outgroup = [0 for _ in range(len(index_outgroup))]
        lineMap = next(map)
        slineMap = lineMap.split()
        sline_reference = next(file_reference).split()
        sline_outgroup = next(file_outgroup).split()
        for line1 in file_mapped:
            sline1 = line1.split()
            chr = int(sline1[0])
            try:
                if int(sline1[1]) > int(slineMap[2]) and int(sline_outgroup[1]) > int(slineMap[2]) and int(slineMap[0]) == chr and \
                getChrom(sline_outgroup[0]) == chr and int(sline_reference[0]) ==  chr:
                    currLine = lineMap[:-1]
                    if total>0 and sum(elt > 0 for elt in total_outgroup) and sum(match_outgroup)>0:
                        average_distance_outgroup = sum([1-a/b for a,b in zip(match_outgroup, total_outgroup)])/len(match_outgroup)
                        
                        if average_distance_outgroup == 0:
                            processed_lines.append(f'{currLine}\tNaN\n') 
                        else:
                            processed_lines.append(f'{currLine}\t{((1-match/total)/average_distance_outgroup):.4f}\n') 
                    else:
                        processed_lines.append(f'{currLine}\n') 
                    lineMap = next(map)
                    slineMap = lineMap.split()    
                    total=0
                    match=0
                    match_outgroup = [0 for _ in range(len(index_outgroup))]
                    total_outgroup = [0 for _ in range(len(index_outgroup))]
                elif int(slineMap[0])<chr:
                    while int(slineMap[0])<chr:
                        processed_lines.append(lineMap)
                        lineMap = next(map)
                        slineMap = lineMap.split() 

            except StopIteration:
                # End of file reached
                break
            if chr >= int(sline_reference[0]) and chr >= getChrom(sline_outgroup[0]):
                try:
                    while (int(sline_reference[1]) < int(sline1[1]) and int(sline_reference[0]) ==  chr) or int(sline_reference[0]) <  chr:
                            sline_reference = next(file_reference).split() 
                    while (int(sline_outgroup[1]) < int(sline1[1]) and getChrom(sline_outgroup[0]) ==  chr) or getChrom(sline_outgroup[0]) <  chr:
                            sline_outgroup = next(file_outgroup).split()
                    
                    if  int(sline1[1])==int(sline_reference[1]) and int(sline_reference[0]) ==  chr: 
                        if sline1[index1][2] != '.' and sline_reference[index_reference][2] != '.': 
                            total += 1
                            gen=[]
                            if sline1[index1][0] == '0' or sline1[index1][2] =='0':
                                gen.append('0')
                            if sline1[index1][0] == '1' or sline1[index1][2] =='1':
                                gen.append('1')
                            if sline_reference[index_reference][0] in gen or sline_reference[index_reference][2] in gen:
                                match += 1
                    if int(sline_outgroup[1])==int(sline_reference[1]) and getChrom(sline_outgroup[0]) ==  int(sline_reference[0]):
                        for i,idx in enumerate(index_outgroup):
                            if sline_outgroup[idx][2] != '.' and sline_reference[index_reference][2] != '.':
                                total_outgroup[i] += 1
                                gen=[]
                                if sline_outgroup[idx][0] == '0' or sline_outgroup[idx][2] =='0':
                                    gen.append('0')
                                if sline_outgroup[idx][0] == '1' or sline_outgroup[idx][2] =='1':
                                    gen.append('1')
                                if sline_reference[index_reference][0] in gen or sline_reference[index_reference][2] in gen:
                                    match_outgroup[i] += 1
                except StopIteration:
                    processed_lines.append(lineMap)
                    for lineMap in map:
                        processed_lines.append(lineMap)
                    break
        with open(decodePath, 'w') as outfile:
            outfile.writelines(processed_lines)
        file_mapped.close()
        file_reference.close()
        map.close()

def extract_chromosome_from_name(s):
    match = re.search(r'chr(\d+)', s)
    if match:
        return int(match.group(1))
    return -1