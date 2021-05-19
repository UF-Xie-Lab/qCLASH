infile1 = open('file.hyb', 'r')
infile2 = open('filteredHuman_mRNAs.txt', 'r')

#input requires a hyb file along with a Ensemble Biomart assembled file for the start and stop locations of each gene id

mergedlist = []
database = {}
final = []
newlist3 = {}

for line in infile1.readlines():
    line = line.strip('\n').split('\t')
    mergedlist.append(line)

for line in infile2.readlines():
    line = line.strip('\n').split('\t')
    if line[0] not in database:
        database[line[0]] = line[2:]
    else:
        database[line[0]] = database[line[0]] + (line[2:])

for x,y in database.items():
    newlist = sorted(list(map(lambda x: int(x), y)))
    max1_list = max(newlist)
    min1_list = min(newlist)

    newlist3[x] = [min(newlist),max(newlist)]

for line in mergedlist:
    if line[9].split('_')[1] in newlist3:
        gene_id = line[9].split('_')[1]
        final.append(line + newlist3[gene_id])

a = 0 #placeholder for 5'UTR
b = 0 #placeholder for 5' CDS
c = 0 #placeholder for CDS
d = 0 #placeholder for 3' CDS
e = 0 #placeholder for 3' UTR

for line in final:
    transcript_start = int(line[12])
    transcript_end = int(line[13])
    CDS_start = int(line[-2])
    CDS_end = int(line[-1])
    if transcript_end < CDS_start:
        a += 1
    elif transcript_start < CDS_start and CDS_end > transcript_end > CDS_start:
        b += 1
    elif transcript_start > CDS_start and CDS_end > transcript_end > CDS_start:
        c += 1
    elif CDS_end > transcript_start > CDS_start and CDS_end < transcript_end:
        d += 1
    elif transcript_start > CDS_end:
        e += 1
    else:
        None

print(infile1.name)        
print(f"5'UTR = {a}\n5'UTR CDS = {b}\nCDS = {c}\n3'UTR CDS = {d}\n3'UTR= {e}")        

infile1.close()
infile2.close()
