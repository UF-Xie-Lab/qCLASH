infile1 = open('.txt', 'r')

#takes an .viennad input file where all six lines of the viennad had been joined to a single line

newlist = []
a = 0 #8mer
b = 0 #7mer-A1
c = 0 #7mer-m8
d = 0 #6mer
e = 0 #6mer-off
f = 0 #mismatch
g = 0 #other

for line in infile1:
    line = line.strip('\n').split(' ')
    newlist.append(line)

for line in newlist:
    try:
        element = line[4].split('\t')[0].rfind(')')
        mRNA_position = int(element) + 1
        if line[4][1:8] == '(((((((' and line[1][mRNA_position]== 'A':
            a += 1
        elif line[4][1:7] == '((((((' and line[1][mRNA_position]== 'A':
            b += 1
        elif line[4][1:8] == '(((((((':
            c += 1
        elif line[4][1:7] == '((((((':
            d += 1
        elif line[4][2:8] == '((((((':
            e += 1
        elif line[4][1:8] == '.((((((' or line[4][1:8] == '(.(((((' or line[4][1:8] == '((.((((' or line[4][1:8] == '(((.(((' or line[4][1:8] == '((((.((' or line[4][1:8] == '(((((.(' or line[4][1:8] == '((((((.':
            f += 1
        else:
            g +=1
    except:
		None
        
print(infile1.name)
print(f'8mer = {a}\n7mer-A1 = {b}\n7mer-m8 = {c}\n6mer = {d}\n6mer-off = {e}\nmismatch = {f}\nother={g}')
print(a+b+c+d+e+f+g)
infile1.close()
