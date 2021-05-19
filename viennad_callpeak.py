## This script uses the Viennad file for call peak analysis.
## Hyb analyse type=mim pref=mim only output miRNA-mRNA chimera, so this script only for for mRNA as a target.
import os, re, pandas as pd, numpy as np, glob, sys, getopt, random
from datetime import datetime
random.seed('lu')
Current_time = str(datetime.now()).split()[0]

####################################################################### how to use this script 
def main(format1): ### Input format, the script needs to input Viennad file , replicates
	try:
		opts, args = getopt.getopt(format1, 'h:i:g:r:m:o:')
		for i in opts:
			if 'h' in i[0]:
				print 'python2 test.py -i <input_Viennad>  -r <replicates> -m <miRNA_name> -o <output_name>'
				sys.exit()
			elif i[0] == '-i':
				Viennad = i[1]
			elif i[0] == '-r':
				replicates1 = i[1]
			elif i[0] == '-m':
				miR_name1 = i[1]
			elif i[0] == '-o':
				output1_name = i[1]
                
	except getopt.GetoptError as err:
		print 'python2 test.py -i <input_Viennad>  -r <replicates> -m <miRNA_name> -o <output_name>'
		sys.exit()

	return Viennad, replicates1, miR_name1, output1_name

####################################################################### integrate viennad files
def genomecov_eachgene():
	global list1, list_basepattern, num_bm1
	t_line2= list1[3].split('\t')
	target_name = '_'.join(map(lambda x: list1[0].split('_')[x], range(-6,-2)))
	dG = list1[4].split('\t')[-1]
	mirna_end1 = int(list1[2].split('\t')[-1])
	mirna_start1 =  int(list1[2].split('\t')[-2])
	motif1 = list1[4].split('\t')[0][mirna_end1:].strip('.')
	motif2 = ''.join(map(lambda x: '\\'+x,motif1))
	rex1 = re.compile(motif2)
	modol_seq = re.search(rex1,list1[4].split('\t')[0])
	motif_seq = list1[1].split('\t')[0][modol_seq.span()[0]:modol_seq.span()[1]] ## extract motif sequence
	list_basepattern.append('>'+ str(num_bm1)  +'@' + list1[0]+'_dG:'+str(dG[1:-1])+'\t'+list1[4].split('\t')[0][mirna_start1-1:mirna_end1]+'\t'+list1[4].split('\t')[0][mirna_end1:].strip('.')+'\t'+motif_seq)
	num_bm1 += 1
	dis1 = (int(t_line2[3])+1)-int(t_line2[2])
	dfm = pd.DataFrame({'gene': [target_name]*dis1,'position': range(int(t_line2[2]),int(t_line2[3])+1),'Number': [1]*dis1})
	dfm = dfm.set_index(['gene', 'position',])
	return dfm

def integrate_genomecov(files):
	global list1, list_basepattern, num_bm1
	num_bm1=0  ## calculate total number of reads in peak position
	dG,F = 0,0 ## Calculate the minimum free energy
	list1, list_basepattern = [],[] ##Temporarily get information about each target
	dis1=0 ##The distance between the beginning and the end of the reads
	f=pd.DataFrame() ## Create empty DataFrame 
	dict1={} ##Integrate genomecov_peak files into the dictionary
	if len((glob.glob(files)))<= 1:
		print 'the total Viennad file is ' +str(len(glob.glob(files)))
	else:
		print 'the total Viennad file are ' +str(len(glob.glob(files)))
	for i in glob.glob(files):
		print i
        	with open(i,'r+') as f1:
			for line1 in f1:
				if ('microRNA' in line1) and ('mRNA' in line1) and (miR_name1.upper().replace('ALL','RNA') in line1.upper()):
					if (list1 != []) and ('microRNA' in list1[0]) and ('mRNA' in list1[0]) and  (miR_name1.upper().replace('ALL','RNA')  in list1[0].upper()):
						dfm = genomecov_eachgene()
						for all_re in range(int(list1[0].split('_')[1])): ## expand the reads
							f= pd.concat([f,dfm],axis=1).sum(axis=1)  ## f indicate the total replicates
					list1=[]
				list1.append(line1.strip())
			if (list1 != []) and ('microRNA' in list1[0]) and ('mRNA' in list1[0]) and (miR_name1.upper().replace('ALL','RNA') in list1[0].upper()):
				dfm = genomecov_eachgene()
				for all_re in range(int(list1[0].split('_')[1])):## expand the reads
					f = pd.concat([f,dfm],axis=1).sum(axis=1).rename(i)
				dict1[i]=f
				list1 = [] ## Empty the list after analyzing the contents of a Viennad file
				f=pd.DataFrame() ## Create empty DataFrame
	f = pd.concat(map(lambda x: dict1[x],dict1.keys()),axis=1)
	Rep1 = f.count(axis=1).rename('Replicate')
	Peak1 = f.sum(axis=1).rename('Peak')
	f = pd.concat([f,Rep1,Peak1],axis=1)
	f = f.sort_index()
	f = f.dropna(axis=0,how='all')
	f = f[f['Replicate']>= int(replicates1)]  ## Enter many replicates you need?
	return f, list_basepattern

#######################################################################  human hOH7.fasta file
def transcript_database():  
	with open ('/blue/mingyi.xie/luli1/hybdb/db/hOH7_20210422.fasta') as f1: 
##	with open ('hOH7.fasta') as f1:
		dict1 = {}
		name1,seq1='',''
		for line1 in f1:
			if line1[0] == '>':
				if name1 != '':
					dict1[name1] = seq1
				name1,seq1='',''
				name1= line1.strip()[1:]
			else:
				seq1 += line1.strip().upper()
		dict1[name1] = seq1
	return dict1

#######################################################################  Call_Peak
def call_peak():
	f_out = open(str(Current_time)+'_'+str(output1_name)+'/Call_peak.txt','w+')
	f_out.write('gene name'+'\t'+'start'+'\t'+'end'+'\t'+'Replicate'+'\t'+'Peak'+'\t'+ 'sequence'+'\n')
	gene = '' ## Temporary gene when analyzing replicate
	list_Replicate, list_Peak = [], []
	start1,end1 = 0,0 ## The beginning and end of Replicate
	seq1='' ## Extract a partial sequence of a specific mRNA
	for x,y in f1[0].index:
		if (end1+1) == int(y) and gene == x:
			end1 = int(y)
		else:
			if start1 != 0:
				f_out.write(gene+'\t'+str(start1)+'\t'+str(end1)+'\t'+str(int(max(list_Replicate)))+'\t'+str(int(max(list_Peak)))+'\t'+ str(dict_trans[gene][start1-1:end1])+'\n')
			list_Replicate, list_Peak = [], []
			gene = x
			start1, end1 = int(y), int(y)
		list_Replicate.append(f1[0].loc[x,y]['Replicate'])
		list_Peak.append(f1[0].loc[x,y]['Peak'])
	f_out.write(gene+'\t'+str(start1)+'\t'+str(end1)+'\t'+str(int(max(list_Replicate)))+'\t'+str(int(max(list_Peak)))+'\t'+ str(dict_trans[gene][start1-1:end1])+'\n')
	f_out.close()

####################################################################### basepattern and motif result
def output_symbol_motif():  ## output basepattern and motif result
	list1_p1, list1_m1 = [],[] ## create base pattern and motif list 
	with open(str(Current_time)+'_'+str(output1_name)+'/Call_peak.txt','r+') as f_peak:
		for line1 in f_peak:
			line2 = line1.strip().split('\t')
			for i in f1[1]:  ## loop of base pattern and motif information
				if line2[0] in i:  ## Determine if the gene name is in the base pattern 
					beginmotif1,endmotif1 = i.split('\t')[0].split('_')[-3],i.split('\t')[0].split('_')[-2] ## The number of position in listbasepattern
					if set(range(int(beginmotif1),int(endmotif1)+1)) & set(range(int(line2[1]),int(line2[2])+1)) != set(): ## Determine if there are intersect between Peak and base pattern
						list1_p1.append(i.split('\t')[0] + '\t' + i.split('\t')[1])
						list1_m1.append(i.split('\t')[0] + '\n' + i.split('\t')[3])
	f1_symbol = open(str(Current_time)+'_'+str(output1_name)+'/Symbol.txt','w+')
	list_symbol = []
	set1_p1 = set(list1_p1) ## remove redundant base pattern
	set1_m1 = set(list1_m1) ## remove redundant motif
	for p in random.sample(set1_p1,len(set1_p1)):  
		dG1 =  float(p.split('\t')[0].split(':')[-1])
		motif1 = p.split('\t')[1]
		for x in motif1.replace('.','0').replace('(','1'):
			if x == '0':
				list_symbol.append(x)
			else:
				F = -dG1-11
				if F>=5:
					list_symbol.append('5')
				elif F <= 0:
					list_symbol.append('0')
				else:
					list_symbol.append(str(F))
		for all_p_re in range(int(p.split('_')[1])): ## expand base pattern
			f1_symbol.write(p.split('\t')[0]+'\t'+'\t'.join(list_symbol)+'\n')
		list_symbol = []

	f1_motif = open(str(Current_time)+'_'+str(output1_name)+'/motif.txt','w+')
	for m in set1_m1:
		for all_m_re in range(int(m.split('_')[1])): ##expand motif
			f1_motif.write(m+'\n')

	f1_motif.close()
	f1_symbol.close()
    
#######################################################################
if __name__ == '__main__':
	Viennad, replicates1, miR_name1, output1_name = main(sys.argv[1:])   ### remember write the number of replicate
	os.mkdir(str(Current_time)+'_'+str(output1_name))
	dict_trans = transcript_database()  ## good
	f1 = integrate_genomecov(Viennad)
	call_peak()
	output_symbol_motif()
