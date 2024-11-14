import numpy as np
import sys

##open files
zea = sys.argv[1]
vcf = sys.argv[2]
chrom = sys.argv[3]
pop = sys.argv[4]

###dictionary with the positions of gerp estimated
gpns = open(zea, 'r')
ancestrals = {}
for line in gpns:
	line2 = line.split('\t')
	#print(line2)
	key = line2[0].split(' ')[1]
	#print(key)
	snp = line2[4][0:-1] #use the sorghum to polarize
	ancestrals[key] = [snp]

##read the fixed snps list
fix = open('./Data/{}.fixed.frq'.format(pop))
fixed = {}
for line in fix:
	line2 = line.split('\t')
	if line2[0] == chrom: #keep only the current chromosome
		key = line2[1]
		if line2[4].split(":")[1][0] == '1':
			#print(line2[4].split(":")[1][0])
			fixed[key] = line2[4][0]
		else:
			try:
				fixed[key] = line2[5][0]
			except IndexError:
				fixed[key] = "NA"	
	else:
		pass

file3 = open('fixed_load_{}_{}'.format(chrom, sys.argv[4]), 'w')
file4 = open('segregating_load_{}_{}'.format(chrom, sys.argv[4]), 'w')

#open the vcf file
vcfIO = open(vcf, 'r')
for line in vcfIO:
	if line[1] == '#':
		pass
	elif line[1] == 'C':
		line2 = line.split('\t')
		gp = line2[9::]
		file3.write('POS\n')
		file4.write('POS\n')
	else:
		line2 = line.split('\t')
		key = str(line2[1][0::])
		#print(line2[1], pos)
		ref = line2[3]
		alt = line2[4]
		try:
			if ref == ancestrals[key][0]: #ref the ancestral and alt the derived
				if alt == fixed[key]:
					file3.write('{}\n'.format(key))
				elif ref == fixed[key]:
					pass		
			elif alt == ancestrals[key][0]:
				if ref == fixed[key]:
					file3.write('{}\n'.format(key))
				elif alt == fixed[key]:
					pass
			else:
				pass
		except KeyError:
			try:
				tmp = ancestrals[key][0]
				file4.write('{}\n'.format(key))
			except KeyError:
				pass

file3.close()
file4.close()
sys.exit("Live long and prosper!")