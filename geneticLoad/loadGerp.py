import math
import numpy as np
import math
import sys

##open files
zea = sys.argv[1]
vcf = sys.argv[2]
chrom = sys.argv[3]

###dictionary with the positions of gpn estimated
gpns = open(zea, 'r')
values = {}
for line in gpns:
	line2 = line.split('\t')
	#print(line2)
	key = line2[0].split(' ')[1]
	#print(key)
	snp = line2[4][0:-1] #use the sorghum to polarize
	if snp == '-':
		pass
	else:
		values[key] = [snp, line2[2]]

#open and filter the gerp files
vcfIO = open(vcf, 'r')
file3 = open("./Scores/gerp_scores_{}".format(chrom), 'w')
#open the vcf file
for line in vcfIO:
	if line[1] == '#':
		pass
	elif line[1] == 'C':
		line2 = line.split('\t')
		gp = line2[9::]
		#print(len(gp))
		file3.write('POS\t{}'.format('\t'.join(gp)))
	else:
		line2 = line.split('\t')
		key = str(line2[1][0::])
		#print(line2[1], pos)
		ref = line2[3]
		alt = line2[4]
		gp = line2[9::]
		gpns2 = []
		gerp2 = 0
		try:
			if ref == values[key][0]:
				gerp2 = values[key][1]
				#print(gp)
				for u in gp:
					i = u.split(":")
					#print(i, i[0])
					if i[0] in ['0/0', '0/0\n']:
						gpns2.append(str(0))
					elif i[0] in ['0/1','1/0', '0/1\n','1/0\n']:
						gpns2.append(str(gerp2))
					elif i[0] in ['1/1', '1/1\n']:
						gpns2.append(str(float(gerp2)*2))
					else:
						gpns2.append('NA')
						#break
				#print(ref, key, values[key], gp, gpns2)
				file3.write('{}\t{}\n'.format(key,'\t'.join(gpns2)))
			elif alt == values[key][0]:
				gerp2 = values[key][1]
				for u in gp:
					i = u.split(":")
					if i[0] in ['1/1', '1/1\n']:
						gpns2.append(str(0))
						#print(gpns)
					elif i[0] in ['0/1','1/0', '0/1\n','1/0\n']:
						gpns2.append(str(gerp2))
						#print(values[pos])
					elif i[0] in ['0/0', '0/0\n']:
						gpns2.append(str(2*float(gerp2)))
					else:
						gpns2.append('NA')
				#print(len(gpns2))
				#print(alt, key, values[key], gp, gpns2)
				file3.write('{}\t{}\n'.format(key,'\t'.join(gpns2)))
			else:
				gerp2 = 0
				pass
		except KeyError:
			pass

file3.close()
sys.exit("Live long and prosper!")