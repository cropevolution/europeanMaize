#remove the list of snps I do not need
vcf = './merge_merge_all_chrom_all_edits_with_allel_column_for_masters_vcf.vcf'
snps = 'filteredSNPs'
keeps = []
snps = open(snps, 'r')
for l in snps:
	keeps.append(l[0:-1])
print("filtering the vcf file...")
vcfIO = open(vcf, 'r')
vcfNew = open('maizeEU.filt.python.vcf', 'w')
for line in vcfIO:
	if line[0] == '#':
		vcfNew.write(line)
	else:
		line2 = line.split('\t')
		#print(line2)
		#print(line)
		if line2[2] in keeps:
			vcfNew.write(line)
		else:
			pass
vcfNew.close()