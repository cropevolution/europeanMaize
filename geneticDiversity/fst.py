#modified from Kerstins file

import numpy as np
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas as pd
import matplotlib.pyplot as plt
import allel
import random
import h5py
import bcolz
import csv
from itertools import chain
import itertools
import os
from collections import Counter
from itertools import combinations
import seaborn as sns; sns.set_theme(color_codes=True)
from matplotlib import pyplot


##load the data
callset = allel.read_vcf('./vcf/maizeEU.noDH.filteredQual.recode.vcf')

#get the samples in the .vcf
samples= callset["samples"]
len(samples)

#put all the samples in a Dataframe
sample_subset = pd.DataFrame(data=samples, columns = ['sample_names'])

#get the name of the landraces from the vcf file
mayer_2017_landraces_samples=sample_subset[sample_subset.sample_names.str.contains("_2017_landraces")]
print(len(mayer_2017_landraces_samples))
mayer_2017_35_landraces=mayer_2017_landraces_samples[mayer_2017_landraces_samples.sample_names.str.contains("\\.01")]
print(len(mayer_2017_35_landraces))

#add the populations to the individuals
for index, row in mayer_2017_landraces_samples.iterrows():
    #print(row['sample_names'])
    mayer_string=row['sample_names']#take the sample names row
    #print(mayer_string) 
    mayer_substring=mayer_string[0:2] #take the first two letters of the sample names
    #print(mayer_substring)
    mayer_2017_landraces_samples.at[index,'Population'] = mayer_substring + "_mayer_2017" #add the dataset info and make the Population column

#get all the unique population that you have assigned and check the number 
populations_mayer = mayer_2017_landraces_samples.Population.unique()
print(len(populations_mayer))
print(populations_mayer)

#Get all the combinations of the different populations for the pairwise-Fst
two_row_dataframe = {}
splits = list(mayer_2017_landraces_samples.groupby("Population"))
for i in range(35):
    two_row_dataframe[i] = [] 
    # view splitted dataframe
    #print(splits[i][1])
    two_row_dataframe[i]=splits[i][1]

##
first_row_dataframe={} 
second_row_dataframe={}
for i in range(35):
    first_row_dataframe[i]=[]
    second_row_dataframe[i]=[]
    # view splitted dataframe
    #print(two_row_dataframe[i].sample_names)
    #print(two_row_dataframe[i].Population)
    first_row_dataframe[i]=two_row_dataframe[i].sample_names
    second_row_dataframe[i]=two_row_dataframe[i].Population

###make the combinations by numbers
#make a list of the length 35 from 0 to 34
L=list(range(35))
print(len(L))
print(L)

#make a list of all combinations - **including the reciprocal (O,1) and (1,0)** because they are needed for the heatmap
#combo=list(itertools.permutations(L,2))
##Change tp have only once the comparison
combo=list(itertools.combinations(L,2))
print(len(combo))
#print(combo)

####append the self ones to the combo list 
combo.append((0,0))
combo.append((1,1))
combo.append((2,2))
combo.append((3,3))
combo.append((4,4))
combo.append((5,5))
combo.append((6,6))
combo.append((7,7))
combo.append((8,8))
combo.append((9,9))
combo.append((10,10))
combo.append((11,11))
combo.append((12,12))
combo.append((13,13))
combo.append((14,14))
combo.append((15,15))
combo.append((16,16))
combo.append((17,17))
combo.append((18,18))
combo.append((19,19))
combo.append((20,20))
combo.append((21,21))
combo.append((22,22))
combo.append((23,23))
combo.append((24,24))
combo.append((25,25))
combo.append((26,26))
combo.append((27,27))
combo.append((28,28))
combo.append((29,29))
combo.append((30,30))
combo.append((31,31))
combo.append((32,32))
combo.append((33,33))
combo.append((34,34))

#sort the combo list
combo.sort()

#calculate Fst
combis=combo
combis.sort()

#gets for each comparison which individuals are in that list
combi_list= {}
for i in range(630):
    combis_element=combis[i]
    combis_fist=combis_element[0]
    combis_second=combis_element[1]
    #print(combis_element)
    #print(combis_fist)
    #print(combis_second)
    list_first=(list(first_row_dataframe[combis_fist]))
    list_second=(list(first_row_dataframe[combis_second]))
    #print(list_first)
    #print(list_second)
    combi_list[i]=[]
    combi_list[i]=list_first+list_second
    #print(combi_list)
    #print(len(combi_list))

#set up the lists for the Fst values and the corresponding Fst standard deviations
fst_list=[]
fst_se_list=[]

#set up a list of all Mayer_2017 landraces names for the print statement of the calculations
unique_populations_mayer=pd.unique(populations_mayer)
print(len(unique_populations_mayer))
print(unique_populations_mayer)

#calculate the pairwise Fst values and corresponding deviations 
file = open('fsts_maizeEU.txt', 'w')
for k in range(630):
    combis_element=combis[k]
    combis_fist=combis_element[0]
    first_landrace_name=unique_populations_mayer[combis_fist]
    combis_second=combis_element[1]
    second_landrace_name=unique_populations_mayer[combis_second]
    combi_list_element=combi_list[k]
    #print(combis_element)
    #print(combis_fist)
    #print(combis_second)
    print(first_landrace_name)
    print(second_landrace_name)
    #print(combi_list_element)
    
    #load the data for the two landraces seperatly 
    loc_samples_pop1 = allel.read_vcf('./vcf/maizeEU.noDH.filteredQual.recode.vcf', samples= first_row_dataframe[combis_fist])
    loc_samples_pop2 = allel.read_vcf('./vcf/maizeEU.noDH.filteredQual.recode.vcf', samples= first_row_dataframe[combis_second])
    
    #load the data for the combination of the two samples 
    loc_samples_combi_pops = allel.read_vcf('./vcf/maizeEU.noDH.filteredQual.recode.vcf', samples= combi_list_element)  
    
    #get the genotype data for the two landraces seperatly
    geno_samples_pop1=allel.GenotypeArray(loc_samples_pop1['calldata/GT'])
    geno_samples_pop2=allel.GenotypeArray(loc_samples_pop2['calldata/GT'])
    
    #get the genotype data for the combination
    geno_samples_combi_pop=allel.allel.GenotypeArray(loc_samples_combi_pops['calldata/GT'])
    
    #get the allel count array for the two landraces seperatly
    geno_samples_pop1_allel_count=geno_samples_pop1.count_alleles()
    geno_samples_pop2_allel_count=geno_samples_pop2.count_alleles()
    
    #get the allel count array for the combination of the two landraces
    geno_samples_combi_pop_allel_count=geno_samples_combi_pop.count_alleles()
    
    #boolean to check if the selected individuals segregeate at the SNPs
    is_seg = geno_samples_combi_pop_allel_count.is_segregating()[:]
    
    #apply the boolean to the allel count arrays 
    geno_samples_pop1_allel_count_seg=geno_samples_pop1_allel_count.compress(is_seg, axis=0)
    geno_samples_pop2_allel_count_seg=geno_samples_pop2_allel_count.compress(is_seg, axis=0)
    
    
    #calculate the Fst
    fst, fst_se, _, _ = allel.blockwise_hudson_fst(geno_samples_pop1_allel_count_seg,geno_samples_pop2_allel_count_seg, blen=100000)
    
    print(fst)
    print(fst_se)
    #append fst and fst_se to list
    fst_list.append(fst)
    fst_se_list.append(fst_se)
    
    #print("Fst for %s combination - landraces %s/%s was calculated" % (str(combis_element)))
    print("Fst for %(a)s combination - landraces %(b)s/%(c)s was calculated" % {'a': str(combis_element), 'b': first_landrace_name, 'c':second_landrace_name})
    file.write('{}\t{}\{}\t{}\n'.format(first_landrace_name,second_landrace_name, fst, fst_se))
file.close()


#