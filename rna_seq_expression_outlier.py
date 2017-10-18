#import argparse
import sys
import numpy as np
from scipy import stats

input_files=open(sys.argv[1],'r')

gene_list_file=open(sys.argv[2],'r')
gene_list=[]
for gene in gene_list_file:
	gene_list.append(gene.strip())
gene_list_file.close()

output_file=open(sys.argv[3],'w')
outlier_levels=int(sys.argv[4])
input_file_header=sys.argv[5]
gene_column=int(sys.argv[6])
expr_column=int(sys.argv[7])

gene_dict = {}
for gene in gene_list:
	if gene in gene_dict:
		next
	else:
		gene_dict[gene] = [[],None,None,None,None,None]

sample_list = []

for line in input_files:
	line = line.strip().split()
	sample = line[0]
	print("Sample "+sample)
  	sample_file = open(line[1], 'r')
        if input_file_header:
		sample_file.readline()
  	if sample in sample_list:
  		sys.exit("Sample name "+sample+" appears at least twice in sample file list input file.")
  	else:
  		sample_list.append(sample)
  	sample_dict = {}
	for sline in sample_file:
  		sline = sline.strip().split()
  		gene = sline[gene_column]
		expr = sline[expr_column]
		if gene in sample_dict:
			sys.exit("Gene "+gene+" appears at least twice in sample "+sample+" input file.")
  		else:
			sample_dict[gene] = float(expr)
  	for gene in gene_list:
  		if gene in sample_dict:
  	  		gene_dict[gene][0].append(sample_dict[gene])
    		else:
    			gene_dict[gene].append("NA")
  	sample_file.close()
input_files.close()

gene_count=0.0
for gene in gene_list:
	gene_count += 1
	print("Gene "+str(gene_count/len(gene_list)))
	expr_values = gene_dict[gene][0]
	bc_expr_values, bc_lambda = stats.boxcox( [ x + 1e-9 for x in expr_values ] ) #add 1e-9 for 0-valued expression levels
	bc_cdf = list( stats.norm.cdf( bc_expr_values ) )
	gene_dict[gene][1] = bc_expr_values
	gene_dict[gene][2] = bc_lambda
	gene_dict[gene][3] = bc_cdf
	pct75 = np.percentile(expr_values, 75)
	pct25 = np.percentile(expr_values, 25)
	iqr = pct75-pct25
	overexpression_outlier_level_list = []
  	underexpression_outlier_level_list = []
  	for i in range(outlier_levels):
    		overexpression_outlier_level_list.append(pct75+1.5*iqr*(i+1))
		underexpression_outlier_level_list.append(pct25-1.5*iqr*(i+1))
	gene_dict[gene][4] = overexpression_outlier_level_list
	gene_dict[gene][5] = underexpression_outlier_level_list

output_file.write( '\t'.join(['gene','sample','expression_level','percentile','boxcox_z', 'boxcox_l','boxcox_cdf'])+'\t'+'\t'.join(['overexpression'+str(x) for x in range(1,outlier_levels+1)])+'\t'+'\t'.join(['underexpression'+str(x) for x in range(1,outlier_levels+1)])+'\n')
for gene in gene_list:
	for sample_number in range(len(sample_list)):
		expr = gene_dict[gene][0][sample_number]
		pct = float( [ x <= gene_dict[gene][0][sample_number] for x in gene_dict[gene][0] ].count(True) )/len(sample_list)
		output_list = [gene, sample_list[sample_number], str(expr), str(pct), str(gene_dict[gene][1][sample_number]), str(gene_dict[gene][2]), str(gene_dict[gene][3][sample_number])]
    		for i in range(outlier_levels):
    			output_list.append(str(int(expr > gene_dict[gene][4][i])))
    		for i in range(outlier_levels):
    			output_list.append(str(int(expr < gene_dict[gene][5][i])))
    		output_file.write('\t'.join(output_list)+'\n')
output_file.close()
