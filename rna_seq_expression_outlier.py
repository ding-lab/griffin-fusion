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

gene_dict = {}
for gene in gene_list:
	if gene in gene_dict:
		next
	else:
		gene_dict[gene] = [[],None,None]

sample_list = []

for line in input_files:
        print(line)
	line = line.strip().split()
	sample = line[0]
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
  		gene = sline[1]
  		expr = sline[3]
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

for gene in gene_list:
        print(gene)
	expr_values = gene_dict[gene][0]
	print(expr_values)
	#bc_expr_values, bc_lambda = stats.boxcox(expr_values)
	bc_expr_values = ["NA" for x in range(len(expr_values))]
	bc_lambda = "NA"
	gene_dict[gene][1] = bc_expr_values
	gene_dict[gene][2] = bc_lambda
	pct75 = np.percentile(expr_values, 75)
	pct25 = np.percentile(expr_values, 25)
	iqr = pct75-pct25
	overexpression_outlier_level_list = []
  	underexpression_outlier_level_list = []
  	for i in range(outlier_levels):
    		overexpression_outlier_level_list.append(pct75+1.5*iqr*(i+1))
    		underexpression_outlier_level_list.append(pct25-1.5*iqr*(i+1))

output_file.write( '\t'.join(['gene','sample','expression_level','boxcox_pct'])+'\t'+'\t'.join(['overexpression'+str(x) for x in range(1,outlier_levels+1)])+'\t'+'\t'.join(['underexpression'+str(x) for x in range(1,outlier_levels+1)])+'\n')
for gene in gene_list:
	for sample_number in range(len(sample_list)):
		expr = str(gene_dict[gene][0][sample_number])
		#bc_pct = float(list(gene_dict[gene][1] <= gene_dict[gene][1][sample_number]).count(True))/len(sample_list)
    		bc_pct = "NA"
		output_list = [gene, sample_list[sample_number], expr, bc_pct]
    		for i in range(outlier_levels):
    			output_list.append(str(int(expr > overexpression_outlier_level_list[i])))
    		for i in range(outlier_levels):
    			output_list.append(str(int(expr < underexpression_outlier_level_list[i])))
		print(output_list)
    		output_file.write('\t'.join(output_list)+'\n')
output_file.close()
