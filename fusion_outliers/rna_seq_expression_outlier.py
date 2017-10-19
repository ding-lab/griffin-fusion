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

sample_of_interest=sys.argv[3]
output_file=open(sys.argv[4],'w')
outlier_levels=int(sys.argv[5])
input_file_header=sys.argv[6]
gene_column=int(sys.argv[7])
expr_column=int(sys.argv[8])

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
        gene_dict[gene][0].append("NA")
  sample_file.close()
input_files.close()

#check that sample of interest appears in input files
if not sample_of_interest in sample_list:
  sys.exit("Sample of interest not found in list of input files")

for gene in gene_list:
  if gene in gene_dict:
    expr_values = gene_dict[gene][0]
    if all([ x == "NA" for x in expr_values]):
      gene_dict[gene][1] = expr_values
      gene_dict[gene][2] = "NA"
      gene_dict[gene][3] = expr_values
      gene_dict[gene][4] = "NA"
      gene_dict[gene][5] = "NA"
    else:
      if all([ x == 0 for x in expr_values]):
        bc_expr_values = expr_values
        bc_lambda = "NA"
      else:
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
sample_number=sample_list.index(sample_of_interest)
for gene in sorted(gene_dict.keys()):
  gene_na = all([x == "NA" for x in gene_dict[gene][0]])
  expr = gene_dict[gene][0][sample_number]
  if gene_na:
    pct = "NA"
  else:
    pct = float( [ x <= expr for x in gene_dict[gene][0] ].count(True) )/len(sample_list)
  output_list = [gene, sample_list[sample_number], str(expr), str(pct), str(gene_dict[gene][1][sample_number]), str(gene_dict[gene][2]), str(gene_dict[gene][3][sample_number])]
  if gene_na:
    output_list.extend(["NA"]*outlier_levels*2)
  else:
    for i in range(outlier_levels):
      output_list.append(str(int(expr > gene_dict[gene][4][i])))
    for i in range(outlier_levels):
      output_list.append(str(int(expr < gene_dict[gene][5][i])))
  output_file.write('\t'.join(output_list)+'\n')
output_file.close()
