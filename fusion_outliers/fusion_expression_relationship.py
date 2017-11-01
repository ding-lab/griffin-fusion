import sys
import math
from scipy import stats

import warnings
warnings.simplefilter('error')

input_files = open(sys.argv[1],'r')
outlier_level = sys.argv[2]

input_file_list = []
for input_file in input_files:
  input_file_list.append(input_file.strip())
input_files.close()
n_samples = len(input_file_list)

gene_dict = {}
for i in range(n_samples):
  file=open(input_file_list[i].strip(),'r')
  first_line = file.readline().strip().split()
  outlier_index = first_line.index("overexpression"+outlier_level)
  fusion_index = first_line.index("fusion_info")
  for line in file:
    line = line.strip().split()
    gene = line[0]
    if line[2] == "NA":
      pass
    else:
      expr = float(line[2])
      outlier = int(line[outlier_index])
      fusion = line[fusion_index]
      if gene not in gene_dict:
        gene_dict[gene] = [[None]*n_samples, [None]*n_samples, [None]*n_samples, None, None, None, None]
      gene_dict[gene][0][i] = expr
      gene_dict[gene][1][i] = outlier
      gene_dict[gene][2][i] = fusion
  file.close()

for k,v in gene_dict.items():
  #t-test expression levels between fusion groups
  fusion_exp = []
  nonfusion_exp = []
  for a,b in zip(v[0], v[2]):
    if b == "None_reported":
      nonfusion_exp.append(a)
    elif b == "Fusion_file_NA":
      pass
    else:
      fusion_exp.append(a)
  if len(fusion_exp) < 2 or len(nonfusion_exp) < 2:
    t, gene_dict[k][3] = [float('NaN'), float('NaN')]
    mwu, gene_dict[k][4] = [float('NaN'), float('NaN')]
  else:
    t, gene_dict[k][3] = stats.ttest_ind(fusion_exp, nonfusion_exp, axis=0, equal_var=False)
    #Mann-Whitney U test
    mwu, gene_dict[k][4] = stats.mannwhitneyu(fusion_exp, nonfusion_exp, alternative='two-sided')
  #non-parametric test outlier status and fusion status
  fus_out, fus_notout, notfus_out, notfus_notout = [0]*4
  for a,b in zip(v[1],v[2]):
    if a == 0 and b == "None_reported":
      notfus_notout += 1
    elif a == 1 and b == "None_reported":
      notfus_out += 1
    elif a == 0 and b != "Fusion_file_NA":
      fus_notout += 1
    elif a == 1 and b != "Fusion_file_NA":
      fus_out += 1
  oddsratio, gene_dict[k][5] = stats.fisher_exact([[notfus_notout,notfus_out],[fus_notout,fus_out]])
  gene_dict[k][6] = fus_out

sig_genes = []
for k,v in gene_dict.items():
  if (not math.isnan(v[3]) and v[3] < 0.05) or (not math.isnan(v[4]) and v[4] < 0.05) or (v[5] < 0.05 and gene_dict[k][6] > 0):
    sig_genes.append(k)

for gene in sorted(sig_genes):
  print('\t'.join([gene, str(gene_dict[gene][3]), str(gene_dict[gene][4]), str(gene_dict[gene][5]), str(gene_dict[gene][6])]))

