import sys
import math
import numpy as np
from scipy import stats

#import warnings
#warnings.simplefilter('error')

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
  gene_index = first_line.index("gene")
  expression_index = first_line.index("expression_level")
  percentile_index = first_line.index("percentile")
  over_outlier_index = first_line.index("overexpression"+outlier_level)
  under_outlier_index = first_line.index("underexpression"+outlier_level)
  fusion_index = first_line.index("fusion_info")
  geneA_index = first_line.index("geneA")
  geneB_index = first_line.index("geneB")
  for line in file:
    line = line.strip().split()
    gene = line[gene_index]
    if line[expression_index] == "NA":
      pass
    else:
      expr = float(line[expression_index])
      over_outlier = int(line[over_outlier_index])
      under_outlier = int(line[under_outlier_index])
      fusion = line[fusion_index]
      pct = float(line[percentile_index])
      geneA_status = (gene in line[geneA_index].split(";"))
      geneB_status = (gene in line[geneB_index].split(";"))
      for gene_AB in [ gene+"__5prime", gene+"__3prime" ]:
        if gene_AB not in gene_dict:
          #gene_dict[gene] = [ 0.expression , 1.overexpression, 2.underexpression, 3.fusion, 4.t-test, 5.mann-whitney-u, 
          #                    6.fishers exact over, 7.fishers exact under, 8.num fusion over outliers, 9.num fusion under outliers, 10.total samples with fusion data
          #                    11. percentiles, 12. fusion percentile median , 13. number of fusions, 14 15 geneA gene B status]
          gene_dict[gene_AB] = [[None]*n_samples, [None]*n_samples, [None]*n_samples, [None]*n_samples, None, None, None, None, None, None, None, [None]*n_samples, [None]*n_samples, [None]*n_samples, [None]*n_samples, [None]*n_samples]
        gene_dict[gene_AB][0][i] = expr
        gene_dict[gene_AB][1][i] = over_outlier
        gene_dict[gene_AB][2][i] = under_outlier
        gene_dict[gene_AB][3][i] = fusion
        gene_dict[gene_AB][11][i] = pct
        gene_dict[gene_AB][14][i] = geneA_status
        gene_dict[gene_AB][15][i] = geneB_status
        if fusion == "Fusion_file_NA":
          gene_dict[gene_AB][3][i] = fusion
        elif gene_AB.endswith("__5prime") and geneA_status:
          gene_dict[gene_AB][3][i] = fusion
        elif gene_AB.endswith("__3prime") and geneB_status:
          gene_dict[gene_AB][3][i] = fusion
        else:
          gene_dict[gene_AB][3][i] = "None_reported"
  file.close()

for k,v in gene_dict.items():
  #t-test expression levels between fusion groups
  fusion_exp = []
  nonfusion_exp = []
  for a,b in zip(v[0], v[3]):
    if b == "None_reported":
      nonfusion_exp.append(a)
    elif b == "Fusion_file_NA":
      pass
    else:
      fusion_exp.append(a)
  if len(fusion_exp) < 2 or len(nonfusion_exp) < 2:
    t, gene_dict[k][4] = [float('NaN'), float('NaN')]
    mwu, gene_dict[k][5] = [float('NaN'), float('NaN')]
  else:
    t, gene_dict[k][4] = stats.ttest_ind(fusion_exp, nonfusion_exp, axis=0, equal_var=False)
    #Mann-Whitney U test
    mwu, gene_dict[k][5] = stats.mannwhitneyu(fusion_exp, nonfusion_exp, alternative='two-sided')
  #non-parametric test overexpression outlier status and fusion status
  fus_overout, fus_notoverout, notfus_overout, notfus_notoverout = [0]*4
  for a,b in zip(v[1],v[3]):
    if a == 0 and b == "None_reported":
      notfus_notoverout += 1
    elif a == 1 and b == "None_reported":
      notfus_overout += 1
    elif a == 0 and b != "Fusion_file_NA":
      fus_notoverout += 1
    elif a == 1 and b != "Fusion_file_NA":
      fus_overout += 1
  oddsratio, gene_dict[k][6] = stats.fisher_exact([[notfus_notoverout,notfus_overout],[fus_notoverout,fus_overout]])
  gene_dict[k][8] = fus_overout
  #non-parametric test underexpression outlier status and fusion status
  fus_underout, fus_notunderout, notfus_underout, notfus_notunderout = [0]*4
  for a,b in zip(v[2],v[3]):
    if a == 0 and b == "None_reported":
      notfus_notunderout += 1
    elif a == 1 and b == "None_reported":
      notfus_underout += 1
    elif a == 0 and b != "Fusion_file_NA":
      fus_notunderout += 1
    elif a == 1 and b != "Fusion_file_NA":
      fus_underout += 1
  oddsratio, gene_dict[k][7] = stats.fisher_exact([[notfus_notunderout,notfus_underout],[fus_notunderout,fus_underout]])
  gene_dict[k][9] = fus_underout
  gene_dict[k][10] = notfus_notunderout + notfus_underout + fus_notunderout + fus_underout
  #Fusion status and percertile
  fus_pct = []
  for a,b in zip(v[3],v[11]):
    if a == "None_reported":
      pass
    elif a != "Fusion_file_NA":
      fus_pct.append(b)
  if not fus_pct == []:
    gene_dict[k][12] = float(np.median(fus_pct))
  else:
    gene_dict[k][12] = float('NaN')
  gene_dict[k][13] = len(fus_pct)

multiple_test_corrected_pvalue = 0.05/len(gene_dict.keys())
for gene in sorted(gene_dict.keys()):
  print_gene=False
  v = gene_dict[gene]
  ttest, mwutest, overfisher, underfisher, fusionpct = ["NS"]*5
  #t-test
  if not math.isnan(v[4]) and v[4] < multiple_test_corrected_pvalue:
    ttest = v[4]
    print_gene=True
  #MWU test
  if not math.isnan(v[5]) and v[5] < multiple_test_corrected_pvalue:
    mwutest = v[5]
    print_gene=True
  #Fisher Exact -- overexpression
  if v[6] < multiple_test_corrected_pvalue:
    overfisher = v[6]
    print_gene=True
  #Fisher Exact -- underexpression
  if v[7] < multiple_test_corrected_pvalue:
    underfisher = v[7]
    print_gene=True
  #Fusion percentile
  if not math.isnan(v[12]) and (v[12] >= 0.90 or v[12] <= 0.10):
    fusionpct = v[12]
    print_gene=True

  n_fusions = v[13]
  if print_gene and n_fusions > 2:
    print('\t'.join([str(x) for x in [gene.split("__")[0], gene.split("__")[1], n_fusions, ttest, mwutest, overfisher, underfisher, fusionpct]]))

