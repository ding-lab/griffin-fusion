import sys
import numpy

ensg_gene = open(sys.argv[1], 'r')
abun_file = open(sys.argv[2], 'r')
output_file = open(sys.argv[3],'w')

#set up ensg dict
ensg_dict = {}
gene_dict = {}
gene_list = []
for line in ensg_gene:
  ensg, gene = line.strip().split()
  ensg_dict[ ensg ] = gene
  if gene not in gene_dict:
    gene_dict[gene] = [ ensg, 0 ]
  else:
    gene_dict[gene][0] += ','+ensg
  if gene not in gene_list:
    gene_list.append(gene)
ensg_gene.close()

#add abundances to ensg_dict
abun_file.readline() #skip first line
for line in abun_file:
  ensg = line.strip().split()[0].split('_')[0]
  tpm = float(line.strip().split()[4])
  if ensg in ensg_dict:
    gene = ensg_dict[ ensg ]
    gene_dict[gene][1] += tpm
abun_file.close()

#convert gene level tpm to log10(tpm+1)
for gene in gene_dict.keys():
  gene_dict[gene].append( numpy.log10(gene_dict[gene][1] + 1) )

#print output file
output_file.write( '\t'.join( ['ensg', 'gene', 'tpm','log10tpm'] ) + '\n')
for gene in gene_list:
  output_file.write( '\t'.join( [ gene_dict[gene][0], gene, str(gene_dict[gene][1]), str(gene_dict[gene][2]) ] )+ '\n')
output_file.close()
