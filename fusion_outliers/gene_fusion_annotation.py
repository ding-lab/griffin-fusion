import sys

def fix_ighl(gene):
  #From the file given by STAR-Fusion: gtf=~/GRCh37_gencode_v19_CTAT_lib_July192017/ctat_genome_lib_build_dir/ref_annot.gtf
  #grep -P "IG_V_gene|IG_C_gene|IG_D_gene|IG_J_gene" $gtf | grep "^chr14" | grep -w transcript | cut -f9 | cut -f2,4,10 -d' ' | tr '"; ' '\t' | cut -f2,6,10 | uniq | grep -vP "^ENSGR" | sed 's/\./\t/1' | sed 's/\./\t/1' | cut -f5 | uniq | sort -u > igh_superlocus_genes.txt
  #grep -P "IG_V_gene|IG_C_gene|IG_D_gene|IG_J_gene" $gtf | grep "^chr22" | grep -w transcript | cut -f9 | cut -f2,4,10 -d' ' | tr '"; ' '\t' | cut -f2,6,10 | uniq | grep -vP "^ENSGR" | sed 's/\./\t/1' | sed 's/\./\t/1' | cut -f5 | uniq | sort -u > igl_superlocus_genes.txt
  #These get special treatment in STAR-Fusion output. Fusions may be reported with partner gene IGH@, IGH-@, or IGL@ if they come from the IGH (or IGH- reverse complement) or IGL lists, or may be reported as the particular IGH or IGL gene.
  #The IGH super locus is indicated by IGH@ or IGH-@ (+strand)
  #The IGL super locus is indicated by IGL@ (no opposite strand)
  star_fusion_igh_superlocus=["IGH@","IGH-@","IGHA1","IGHA2","IGHD","IGHD1-1","IGHD1-14","IGHD1-20","IGHD1-26","IGHD1-7","IGHD2-15","IGHD2-2","IGHD2-21","IGHD2-8","IGHD3-10","IGHD3-16","IGHD3-22","IGHD3-3","IGHD3-9","IGHD4-11","IGHD4-17","IGHD4-23","IGHD4-4","IGHD5-12","IGHD5-18","IGHD5-24","IGHD5-5","IGHD6-13","IGHD6-19","IGHD6-25","IGHD6-6","IGHD7-27","IGHE","IGHG1","IGHG2","IGHG3","IGHG4","IGHJ1","IGHJ2","IGHJ3","IGHJ4","IGHJ5","IGHJ6","IGHM","IGHV1-18","IGHV1-2","IGHV1-24","IGHV1-3","IGHV1-45","IGHV1-46","IGHV1-58","IGHV1-69","IGHV1-8","IGHV2-26","IGHV2-5","IGHV2-70","IGHV3-11","IGHV3-13","IGHV3-15","IGHV3-16","IGHV3-20","IGHV3-21","IGHV3-23","IGHV3-30","IGHV3-33","IGHV3-35","IGHV3-38","IGHV3-43","IGHV3-48","IGHV3-49","IGHV3-53","IGHV3-64","IGHV3-66","IGHV3-7","IGHV3-72","IGHV3-73","IGHV3-74","IGHV3-9","IGHV4-28","IGHV4-31","IGHV4-34","IGHV4-39","IGHV4-4","IGHV4-59","IGHV4-61","IGHV5-51","IGHV6-1","IGHV7-81"]
  star_fusion_igl_superlocus=["IGL@","IGLC1","IGLC2","IGLC3","IGLC7","IGLJ1","IGLJ2","IGLJ3","IGLJ4","IGLJ5","IGLJ6","IGLJ7","IGLV1-36","IGLV1-40","IGLV1-44","IGLV1-47","IGLV1-50","IGLV1-51","IGLV10-54","IGLV11-55","IGLV2-11","IGLV2-14","IGLV2-18","IGLV2-23","IGLV2-33","IGLV2-8","IGLV3-1","IGLV3-10","IGLV3-12","IGLV3-16","IGLV3-19","IGLV3-21","IGLV3-22","IGLV3-25","IGLV3-27","IGLV3-32","IGLV3-9","IGLV4-3","IGLV4-60","IGLV4-69","IGLV5-37","IGLV5-45","IGLV5-48","IGLV5-52","IGLV6-57","IGLV7-43","IGLV7-46","IGLV8-61","IGLV9-49"]
  if gene in star_fusion_igh_superlocus:
    g = "IGH@"
  elif gene in star_fusion_igl_superlocus:
    g = "IGL@"
  else:
    g = gene
  return(g)

def parse_fusion_file(fusion_file_path, fusion_tool):
  fusion_file = open(fusion_file_path, 'r')
  fusion_dict = {}
  if fusion_tool == "star-fusion":
    fusion_file.readline() #skip the header
    for line in fusion_file:
      line = line.strip().split('\t')
      FusionName, JunctionReadCount, SpanningFragCount, SpliceType, LeftGene, LeftBreakpoint, RightGene, RightBreakpoint = line[0:8]
      geneA, geneB = FusionName.split("--")
      g1 = fix_ighl(geneA)
      g2 = fix_ighl(geneB)
      if g1 not in fusion_dict:
        fusion_dict[g1] = [ [] for x in range(10) ]
      if g2 not in fusion_dict:
        fusion_dict[g2] = [ [] for x in range(10) ]
      if g1 in fusion_dict:
        fusion_dict[g1][0].append(geneA)
        fusion_dict[g1][1].append(geneB)
        fusion_dict[g1][2].append(FusionName)
        fusion_dict[g1][3].append(JunctionReadCount)
        fusion_dict[g1][4].append(SpanningFragCount)
        fusion_dict[g1][5].append(LeftGene)
        fusion_dict[g1][6].append(LeftBreakpoint)
        fusion_dict[g1][7].append(RightGene)
        fusion_dict[g1][8].append(RightBreakpoint)
        fusion_dict[g1][9].append('#'.join(line))
      if g2 in fusion_dict:
        fusion_dict[g2][0].append(geneA)
        fusion_dict[g2][1].append(geneB)
        fusion_dict[g2][2].append(FusionName)
        fusion_dict[g2][3].append(JunctionReadCount)
        fusion_dict[g2][4].append(SpanningFragCount)
        fusion_dict[g2][5].append(LeftGene)
        fusion_dict[g2][6].append(LeftBreakpoint)
        fusion_dict[g2][7].append(RightGene)
        fusion_dict[g2][8].append(RightBreakpoint)
        fusion_dict[g2][9].append('#'.join(line))
  else:
    fusion_file.close()
    sys.exit("Fusion tool "+tool+" not yet supported.")
  fusion_file.close()
  return(fusion_dict)

outlier_file_path=sys.argv[1]
fusion_outlier_file_path=sys.argv[2]
gene_list_file=open(sys.argv[3],'r')
fusion_tool=sys.argv[4]
fusion_file_path=sys.argv[5]

gene_list=[]
for gene in gene_list_file:
  gene_list.append(gene.strip())
gene_list_file.close()

fusion_dict = parse_fusion_file(fusion_file_path, fusion_tool)
fusion_outlier_file = open(fusion_outlier_file_path, 'w')
outlier_file = open(outlier_file_path, 'r')
line = outlier_file.readline()
line = line.strip().split()
line.extend(['fusion_tool','geneA','geneB','fusion','junction_read_count','spanning_read_count','left_gene','left_bp','right_gene','right_bp','fusion_info'])
fusion_outlier_file.write('\t'.join(line)+'\n')
for line in outlier_file:
  line = line.strip().split()
  line.append(fusion_tool)
  gene = fix_ighl(line[0])
  if gene in fusion_dict:
    for element in fusion_dict[gene]:
      line.append(';'.join(element))
  else:
    line.append(fusion_tool)
    line.extend(["NA"]*10)
  fusion_outlier_file.write('\t'.join(line)+'\n')

outlier_file.close()
fusion_outlier_file.close()

