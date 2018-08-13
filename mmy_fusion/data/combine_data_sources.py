regenerate_expr_file = False #warning: takes 3+ hours if regenerate_expr_file = True and test = False
test = False

import sys
import numpy as np

#gene lists
f = open("oncogene.tsv", "r")
oncogene_list = []
for line in f:
  oncogene_list.append(line.strip())
f.close()
f = open("tsg.tsv", "r")
tsg_list = []
for line in f:
  tsg_list.append(line.strip())
f.close()
f = open("kinase.tsv", "r")
kinase_list = []
for line in f:
  kinase_list.append(line.strip())
f.close()
f = open("mmy_known.tsv", "r")
mmy_known_list = []
for line in f:
  mmy_known_list.append(line.strip())
f.close()
f = open("driver.tsv", "r")
driver_list = []
for line in f:
  driver_list.append(line.strip())
f.close()

#Druggability DEPO
f = open("DEPO_final_20170206.txt")
f.readline()
depo_dict = {}
for line in f:
  Gene, Disease, Variant, Effect, Source, Drug, Drug_Class, Evidence_Level, PubMed_ID, Gene_Class, Variant_Type = line.strip().split('\t')
  Gene.replace("-","--")
  if Effect == "Sensitive" and Variant == "any" and Variant_Type == "Fusion":
    if Gene not in depo_dict:
      depo_dict[Gene] = ''
f.close()

#read in sample list, assign primary secondary etc
f = open("sample_list.with_file_names.txt","r")
w = open("sample_list.primary.txt","w")
sample_dict = {}
srr_dict = {}
for line in f:
  mmrf, primary_srr, other_srr, sf, sf2, es, es2, i, i2, pd, pd2, fc, fc2, wgstn_srr, delly, manta, lumpy, cnv, tpm = line.strip().split()
  w.write(mmrf+"\t"+primary_srr+"\n")
  sample_dict[mmrf] = [primary_srr]
  srr_dict[primary_srr] = 1
  if other_srr != "NA":
    count = 1
    for srr in other_srr.split(","):
      count += 1
      sample_dict[mmrf].append(srr)
      srr_dict[srr] = count
f.close()
w.close()

#read in Filtered_Fusions.tsv and create filtered_fusions_dict
genes_with_fusions = {}
fusion_recurrence = {}
f = open("Filtered_Fusions_v2.tsv","r")
filtered_fusions_dict = {}
f.readline()
for line in f:
  FusionName, LeftBreakpoint, RightBreakpoint, Cancer, Sample, JunctionReadCount, SpanningFragCount, FFPM, PROT_FUSION_TYPE, GTEx, Callers, CallerN = line.strip().split()
  chrA, posA, strandA = LeftBreakpoint.split(":")
  chrB, posB, strandB = RightBreakpoint.split(":")
  fusion_key = Cancer+":"+Sample+":"+FusionName
  geneA, geneB = FusionName.split("--")
  called_by = [ str(int(x in Callers)) for x in ["E","F","I","P","S"] ] # determine what callers called this fusion Ericscript FusionCatcher Integrate Prada Star-Fusion
  if geneA not in genes_with_fusions:
    genes_with_fusions[geneA] = ''
  if geneB not in genes_with_fusions:
    genes_with_fusions[geneB] = ''
  if fusion_key in filtered_fusions_dict:
    sys.exit(fusion_key + " already in filtered_fusions_dict")
  else:
    filtered_fusions_dict[fusion_key] = [fusion_key, FusionName, geneA, geneB, LeftBreakpoint, RightBreakpoint, chrA, posA, strandA, chrB, posB, strandB, Cancer, Sample, JunctionReadCount, SpanningFragCount, FFPM, PROT_FUSION_TYPE, Callers, CallerN, called_by[0], called_by[1], called_by[2], called_by[3], called_by[4]]
    if FusionName in fusion_recurrence:
      fusion_recurrence[FusionName] += 1
    else:
      fusion_recurrence[FusionName] = 1
f.close()

#read in SeqFISH.20180104.csv and create seqfish_dict
f = open("SeqFISH.20180104.csv","r")
w = open("seqfish_df.tsv","w")
w.write("\t".join(["MMRF","Study_Visit_ID","seqfish_category","has_event"])+"\n")
seqfish_dict = {}
f.readline()
for line in f:
  Study_Visit_ID, CN_del_13q14, CN_del_13q34, CN_del_17p13, CN_gain_1q21, Hyperdiploidy, Translocation_WHSC1_4_14, Translocation_CCND3_6_14, Translocation_MYC_8_14, Translocation_MAFA_8_14, Translocation_CCND1_11_14, Translocation_CCND2_12_14, Translocation_MAF_14_16, Translocation_MAFB_14_20 = line.strip().split(",")
  sample_key = Study_Visit_ID[0:9]
  if sample_key in seqfish_dict:
    #seqfish_dict[sample_key].append([sample_key, Study_Visit_ID, CN_del_13q14, CN_del_13q34, CN_del_17p13, CN_gain_1q21, Hyperdiploidy, Translocation_WHSC1_4_14, Translocation_CCND3_6_14, Translocation_MYC_8_14, Translocation_MAFA_8_14, Translocation_CCND1_11_14, Translocation_CCND2_12_14, Translocation_MAF_14_16, Translocation_MAFB_14_20])
    seqfish_dict[sample_key].append([x if x else "NA" for x in [sample_key, Study_Visit_ID, CN_del_13q14, CN_del_13q34, CN_del_17p13, CN_gain_1q21, Hyperdiploidy, Translocation_WHSC1_4_14, Translocation_CCND3_6_14, Translocation_MYC_8_14, Translocation_MAFA_8_14, Translocation_CCND1_11_14, Translocation_CCND2_12_14, Translocation_MAF_14_16, Translocation_MAFB_14_20]])
  else:
    seqfish_dict[sample_key] = [ [x if x else "NA" for x in [sample_key, Study_Visit_ID, CN_del_13q14, CN_del_13q34, CN_del_17p13, CN_gain_1q21, Hyperdiploidy, Translocation_WHSC1_4_14, Translocation_CCND3_6_14, Translocation_MYC_8_14, Translocation_MAFA_8_14, Translocation_CCND1_11_14, Translocation_CCND2_12_14, Translocation_MAF_14_16, Translocation_MAFB_14_20]]]
  list_of_indicators = ["sample_key", "Study_Visit_ID", "CN_del_13q14", "CN_del_13q34", "CN_del_17p13", "CN_gain_1q21", "Hyperdiploidy", "Translocation_WHSC1_4_14", "Translocation_CCND3_6_14", "Translocation_MYC_8_14", "Translocation_MAFA_8_14", "Translocation_CCND1_11_14", "Translocation_CCND2_12_14", "Translocation_MAF_14_16", "Translocation_MAFB_14_20"]
  for c in range(2,len(list_of_indicators)):
    w.write('\t'.join([sample_key, Study_Visit_ID, list_of_indicators[c], seqfish_dict[sample_key][-1][c]])+'\n')
f.close()
w.close()

#read in clinical_data.20171201.csv and create clinical_dict 
f = open("clinical_data.20171201.csv","r")
w = open("clinical_df.tsv","w")
w.write("\t".join(["MMRF", "Clinical_category", "Clinical_value"])+"\n")
clinical_dict = {}
f.readline()
for line in f:
  if line.strip(): #necessary to deal with empty line at end of input file
    PUBLIC_ID, Age, age_ge_66, Female, Race_White, Race_Black, Race_Other, race, ECOG, BM_Plasma_Cell_Percent, ISS_Stage, LDH, Bone_lesions, Plamacytoma, D_PT_deathdy, D_PT_lstalive, D_PT_pddy, TTPD, EFS, EFS_censor = line.strip().split(",")
    sample_key = PUBLIC_ID
    if sample_key in clinical_dict:
      sys.exit(sample_key + " already in clinical_dict")
    else:
      clinical_dict[sample_key] = [ str(x).replace(" ","_") if x else "NA" for x in [sample_key, PUBLIC_ID, Age, age_ge_66, Female, Race_White, Race_Black, Race_Other, race, ECOG, BM_Plasma_Cell_Percent, ISS_Stage, LDH, Bone_lesions, Plamacytoma, D_PT_deathdy, D_PT_lstalive, D_PT_pddy, TTPD, EFS, EFS_censor] ]
    list_of_indicators = ["sample_key", "PUBLIC_ID", "Age", "age_ge_66", "Female", "Race_White", "Race_Black", "Race_Other", "race", "ECOG", "BM_Plasma_Cell_Percent", "ISS_Stage", "LDH", "Bone_lesions", "Plamacytoma", "D_PT_deathdy", "D_PT_lstalive", "D_PT_pddy", "TTPD", "EFS", "EFS_censor"]
    for c in range(2, len(list_of_indicators)):
      w.write("\t".join([sample_key, list_of_indicators[c], clinical_dict[sample_key][c] ])+"\n")
f.close()
w.close()

#read in fusion_evidence_discordant_reads.100000.txt and create discordant_dict
f = open("fusion_evidence_discordant_reads.100000.txt","r")
discordant_dict = {}
f.readline()
for line in f:
  FusionName, LeftBreakpoint, RightBreakpoint, Cancer, Sample, JunctionReadCount, SpanningFragCount, FFPM, PROT_FUSION_TYPE, GTEx, CallerN, Callers, Overlap, bpRangeA, bpRangeB, depthA, depthB, n_discordant, discordant_reads, wgs_bam = line.strip().split()
  fusion_key = Cancer+":"+Sample+":"+FusionName
  if fusion_key in discordant_dict:
    sys.exit(fusion_key + " already in discordant_dict")
  else:
    discordant_dict[fusion_key] = [fusion_key, Overlap, bpRangeA, bpRangeB, depthA, depthB, n_discordant, discordant_reads, wgs_bam]
f.close()

#read in Filtered_Fusions_100000_delly_manta_20180813.txt and create sv_dict
f = open("Filtered_Fusions_100000_delly_manta_20180813.txt","r")
sv_dict = {}
f.readline()
for line in f:
  FusionName, LeftBreakpoint, RightBreakpoint, Cancer, Sample, JunctionReadCount, SpanningFragCount, FFPM, PROT_FUSION_TYPE, GTEx, Callers, CallerNumber, sv_type, delly_support_level, delly_evidence, manta_support_level, manta_evidence = line.strip().split()
  fusion_key = Cancer+":"+Sample+":"+FusionName
  # TODO need to fix multiple time point entries and change NAs to None
  #if delly_evidence == "NA":
  #  delly_evidence = "None"
  #if manta_evidence == "NA":
  #  manta_evidence == "None"
  if fusion_key in sv_dict:
    sys.exit(fusion_key + " already in sv_dict")
  else:
    sv_dict[fusion_key] = [fusion_key, sv_type, delly_support_level, delly_evidence, manta_support_level, manta_evidence]
f.close()

if regenerate_expr_file:
  #read in combined_cnv_results.txt and create multi-layered dictionary cnv_dict
  if test:
    f = open("combined_cnv_results.test.txt","r")
  else:
    f = open("combined_cnv_results.txt","r")
  cnv_dict = {}
  f.readline()
  for line in f:
    if line.strip().split()[6] == "0":
      continue
    mmrf, chr_bed, start_bed, stop_bed, ensg_bed, gene_bed, chr_cnv, start_cnv, stop_cnv, log2ratio_cnv, bp_overlap = line.strip().split()
    if mmrf not in cnv_dict:
      cnv_dict[mmrf] = {}
    if gene_bed not in cnv_dict[mmrf]:
      cnv_dict[mmrf][gene_bed] = []
    print(mmrf, gene_bed)
    cnv_dict[mmrf][gene_bed].append([mmrf, chr_bed, start_bed, stop_bed, ensg_bed, gene_bed, chr_cnv, start_cnv, stop_cnv, log2ratio_cnv, bp_overlap])
  f.close()

  #read in mmy_gene_tpm_table.tsv and create expr_dict
  if test:
    f = open("mmy_gene_tpm_table.test.tsv","r")
  else:
    f = open("mmy_gene_tpm_table.tsv","r")
  gene_expr_dict = {}
  f.readline()
  for line in f:
    mmrf, srr, ensg, gene, tpm, log10tpm = line.strip().split()
    if srr in srr_dict and srr_dict[srr] == 1:
      if gene in gene_expr_dict:
        gene_expr_dict[gene].append(float(tpm))
      else:
        gene_expr_dict[gene] = [float(tpm)]
  f.close()

  #re-read mmy_gene_tpm_table.tsv and calculate percentiles at fusion genes
  if test:
    f = open("mmy_gene_tpm_table.test.tsv","r")
    w = open("mmy_gene_expr_with_fusions.test.tsv","w")
  else:
    f = open("mmy_gene_tpm_table.tsv","r")
    w = open("mmy_gene_expr_with_fusions.tsv","w")
  w.write('\t'.join(["expr_key", "mmrf", "srr", "ensg", "gene", "tpm", "log10tpm", "pct", "pct75_tpm", "pct25_tpm", "iqr_tpm", "pct75_log10tpm", "pct25_log10tpm", "iqr_log10tpm", "outlier_over_tpm", "outlier_under_tpm", "outlier_over_log10tpm", "outlier_under_log10tpm", "gene_avg_cnv"])+"\n")
  total_expr_lines = len(srr_dict.keys())*len(genes_with_fusions.keys())
  count_up = 0
  expr_dict = {}
  f.readline()
  for line in f:
    mmrf, srr, ensg, gene, tpm, log10tpm = line.strip().split()
    if gene in genes_with_fusions:
      count_up += 1
      print(1.0*count_up/total_expr_lines)
      expr_key = mmrf+":"+srr+":"+gene
      if expr_key in expr_dict:
        sys.exit(expr_key + " already in expr_dict")
      else:
        gene_expr_list = gene_expr_dict[gene]
        pct = 1.0*sum([float(tpm) >= float(x) for x in gene_expr_list])/len(gene_expr_list)
        pct75_tpm = np.percentile(gene_expr_list, 75)
        pct25_tpm = np.percentile(gene_expr_list, 25)
        iqr_tpm = pct75_tpm-pct25_tpm
        pct75_log10tpm = np.percentile([ np.log10(x+1) for x in gene_expr_list ], 75)
        pct25_log10tpm = np.percentile([ np.log10(x+1) for x in gene_expr_list ], 25)
        iqr_log10tpm = pct75_log10tpm-pct25_log10tpm
        outlier_over_tpm = int(float(tpm) >= (float(pct75_tpm) + 1.5*float(iqr_tpm)))
        outlier_under_tpm = int(float(tpm) <= (float(pct25_tpm) - 1.5*float(iqr_tpm)))
        outlier_over_log10tpm = int(float(log10tpm) >= (float(pct75_log10tpm) + 1.5*float(iqr_log10tpm)))
        outlier_under_log10tpm = int(float(log10tpm) <= (float(pct25_log10tpm) - 1.5*float(iqr_log10tpm)))
        #add info from combined_cnv_results.txt
        if mmrf in cnv_dict and gene in cnv_dict[mmrf]:
          n_obs = len(cnv_dict[mmrf][gene])
          gene_cnv = []
          gene_bp = []
          num = 0
          den = 0
          for i in range(n_obs):
            if cnv_dict[mmrf][gene][i][9] != ".":
              gene_cnv.append(cnv_dict[mmrf][gene][i][9])  #.append([mmrf, chr_bed, start_bed, stop_bed, ensg_bed, gene_bed, chr_cnv, start_cnv, stop_cnv, log2ratio_cnv, bp_overlap])
              gene_bp.append(cnv_dict[mmrf][gene][i][10])
              num += float(cnv_dict[mmrf][gene][i][9])*float(cnv_dict[mmrf][gene][i][10])
              den += float(cnv_dict[mmrf][gene][i][10])
          if den > 0:
            gene_avg_cnv = 1.0*num/den
          else:
            gene_avg_cnv = "NA"
        else:
          gene_avg_cnv = "NA"
        expr_dict[expr_key] = [expr_key, mmrf, srr, ensg, gene, tpm, log10tpm, pct, pct75_tpm, pct25_tpm, iqr_tpm, pct75_log10tpm, pct25_log10tpm, iqr_log10tpm, outlier_over_tpm, outlier_under_tpm, outlier_over_log10tpm, outlier_under_log10tpm, gene_avg_cnv]        
        w.write('\t'.join([str(x) for x in expr_dict[expr_key]])+"\n")

  w.close()
  f.close()

#Generate expr_dict based on mmy_gene_expr_with_fusions.tsv
if test:
  f = open("mmy_gene_expr_with_fusions.test.tsv", "r")
else:
  f = open("mmy_gene_expr_with_fusions.tsv","r")
f.readline()
expr_dict = {}
for line in f:
  expr_key, mmrf, srr, ensg, gene, tpm, log10tpm, pct, pct75_tpm, pct25_tpm, iqr_tpm, pct75_log10tpm, pct25_log10tpm, iqr_log10tpm, outlier_over_tpm, outlier_under_tpm, outlier_over_log10tpm, outlier_under_log10tpm, gene_avg_cnv = line.strip().split()
  expr_dict[expr_key] = [expr_key, mmrf, srr, ensg, gene, tpm, log10tpm, pct, pct75_tpm, pct25_tpm, iqr_tpm, pct75_log10tpm, pct25_log10tpm, iqr_log10tpm, outlier_over_tpm, outlier_under_tpm, outlier_over_log10tpm, outlier_under_log10tpm, gene_avg_cnv]
f.close()

#write to files
if test:
  w = open("fusion_df.test.txt","w")
else:
  w = open("fusion_df.txt","w")

column_labels = ["mmrf"]
column_labels.extend(["srr"])
column_labels.extend(["fusion"])
column_labels.extend(["sample_number"])
column_labels.extend(["has_secondary"])
column_labels.extend(["geneA", "geneB", "LeftBreakpoint", "RightBreakpoint",  "chrA", "posA", "strandA", "chrB", "posB", "strandB", "JunctionReadCount", "SpanningFragCount", "FFPM", "PROT_FUSION_TYPE", "Callers", "CallerN", "called_by_Eriscript", "called_by_FusionCatcher", "called_by_Integrate", "called_by_Prada", "called_by_StarFusion"])
column_labels.extend(["seqfish_Study_Visit_ID", "seqfish_CN_del_13q14", "seqfish_CN_del_13q34", "seqfish_CN_del_17p13", "seqfish_CN_gain_1q21", "seqfish_Hyperdiploidy", "seqfish_Translocation_WHSC1_4_14", "seqfish_Translocation_CCND3_6_14", "seqfish_Translocation_MYC_8_14", "seqfish_Translocation_MAFA_8_14", "seqfish_Translocation_CCND1_11_14", "seqfish_Translocation_CCND2_12_14", "seqfish_Translocation_MAF_14_16", "seqfish_Translocation_MAFB_14_20"])
column_labels.extend(["Age", "age_ge_66", "Female", "Race_White", "Race_Black", "Race_Other", "race", "ECOG", "BM_Plasma_Cell_Percent", "ISS_Stage", "LDH", "Bone_lesions", "Plamacytoma", "D_PT_deathdy", "D_PT_lstalive", "D_PT_pddy", "TTPD", "EFS", "EFS_censor"])
column_labels.extend(["R-ISS"])
column_labels.extend(["Overlap", "bpRangeA", "bpRangeB", "depthA", "depthB", "n_discordant", "discordant_reads","wgs_bam"])
column_labels.extend(["sv_type", "delly_support_level", "delly_evidence", "manta_support_level", "manta_evidence"])
column_labels.extend(["geneA_tpm", "geneA_log10tpm", "geneA_pct", "geneA_pct75_tpm", "geneA_pct25_tpm", "geneA_iqr_tpm", "geneA_pct75_log10tpm", "geneA_pct25_log10tpm", "geneA_iqr_log10tpm", "geneA_outlier_over_tpm", "geneA_outlier_under_tpm", "geneA_outlier_over_log10tpm", "geneA_outlier_under_log10tpm", "geneA_log2ratio_cnv"])
column_labels.extend(["geneB_tpm", "geneB_log10tpm", "geneB_pct", "geneB_pct75_tpm", "geneB_pct25_tpm", "geneB_iqr_tpm", "geneB_pct75_log10tpm", "geneB_pct25_log10tpm", "geneB_iqr_log10tpm", "geneB_outlier_over_tpm", "geneB_outlier_under_tpm", "geneB_outlier_over_log10tpm", "geneB_outlier_under_log10tpm", "geneB_log2ratio_cnv"])
column_labels.extend(["geneA_oncogene", "geneA_tsg", "geneA_kinase", "geneA_mmy_known", "geneA_driver", "geneB_oncogene", "geneB_tsg", "geneB_kinase", "geneB_mmy_known", "geneB_driver", "fusion_recurrence"])
column_labels.extend(["drug_fusion", "drug_geneA", "drug_geneB"])
#column_labels.extend([""])


w.write('\t'.join([str(x) for x in column_labels])+"\n")
for fusion_key in sorted(filtered_fusions_dict.keys()):
  print_list = []
  mmrf, srr, fus = fusion_key.split(":")
  geneA, geneB = fus.split("--")
  print_list.append(mmrf) #mmrf
  print_list.append(srr) #srr
  print_list.append(fus) #geneA--geneB
  print_list.append(str(srr_dict[srr])) #1 2 3 4 sample
  if len(sample_dict[mmrf])>1:
    print_list.append("1") #has secondary
  else:
    print_list.append("0")
  print_list.extend([ filtered_fusions_dict[fusion_key][i] for i in [2,3,4,5,6,7,8,9,10,11,14,15,16,17,18,19,20,21,22,23,24] ]) #geneA, geneB, LeftBreakpoint, RightBreakpoint, chrA, posA, strandA, chrB, posB, strandB, JunctionReadCount, SpanningFragCount, FFPM, PROT_FUSION_TYPE, Callers, CallerN, called_by_X
  #seqfish
  if srr_dict[srr] == 1 and mmrf in seqfish_dict:
    print_list.extend( seqfish_dict[mmrf][0][1:] ) #Study_Visit_ID, CN_del_13q14, CN_del_13q34, CN_del_17p13, CN_gain_1q21, Hyperdiploidy, Translocation_WHSC1_4_14, Translocation_CCND3_6_14, Translocation_MYC_8_14, Translocation_MAFA_8_14, Translocation_CCND1_11_14, Translocation_CCND2_12_14, Translocation_MAF_14_16, Translocation_MAFB_14_20
  else: #Study_Visit_ID can end in PB or BM, both are cancer samples
    print_list.extend( ["NA"]*14 )
  
  #clinical
  print_list.extend( clinical_dict[mmrf][2:] ) #Age, age_ge_66, Female, Race_White, Race_Black, Race_Other, race, ECOG, BM_Plasma_Cell_Percent, ISS_Stage, LDH, Bone_lesions, Plamacytoma, D_PT_deathdy, D_PT_lstalive, D_PT_pddy, TTPD, EFS, EFS_censor
  #R-ISS
  if mmrf in seqfish_dict and mmrf in clinical_dict and srr_dict[srr] == 1: # only report for primary samples
    iss_stage = clinical_dict[mmrf][11]
    if clinical_dict[mmrf][12] != "NA":
      ldh = float(clinical_dict[mmrf][12])
    else:
      ldh = "NA"
    del_17p = seqfish_dict[mmrf][0][4]
    t414 = seqfish_dict[mmrf][0][7]
    t1416 = seqfish_dict[mmrf][0][13]
    # R-ISS Stage 1
    if iss_stage == "NA" or ldh == "NA" or del_17p == "NA" or t414 == "NA" or t1416 == "NA":
      riss = "NA"
    elif iss_stage == "1" and ldh < 300 and del_17p == "0" and t414 == "0" and t1416 == "0":
      riss = "1"
    elif iss_stage == "3" and ( ldh > 300 or del_17p == "1" or t414 == "1" or t1416 == "1" ):
      riss = "3"
    else:
      riss = "2"
    print_list.extend([riss])
  else:
    print_list.extend(["NA"])
  #discordant read validation
  if fusion_key in discordant_dict and srr_dict[srr] == 1: # only report for primary samples
    print_list.extend( discordant_dict[ fusion_key ][1:] ) #Overlap, bpRangeA, bpRangeB, depthA, depthB, n_discordant, discordant_reads, wgs_bam
  else:
    print_list.extend( ["NA"]*8 )
  # SV caller validation
  if fusion_key in sv_dict and srr_dict[srr] == 1: # only report for primary samples
    print_list.extend( sv_dict[ fusion_key ][1:] ) # sv_type, delly_support_level, delly_evidence, manta_support_level, manta_evidence
  else:
    print_list.extend( ["NA"]*5 )
  #expression geneA and geneB
  if mmrf+":"+srr+":"+geneA in expr_dict:
    print_list.extend( expr_dict[ mmrf+":"+srr+":"+geneA ][5:] ) #tpm, log10tpm, pct, pct75_tpm, pct25_tpm, iqr_tpm, pct75_log10tpm, pct25_log10tpm, iqr_log10tpm, outlier_over_tpm, outlier_under_tpm, outlier_over_log10tpm, outlier_under_log10tpm, gene_avg_cnv
  else:
    print_list.extend( ["NA"]*14)
  if mmrf+":"+srr+":"+geneB in expr_dict:
    print_list.extend( expr_dict[ mmrf+":"+srr+":"+geneB ][5:] ) #tpm, log10tpm, pct, pct75_tpm, pct25_tpm, iqr_tpm, pct75_log10tpm, pct25_log10tpm, iqr_log10tpm, outlier_over_tpm, outlier_under_tpm, outlier_over_log10tpm, outlier_under_log10tpm, gene_avg_cnv
  else:
    print_list.extend( ["NA"]*14)
  #gene lists
  if geneA in oncogene_list:
    geneA_oncogene = 1
  else:
    geneA_oncogene = 0
  if geneA in tsg_list:
    geneA_tsg = 1
  else:
    geneA_tsg = 0
  if geneA in kinase_list:
    geneA_kinase = 1
  else:
    geneA_kinase = 0
  if geneA in mmy_known_list:
    geneA_mmy_known = 1
  else:
    geneA_mmy_known = 0
  if geneA in driver_list:
    geneA_driver = 1
  else:
    geneA_driver = 0
  if geneB in oncogene_list:
    geneB_oncogene = 1
  else:
    geneB_oncogene = 0
  if geneB in tsg_list:
    geneB_tsg = 1
  else:
    geneB_tsg = 0
  if geneB in kinase_list:
    geneB_kinase = 1
  else:
    geneB_kinase = 0
  if geneB in mmy_known_list:
    geneB_mmy_known = 1
  else:
    geneB_mmy_known = 0
  if geneB in driver_list:
    geneB_driver = 1
  else:
    geneB_driver = 0
  #druggability
  if fus in depo_dict:
    drug_fusion = 1
  else:
    drug_fusion = 0
  if geneA in depo_dict:
    drug_geneA = 1
  else:
    drug_geneA = 0
  if geneB in depo_dict:
    drug_geneB = 1
  else:
    drug_geneB = 0
  print_list.extend([str(x) for x in [geneA_oncogene, geneA_tsg, geneA_kinase, geneA_mmy_known, geneA_driver, geneB_oncogene, geneB_tsg, geneB_kinase, geneB_mmy_known, geneB_driver, fusion_recurrence[fus], drug_fusion, drug_geneA, drug_geneB]])

  #print it out
  w.write('\t'.join( [ str(x).replace(" ","_") if x else "NA" for x in print_list] )+"\n") #
  
w.close()
